"""
高光谱相机波长标定程序
======================================================================
核心架构：双路并行+融合
  路A（全局多峰拟合）：输入原始光谱，输出每个候选峰的亚像素精确位置
  路B（鲁棒匹配）：输入候选峰粗位置，输出峰与标准波长的对应关系
  融合层：按索引对位，用精确位置替换粗位置，跳过重叠峰
  全局标定：三次多项式拟合，输出最终波长标定系数

图表方案：xlsxwriter 原生图表 + worksheet.insert_textbox 浮动文本框（无边框黑色字体）
异常点判定：D列绝对残差 > 色散值(|趋势线斜率|)
依赖库：numpy, scipy, xlsxwriter
"""
import os
import re
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import xlsxwriter
from datetime import datetime


# ==============================================================
# 配置中心：所有参数统一管理
# ==============================================================
class Config:
    """全局配置类"""
    STD_WL = np.array([
        404.6563, 435.8328, 546.0735,
        594.2000, 626.6495, 650.6528, 659.8953, 667.8276,
        696.5431,
        703.2413, 724.5167, 743.6000,
        763.5106, 772.4207, 794.8176, 811.5311, 826.4522,
        842.4648, 852.1442, 912.2967, 922.4499, 965.7786
    ])

    LAMP_IDX = {
        "HG": [0, 1, 2],
        "NE": [3, 4, 5, 6, 7, 9, 10, 11],
        "AR": [8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
    }

    NE_JUNC_IDX = 9

    PROFILE_WINDOW = 5
    BASE_ITER = 10

    PEAK_H_RATIO = 0.05
    PEAK_P_RATIO = 0.03
    PEAK_MIN_DIST = 3
    OVERLAP_RATIO = 0.3

    FIT_SIGMA_MIN = 0.5
    FIT_SIGMA_MAX = 20.0
    FIT_CENTER_WIN = 8
    OVERLAP_SEP_RATIO = 1.5
    OVERLAP_CONTRIB = 0.10

    RANSAC_ITER = 2000
    RANSAC_TOL = 8.0
    RANSAC_MIN = 3
    SEED = 42

    OUTLIER_SIGMA = 2.5
    OUTLIER_ITER = 3
    MAX_RMSE = 1.5
    MIN_R2 = 0.999

    OUT_DIR = os.path.join(os.path.expanduser("~"), "Desktop")
    HEADER = ["像素位置", "标准波长(nm)", "计算波长(nm)", "绝对残差(nm)"]
    COEFF_NAME = ["a3 (三次项)", "a2 (二次项)", "a1 (一次项)", "a0 (常数项)"]

    def __init__(self):
        self.root = ""
        self.offset = 0
        self.out_name = ""
        self.debug = False


# ==============================================================
# 工具层：ENVI数据读取与预处理
# ==============================================================
class SpecTool:
    _dtype_map = {1: np.uint8, 2: np.int16, 3: np.int32,
                  4: np.float32, 5: np.float64, 12: np.uint16}

    @classmethod
    def read_center(cls, hdr_path, window_size=5):
        if not os.path.exists(hdr_path):
            raise FileNotFoundError(f"HDR文件不存在: {hdr_path}")

        params = cls._parse_hdr(hdr_path)
        sx, sy, bands = params['samples'], params['lines'], params['bands']
        dtype = cls._dtype_map.get(params['data type'], np.float32)
        interleave = params['interleave']
        wl = params['wavelength']

        hdr_dir = os.path.dirname(hdr_path)
        hdr_name = os.path.splitext(os.path.basename(hdr_path))[0]
        raw_path = os.path.join(hdr_dir, hdr_name)
        if not os.path.exists(raw_path):
            raw_path = os.path.join(hdr_dir, hdr_name + '.raw')
            if not os.path.exists(raw_path):
                raise FileNotFoundError(f"RAW数据不存在: {raw_path}")

        data = np.fromfile(raw_path, dtype=dtype)
        sx, sy, bands = cls._auto_correct_dim(data, sx, sy, bands)
        cube = cls._reshape_cube(data, sx, sy, bands, interleave)

        cx, cy = sx // 2, sy // 2
        half = max(1, window_size // 2)
        roi = cube[cx-half:cx+half+1, cy-half:cy+half+1, :]
        spec = np.mean(roi.reshape(-1, bands), axis=0).astype(np.float64)

        if wl is not None and len(wl) == bands and wl[0] > wl[-1]:
            wl, spec = wl[::-1], spec[::-1]

        return spec, wl

    @classmethod
    def remove_baseline(cls, spec, iter=10):
        x = np.arange(len(spec))
        baseline = spec.copy()
        for _ in range(iter):
            coeff = np.polyfit(x, baseline, 2)
            baseline = np.minimum(baseline, np.polyval(coeff, x))
        out = spec - baseline
        return out - np.min(out), baseline

    @classmethod
    def _parse_hdr(cls, hdr_path):
        params = {}
        with open(hdr_path, 'r') as f:
            text = f.read()

        def _get(key):
            m = re.search(rf'^{key}\s*=\s*([^\n;]+)', text, re.IGNORECASE | re.MULTILINE)
            return m.group(1).strip() if m else None

        params['samples'] = int(_get('samples'))
        params['lines'] = int(_get('lines'))
        params['bands'] = int(_get('bands'))
        params['data type'] = int(_get('data type'))
        params['interleave'] = _get('interleave').lower()

        wl_match = re.search(r'^wavelength\s*=\s*\{([^}]+)\}',
                             text, re.IGNORECASE | re.MULTILINE | re.DOTALL)
        params['wavelength'] = np.array([
            float(x) for x in wl_match.group(1).replace('\n', '').replace(' ', '').split(',') if x
        ]) if wl_match else None
        return params

    @classmethod
    def _auto_correct_dim(cls, data, sx, sy, bands):
        total = data.size
        expected = sx * sy * bands
        if total == expected:
            return sx, sy, bands
        if total % (sx * bands) == 0:
            sy = total // (sx * bands)
            print(f"  [警告] HDR lines参数错误，自动校正为 {sy}")
        elif total % (sy * bands) == 0:
            sx = total // (sy * bands)
            print(f"  [警告] HDR samples参数错误，自动校正为 {sx}")
        elif total % (sx * sy) == 0:
            bands = total // (sx * sy)
            print(f"  [警告] HDR bands参数错误，自动校正为 {bands}")
        else:
            raise ValueError(f"数据维度不匹配：期望{expected}，实际{total}")
        return sx, sy, bands

    @classmethod
    def _reshape_cube(cls, data, sx, sy, bands, interleave):
        try_formats = [interleave] + [f for f in ['bil', 'bsq', 'bip'] if f != interleave]
        for fmt in try_formats:
            try:
                d = data.copy()
                if fmt == 'bsq':
                    return d.reshape((bands, sy, sx)).transpose(2, 1, 0)
                elif fmt == 'bil':
                    return d.reshape((sy, bands, sx)).transpose(2, 0, 1)
                elif fmt == 'bip':
                    return d.reshape((sy, sx, bands)).transpose(1, 0, 2)
            except Exception:
                continue
        raise ValueError("数据重组失败，已尝试所有存储格式")


# ==============================================================
# 业务核心：双路并行峰匹配 + 融合
# ==============================================================
class PeakMatcher:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.SEED)
        self.min_disp = 0.8
        self.max_disp = 1.2

    def run(self, spec, lamp, max_first_band=None):
        cfg = self.cfg
        lamp_idx = cfg.LAMP_IDX[lamp]
        std_wl = cfg.STD_WL[lamp_idx]

        n_bands = len(spec)
        nominal_disp = (cfg.STD_WL.max() - cfg.STD_WL.min()) / n_bands
        self.min_disp = nominal_disp * 0.75
        self.max_disp = nominal_disp * 1.25

        candidates = self._detect_peaks(spec)
        if len(candidates) < cfg.RANSAC_MIN:
            print(f"    ⚠ 候选峰不足，无法匹配")
            return {}, None
        coarse_pos = np.array([c["pos"] for c in candidates], dtype=float)
        print(f"    候选峰：{len(candidates)} 个")

        precise_pos, overlap_mask = self._global_fit(spec, candidates)

        best_coeff, best_inliers = self._ransac(coarse_pos, std_wl, max_first_band)
        if best_coeff is None:
            print(f"    ⚠ 未找到有效匹配模型")
            return {}, None
        matches = self._monotonic_match(coarse_pos, std_wl, best_coeff, max_first_band)
        print(f"    匹配完成：{len(matches)} 个峰")

        if len(matches) < 3:
            print(f"    ⚠ 匹配点不足3个，结果不可靠")
            return {}, None

        res = {}
        skip = 0
        for pi, sj in matches:
            if overlap_mask[pi]:
                skip += 1
                continue
            res[lamp_idx[sj]] = float(precise_pos[pi])
        print(f"    融合完成：有效 {len(res)} 个，跳过重叠 {skip} 个")

        return res, best_coeff

    def _detect_peaks(self, spec):
        cfg = self.cfg
        n = len(spec)
        smax = np.max(spec)
        if smax <= 0:
            return []

        peaks, _ = find_peaks(
            spec,
            height=smax * cfg.PEAK_H_RATIO,
            distance=cfg.PEAK_MIN_DIST,
            prominence=smax * cfg.PEAK_P_RATIO
        )
        if len(peaks) == 0:
            return []

        widths = []
        for p in peaks:
            h_half = spec[p] * 0.5
            left = p
            while left > 0 and spec[left] > h_half:
                left -= 1
            right = p
            while right < n - 1 and spec[right] > h_half:
                right += 1
            widths.append(right - left)

        res = []
        for i, p in enumerate(peaks):
            h, w = spec[p], widths[i]
            overlap = False
            if i > 0:
                dl = p - peaks[i-1]
                vl = np.min(spec[peaks[i-1]:p])
                if dl < cfg.PEAK_MIN_DIST * 1.5 and vl > h * cfg.OVERLAP_RATIO:
                    overlap = True
            if i < len(peaks)-1:
                dr = peaks[i+1] - p
                vr = np.min(spec[p:peaks[i+1]])
                if dr < cfg.PEAK_MIN_DIST * 1.5 and vr > h * cfg.OVERLAP_RATIO:
                    overlap = True
            res.append({"pos": int(p), "h": h, "w": w, "overlap": overlap})
        return res

    def _global_fit(self, spec, candidates):
        cfg = self.cfg
        if not candidates:
            return np.array([]), []

        x = np.arange(len(spec), dtype=float)
        y = np.asarray(spec, dtype=float)

        p0 = [float(np.min(y)), 0.0]
        lo, hi = [-np.inf, -np.inf], [np.inf, np.inf]
        for c in candidates:
            p = int(c["pos"])
            h = float(max(c["h"], 1e-6))
            s0 = float(np.clip(c["w"] / 2.355, cfg.FIT_SIGMA_MIN, cfg.FIT_SIGMA_MAX))
            p0.extend([h, float(p), s0])
            lo.extend([0.0, p - cfg.FIT_CENTER_WIN, cfg.FIT_SIGMA_MIN])
            hi.extend([np.inf, p + cfg.FIT_CENTER_WIN, cfg.FIT_SIGMA_MAX])

        try:
            popt, _ = curve_fit(
                self._gauss_model, x, y, p0=p0, bounds=(lo, hi),
                maxfev=5000, method='trf'
            )
        except Exception:
            return np.array([c["pos"] for c in candidates]), [False] * len(candidates)

        n = len(candidates)
        precise_pos = np.array([float(popt[2 + 3*i + 1]) for i in range(n)])
        fit_params = [(float(popt[2+3*i]), float(popt[2+3*i+1]), float(popt[2+3*i+2]))
                      for i in range(n)]

        overlap_mask = self._mark_overlap(fit_params)
        return precise_pos, overlap_mask

    @staticmethod
    def _gauss_model(x, b0, b1, *peak_params):
        y = b0 + b1 * x
        for i in range(0, len(peak_params), 3):
            a, x0, s = peak_params[i], peak_params[i+1], peak_params[i+2]
            y += a * np.exp(-(x - x0)**2 / (2.0 * s * s))
        return y

    def _mark_overlap(self, fit_params):
        cfg = self.cfg
        n = len(fit_params)
        overlap = [False] * n
        for i in range(n):
            ai, x0i, si = fit_params[i]
            for j in range(i+1, n):
                aj, x0j, sj = fit_params[j]
                sep = abs(x0i - x0j) / (si + sj + 1e-9)
                c_ij = aj * np.exp(-(x0i - x0j)**2 / (2*sj*sj)) / (ai + 1e-9)
                c_ji = ai * np.exp(-(x0j - x0i)**2 / (2*si*si)) / (aj + 1e-9)
                if sep < cfg.OVERLAP_SEP_RATIO or c_ij > cfg.OVERLAP_CONTRIB or c_ji > cfg.OVERLAP_CONTRIB:
                    overlap[i] = overlap[j] = True
        return overlap

    def _ransac(self, peak_pos, std_wl, max_first_band=None):
        cfg = self.cfg
        n_p, n_s = len(peak_pos), len(std_wl)
        best_in, best_coeff, best_score = [], None, -1

        for _ in range(cfg.RANSAC_ITER):
            i1, i2 = self.rng.choice(n_p, 2, replace=False)
            j1, j2 = self.rng.choice(n_s, 2, replace=False)
            p1, p2 = peak_pos[i1], peak_pos[i2]
            w1, w2 = std_wl[j1], std_wl[j2]

            if (p2 - p1) * (w2 - w1) <= 0:
                continue

            a = (w2 - w1) / (p2 - p1)
            b = w1 - a * p1

            if not (self.min_disp <= a <= self.max_disp):
                continue

            inliers = []
            last_j = -1
            for i in range(n_p):
                if last_j + 1 >= n_s:
                    break
                if last_j == -1 and max_first_band is not None and peak_pos[i] >= max_first_band:
                    continue
                w_est = a * peak_pos[i] + b
                j_off = np.argmin(np.abs(std_wl[last_j+1:] - w_est))
                j = last_j + 1 + j_off
                if abs(std_wl[j] - w_est) <= cfg.RANSAC_TOL:
                    inliers.append((i, j))
                    last_j = j

            if len(inliers) > best_score and len(inliers) >= cfg.RANSAC_MIN:
                best_score = len(inliers)
                best_coeff = (a, b)
                best_in = inliers
        return best_coeff, best_in

    def _monotonic_match(self, peak_pos, std_wl, coeff, max_first_band=None):
        cfg = self.cfg
        a, b = coeff
        res = []
        last_j = -1
        n_std = len(std_wl)

        for i, p in enumerate(peak_pos):
            if last_j + 1 >= n_std:
                break
            if last_j == -1 and max_first_band is not None and p >= max_first_band:
                continue
            w_est = a * p + b
            j_off = np.argmin(np.abs(std_wl[last_j+1:] - w_est))
            j = last_j + 1 + j_off
            if abs(std_wl[j] - w_est) <= cfg.RANSAC_TOL:
                res.append((i, j))
                last_j = j
        return res


# ==============================================================
# 标定计算：全局三次多项式标定
# ==============================================================
class Calibrator:
    def __init__(self, cfg: Config):
        self.cfg = cfg

    def run(self, bands_std, mask):
        cfg = self.cfg
        n_std = len(cfg.STD_WL)

        bands_std, mask = self._monotonic_check(bands_std, mask)

        valid_idx = np.where(mask)[0]
        if len(valid_idx) < 5:
            raise ValueError(f"有效点仅{len(valid_idx)}个，至少需要5个")

        x, y = bands_std[valid_idx].copy(), cfg.STD_WL[valid_idx].copy()

        init_coeff = np.polyfit(x, y, 1)
        init_res = np.abs(np.polyval(init_coeff, x) - y)
        x, y = x[init_res < 3.0], y[init_res < 3.0]
        if len(x) < 5:
            raise ValueError("预过滤后有效点不足5个")

        inlier = np.ones(len(x), dtype=bool)
        for _ in range(cfg.OUTLIER_ITER):
            if np.sum(inlier) < 5:
                break
            c = np.polyfit(x[inlier], y[inlier], 3)
            res = np.abs(np.polyval(c, x) - y)
            sigma = np.std(res[inlier])
            new_mask = res < max(cfg.OUTLIER_SIGMA * sigma, 0.5)
            if np.array_equal(new_mask, inlier):
                break
            inlier = new_mask

        if np.sum(inlier) < 5:
            raise ValueError("剔除异常后有效点不足5个")

        xf, yf = x[inlier], y[inlier]
        c3 = np.polyfit(xf, yf, 3)

        pred = np.full(n_std, np.nan)
        res = np.full(n_std, np.nan)
        for i in range(n_std):
            if mask[i] and not np.isnan(bands_std[i]):
                pred[i] = np.polyval(c3, bands_std[i])
                res[i] = abs(cfg.STD_WL[i] - pred[i])

        rf = np.abs(np.polyval(c3, xf) - yf)
        rmse = float(np.sqrt(np.mean(rf ** 2)))
        r2 = float(1 - np.sum(rf ** 2) / np.sum((yf - np.mean(yf)) ** 2))
        passed = rmse < cfg.MAX_RMSE and r2 > cfg.MIN_R2

        return {
            "coeff": c3, "bands": bands_std, "pred": pred, "res": res,
            "rmse": rmse, "r2": r2, "passed": passed,
            "n_valid": int(np.sum(inlier)), "n_out": int(len(x) - np.sum(inlier))
        }

    def _monotonic_check(self, bands_std, mask):
        std_wl = self.cfg.STD_WL
        valid_idx = np.where(mask)[0]
        if len(valid_idx) < 3:
            return bands_std, mask

        sorted_idx = valid_idx[np.argsort(std_wl[valid_idx])]
        sorted_bands = bands_std[sorted_idx]

        i = 0
        while i < len(sorted_bands) - 1:
            if sorted_bands[i] >= sorted_bands[i + 1]:
                wl_i, wl_j = std_wl[sorted_idx[i]], std_wl[sorted_idx[i+1]]
                a_local = (wl_j - wl_i) / (sorted_bands[i+1] - sorted_bands[i]) if sorted_bands[i+1] != sorted_bands[i] else 0
                b_local = wl_i - a_local * sorted_bands[i]
                res_i = abs(wl_i - (a_local * sorted_bands[i] + b_local))
                res_j = abs(wl_j - (a_local * sorted_bands[i+1] + b_local))

                bad_idx = sorted_idx[i] if res_i > res_j else sorted_idx[i+1]
                mask[bad_idx] = False
                bands_std[bad_idx] = np.nan
                print(f"  [校验] 剔除逆序错配点：{std_wl[bad_idx]:.4f}nm")

                valid_idx = np.where(mask)[0]
                sorted_idx = valid_idx[np.argsort(std_wl[valid_idx])]
                sorted_bands = bands_std[sorted_idx]
                i = 0
            else:
                i += 1
        return bands_std, mask


# ==============================================================
# 输出层：xlsxwriter 报表生成
# ==============================================================
class ExcelWriter:
    """
    Excel结果输出器（xlsxwriter）
    - 每个数据列使用专用 num_format，保证显示位数固定
    - 文本框：无边框、无填充、黑色字体
    - 异常点判定：D列残差 > 色散值(|趋势线斜率|) → 整行红字
    """

    _LAMP_BG = {"HG": "#FFFF00", "NE": "#00CCFF", "AR": "#FF9900"}
    _GRAY_BG = "#F2F2F2"

    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.wb = None

    # ---------- 生命周期 ----------
    def start(self, path):
        self.wb = xlsxwriter.Workbook(path, {'nan_inf_to_errors': True})

    def close(self):
        if self.wb is not None:
            self.wb.close()
            self.wb = None

    # ---------- 格式工厂 ----------
    def _fmt_base(self, **extra):
        d = {'align': 'center', 'valign': 'vcenter',
             'border': 1, 'border_color': '#BFBFBF'}
        d.update(extra)
        return self.wb.add_format(d)

    def _make_formats(self):
        """
        一次性创建本 sheet 所需的全部格式。
        关键：每个含数字的格式都带 num_format，Excel 才会按固定位数显示。
        """
        f = {}
        # 表头 / 灰底 / 系数
        f['header'] = self._fmt_base(bold=True, bg_color='#D6EAF8')
        f['gray'] = self._fmt_base(bg_color=self._GRAY_BG)
        f['coeff_name'] = self._fmt_base(bg_color=self._GRAY_BG)
        # 系数值：科学计数法（覆盖 a3~a0 跨度大的场景）
        f['coeff_val'] = self._fmt_base(bg_color=self._GRAY_BG,
                                        num_format='0.00E+00')

        # ---- 4 列数据：黑字版 / 红字版，各自带 num_format ----
        # A列 像素位置：2 位小数
        f['pixel']     = self._fmt_base(num_format='0.00')
        f['pixel_red'] = self._fmt_base(num_format='0.00', font_color='#FF0000')
        # C列 计算波长：2 位小数
        f['pred']      = self._fmt_base(num_format='0.00')
        f['pred_red']  = self._fmt_base(num_format='0.00', font_color='#FF0000')
        # D列 绝对残差：3 位小数
        f['res']       = self._fmt_base(num_format='0.000')
        f['res_red']   = self._fmt_base(num_format='0.000', font_color='#FF0000')

        # B列：灯种底色 + 4 位小数（黑字 / 红字两版）
        for lamp, bg in self._LAMP_BG.items():
            f[f'b_{lamp}']     = self._fmt_base(bg_color=bg, num_format='0.0000')
            f[f'b_{lamp}_red'] = self._fmt_base(bg_color=bg, num_format='0.0000',
                                                 font_color='#FF0000')
        return f

    # ---------- 主流程 ----------
    def write_sheet(self, res, bm, offset):
        cfg = self.cfg
        n_std = len(cfg.STD_WL)
        sheet_name = f"offset{offset}-{bm}"
        ws = self.wb.add_worksheet(sheet_name)
        fmt = self._make_formats()

        # 列宽
        ws.set_column(0, 0, 14)
        ws.set_column(1, 1, 16)
        ws.set_column(2, 2, 16)
        ws.set_column(3, 3, 14)
        ws.set_column(4, 4, 16)
        ws.set_column(5, 5, 22)

        # 索引→灯种
        idx_to_lamp = {}
        for lamp, indices in cfg.LAMP_IDX.items():
            for idx in indices:
                idx_to_lamp[idx] = lamp

        # 1. 表头
        for col, name in enumerate(cfg.HEADER):
            ws.write(0, col, name, fmt['header'])

        # 2. 有效点线性拟合 → 色散值
        valid_idx = [i for i in range(n_std) if not np.isnan(res["bands"][i])]
        c1 = None
        disp = None
        if len(valid_idx) >= 2:
            xv = np.array([res["bands"][i] for i in valid_idx], dtype=float)
            yv = np.array([cfg.STD_WL[i] for i in valid_idx], dtype=float)
            c1 = np.polyfit(xv, yv, 1)
            disp = abs(c1[0])

        # 3. 主数据区
        for i in range(n_std):
            row = i + 1
            has_band = not np.isnan(res["bands"][i])
            has_pred = not np.isnan(res["pred"][i])
            is_out = (has_pred and disp is not None and res["res"][i] > disp)

            # A 列 像素位置：数值本身已是 float，仅靠 num_format 显示 2 位
            if not has_band:
                ws.write(row, 0, "不可用", fmt['gray'])
            else:
                ws.write_number(row, 0, float(res["bands"][i] + 1),
                                fmt['pixel_red'] if is_out else fmt['pixel'])

            # B 列 标准波长 + 灯种底色
            lamp = idx_to_lamp.get(i, "AR")
            wl_val = float(cfg.STD_WL[i])
            if not has_band:
                ws.write_number(row, 1, wl_val, fmt['gray'])
            else:
                key = f'b_{lamp}_red' if is_out else f'b_{lamp}'
                ws.write_number(row, 1, wl_val, fmt[key])

            # C 列 计算波长 / D 列 绝对残差
            if has_pred:
                ws.write_number(row, 2, float(res["pred"][i]),
                                fmt['pred_red'] if is_out else fmt['pred'])
                ws.write_number(row, 3, float(res["res"][i]),
                                fmt['res_red'] if is_out else fmt['res'])

        # 4. 多项式系数
        for i, (name, val) in enumerate(zip(cfg.COEFF_NAME, res["coeff"])):
            ws.write(i, 4, name, fmt['coeff_name'])
            ws.write_number(i, 5, float(val), fmt['coeff_val'])

        # 5. 图表 + 文本框
        if len(valid_idx) >= 5 and c1 is not None:
            self._add_chart(ws, sheet_name, res, valid_idx, bm, offset, n_std, c1)

    # ---------- 图表 + 文本框 ----------
    def _add_chart(self, ws, sheet_name, res, valid_idx, bm, offset, n_std, c1):
        cfg = self.cfg

        data_sheet = f"_chart_data_{bm}"
        ws_data = self.wb.add_worksheet(data_sheet)
        ws_data.hide()

        n_pt = len(valid_idx)
        x_arr = [float(res["bands"][i] + 1) for i in valid_idx]
        y_arr = [float(cfg.STD_WL[i]) for i in valid_idx]
        for k, (xa, ya) in enumerate(zip(x_arr, y_arr)):
            ws_data.write_number(k + 1, 0, xa)
            ws_data.write_number(k + 1, 1, ya)

        x_fit = np.linspace(min(x_arr), max(x_arr), 100)
        y_fit = np.polyval(c1, x_fit)
        for k in range(len(x_fit)):
            ws_data.write_number(k + 1, 3, float(x_fit[k]))
            ws_data.write_number(k + 1, 4, float(y_fit[k]))

        # 创建散点图
        chart = self.wb.add_chart({'type': 'scatter'})

        # 系列1：标定点（圆点 + 黑色细实线）
        chart.add_series({
            'categories': [data_sheet, 1, 0, n_pt, 0],
            'values':     [data_sheet, 1, 1, n_pt, 1],
            'marker': {'type': 'circle', 'size': 5,
                       'border': {'color': 'black'},
                       'fill': {'color': 'black'}},
            'line': {'color': 'black', 'width': 1.0},
            'name': '',
        })

        # 系列2：趋势线（蓝色虚线，无标记）
        chart.add_series({
            'categories': [data_sheet, 1, 3, 100, 3],
            'values':     [data_sheet, 1, 4, 100, 4],
            'marker': {'type': 'none'},
            'line': {'color': '#0070C0', 'width': 1.5, 'dash_type': 'dash'},
            'name': '',
        })

        chart.set_title({'name': f'offset{offset}-{bm}波长标定曲线'})
        chart.set_x_axis({'name': 'Pixel Location', 'min': 0})
        chart.set_y_axis({'name': 'Wavelength(nm)'})
        chart.set_size({'width': 800, 'height': 500})
        chart.set_legend({'none': True})

        chart_top_row_1based = n_std + 3
        ws.insert_chart('G2', chart,
                        {'x_offset': 0, 'y_offset': 0})

        # ★ 文本框：无边框、无填充、黑色字体
        formula_text = f"y = {c1[0]:.4f}x + {c1[1]:.2f}"
        tb_row_0based = 8
        tb_col_0based = 10
        ws.insert_textbox(tb_row_0based, tb_col_0based, formula_text, {
            'width': 200, 'height': 32,
            'x_offset': 18, 'y_offset': 10,
            'font': {'name': 'Consolas', 'size': 11, 'color': '#000000'},
            'align': {'horizontal': 'center', 'vertical': 'center'},
            'line': {'none': True},   # 彻底隐藏边框线
            'fill': {'none': True},   # 彻底隐藏背景填充
        })


# ==============================================================
# 主控层：流程编排
# ==============================================================
class App:
    def __init__(self):
        self.cfg = Config()
        self.spec_tool = SpecTool()
        self.matcher = PeakMatcher(self.cfg)
        self.cal = Calibrator(self.cfg)
        self.writer = ExcelWriter(self.cfg)
        self.warn = []

    def run(self):
        self._input()
        subs, bins = self._scan_folders()
        print(f"识别到 {len(subs)} 个子文件夹，bin模式：{', '.join(bins)}\n")

        n_std = len(self.cfg.STD_WL)
        NE_703_IDX = self.cfg.NE_JUNC_IDX

        out_file = os.path.join(
            self.cfg.OUT_DIR,
            f"{self.cfg.out_name}.xlsx"
        )

        self.writer.start(out_file)

        try:
            for bm in bins:
                print(f"{'='*22} 处理 {bm} {'='*22}")
                bands_std = np.full(n_std, np.nan)
                mask = np.zeros(n_std, dtype=bool)

                ne_703_band = None
                for sf in subs:
                    if sf["bin"] == bm and sf["lamp"] == "NE":
                        self._process_lamp(sf, bands_std, mask)
                        if mask[NE_703_IDX]:
                            ne_703_band = bands_std[NE_703_IDX]
                            print(f"  [交界参考] NE灯703.24nm像素：{ne_703_band+1:.1f}")

                for sf in subs:
                    if sf["bin"] == bm and sf["lamp"] == "AR":
                        self._process_lamp(sf, bands_std, mask,
                                           max_first_band=ne_703_band)

                for sf in subs:
                    if sf["bin"] == bm and sf["lamp"] == "HG":
                        self._process_lamp(sf, bands_std, mask)

                n_valid = int(np.sum(mask))
                print(f"\n---- {bm} 全局标定 ----")
                print(f"有效匹配点：{n_valid} / {n_std}")

                if n_valid < 5:
                    self.warn.append(f"{bm} 有效点不足5个")
                    print("× 有效点不足，跳过标定\n")
                    continue

                try:
                    cal_res = self.cal.run(bands_std, mask)
                except Exception as e:
                    self.warn.append(f"{bm} 标定失败：{e}")
                    print(f"× 标定失败：{e}\n")
                    continue

                status = "✅ 合格" if cal_res["passed"] else "❌ 不合格"
                print(f"RMSE：{cal_res['rmse']:.4f} nm")
                print(f"R²：{cal_res['r2']:.6f}")
                print(f"质量判定：{status}")

                self.writer.write_sheet(cal_res, bm, self.cfg.offset)
                print(f"✅ {bm} 工作表已添加\n")
        finally:
            self.writer.close()

        print(f"✅ 所有结果已输出：{out_file}")

        if self.warn:
            print("\n" + "-"*45)
            print("⚠ 告警汇总：")
            for w in self.warn:
                print(f"  · {w}")

        input("\n按回车键退出")

    def _input(self):
        print("=" * 55)
        print("高光谱波长标定程序")
        print("=" * 55)
        while True:
            p = input("请输入数据根文件夹路径：").strip().strip('"').strip("'")
            if os.path.isdir(p):
                self.cfg.root = p
                break
            print("× 路径不存在，请重新输入")
        while True:
            try:
                self.cfg.offset = int(input("请输入offset整数值：").strip())
                break
            except ValueError:
                print("× 请输入整数")
        # 输出文件名
        default_name = f"{datetime.now():%Y%m%d}_波长标定结果"#默认文件名
        name = input(f"请输入输出文件名（不含扩展名，直接回车默认 {default_name}）：").strip()
        # 去掉用户可能误输入的后缀、路径分隔符
        name = name.replace(".xlsx", "").replace("\\", "").replace("/", "")
        self.cfg.out_name = name if name else default_name


        print("\n开始处理...\n")

    def _scan_folders(self):
        root = self.cfg.root
        subs, bins = [], set()
        for name in os.listdir(root):
            fp = os.path.join(root, name)
            if not os.path.isdir(fp):
                continue
            lamp = None
            for L in ["HG", "NE", "AR"]:
                if L in name.upper():
                    lamp = L
                    break
            if not lamp:
                continue
            m = re.search(r"\d+bin", name, re.IGNORECASE)
            if not m:
                continue
            bm = m.group().lower()
            subs.append({"name": name, "path": fp, "lamp": lamp, "bin": bm})
            bins.add(bm)
        if not bins:
            raise ValueError("未识别到有效子文件夹，请检查命名格式")
        return subs, sorted(bins)

    def _process_lamp(self, sf, bands_std, mask, max_first_band=None):
        print(f"\n▶ {sf['lamp']} 灯段：{sf['name']}")
        hdr_list = [f for f in os.listdir(sf["path"]) if f.lower().endswith(".hdr")]
        if not hdr_list:
            self.warn.append(f"{sf['name']} 无HDR文件")
            print("  × 无HDR文件，跳过")
            return

        try:
            hdr_path = os.path.join(sf["path"], hdr_list[0])
            spec, _ = self.spec_tool.read_center(hdr_path, self.cfg.PROFILE_WINDOW)
            print(f"  光谱尺寸：{len(spec)} 波段")
            print(f"  基线校正...", end="")
            spec, _ = self.spec_tool.remove_baseline(spec, self.cfg.BASE_ITER)
            print("完成")
        except Exception as e:
            self.warn.append(f"{sf['name']} 读取失败：{e}")
            print(f"  × 读取失败：{e}")
            return

        try:
            match_res, _ = self.matcher.run(spec, sf["lamp"], max_first_band)
            for k, v in match_res.items():
                bands_std[k] = v
                mask[k] = True
        except Exception as e:
            self.warn.append(f"{sf['lamp']}匹配失败：{e}")
            print(f"  × 匹配失败：{e}")


if __name__ == "__main__":
    App().run()