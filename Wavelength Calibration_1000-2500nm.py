# -*- coding: utf-8 -*-
"""
高光谱相机波长标定程序（1000-2500nm）
======================================================================
核心流程：
  1. 读 HDR/RAW，取中心区均值光谱
  2. 多项式迭代基线校正
  3. find_peaks 检测候选峰
  4. 全局多峰高斯拟合精化峰位 + 标记重叠峰
  5. RANSAC 求初始色散模型 → 双指针单调贪心匹配
  6. 对未匹配标准线做局部二次拟合补全（重叠区跳过）
  7. 合并全部灯段，迭代 sigma 剔除 → 线性拟合
  8. 输出 Excel（数据表 + 标定曲线 + 浮动公式文本框）

重叠峰判据（本轮修正）：
  距离近 AND 光强相当，二者同时满足才判为重叠。
  一强一弱即使相距很近，弱峰仍可独立辨位，不标重叠。

依赖：numpy, scipy, xlsxwriter
"""

import argparse
import os
import re
import sys
from datetime import datetime

import numpy as np
import xlsxwriter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths


# ==============================================================
# 配置中心：所有可调参数集中管理
# ==============================================================
class Config:
    # ---- 标准波长表（按灯种分组，LAMP_IDX 给出各灯对应的下标）----
    STD_WL = np.array([
        1046.7177, 1148.8109, 1152.2745, 1176.67924, 1181.9377,
        1211.2326, 1248.7663, 1270.2281, 1295.6659, 1350.4191,
        1363.422, 1371.8577, 1409.364, 1442.6793, 1473.4436,
        1504.4191, 1523.9615, 1689.0441, 1694.058, 1816.7315,
        1827.6642, 1859.7698, 2061.623, 2104.127, 2170.811,
        2190.2513, 2253.038, 2337.296, 2363.648, 2395.14, 2436.501,
    ])
    LAMP_IDX = {
        "AR": [0, 1, 5, 6, 7, 8, 9, 11, 12, 15, 18, 22],
        "NE": [2, 3, 20, 21, 23, 24, 26, 27, 28, 29, 30],
        "KR": [4, 10, 13, 14, 16, 17, 19, 25],
    }

    # ---- 光谱读取与基线校正 ----
    PROFILE_WINDOW = 5          # 中心区采样窗口边长
    BASE_ITER = 10              # 基线迭代次数

    # ---- 候选峰检测 ----
    PEAK_H_RATIO = 0.015        # 绝对高度阈值（相对全谱最大值）
    PEAK_P_RATIO = 0.010        # 峰突出度阈值（相对全谱最大值）
    PEAK_MIN_DIST = 2           # 相邻峰最小间隔（像素）
    OVERLAP_RATIO = 0.30        # 粗糙重叠判定用（检测阶段）

    # ---- 全局多峰高斯拟合 ----
    FIT_SIGMA_MIN = 0.5
    FIT_SIGMA_MAX = 12.0
    FIT_CENTER_WIN = 2.0        # 拟合中心允许偏离粗峰位的最大像素

    # ---- 重叠峰判据（距离近 AND 光强相当）----
    OVERLAP_SEP_RATIO = 1.5         # 分离度阈值 |Δx|/(σi+σj)
    OVERLAP_INTENSITY_RATIO = 0.30  # 光强比阈值 min(ai,aj)/max(ai,aj)

    # ---- RANSAC 匹配 ----
    RANSAC_ITER = 5000
    RANSAC_TOL = 2.5            # 初始匹配容差上限（nm）
    RANSAC_MIN = 3
    DISP_TOL = 0.35             # 色散允许偏离标称值比例
    SPACE_CHECK_TOL = 0.35      # 相邻对局部色散一致性阈值
    SEED = 42

    # ---- 局部二次拟合补全 ----
    LOCAL_FIT_WINDOW = 15           # 局部窗口半宽（像素）
    LOCAL_FIT_CENTER_WIN = 3.0      # 拟合中心距预期像素的最大距离
    LOCAL_FIT_OVERLAP_SKIP = 8.0    # 预期位置多少像素内有重叠峰则跳过
    LOCAL_FIT_MIN_SNR = 3.0         # 拟合峰高 / 残差 RMS 下限
    LOCAL_FIT_MIN_REL = 0.01        # 拟合峰高 / 全谱最大值下限

    # ---- 异常剔除与质量判据 ----
    OUTLIER_SIGMA = 2.5
    OUTLIER_ITER = 3
    MAX_RMSE = 3.0
    MIN_R2 = 0.999

    OUT_DIR = os.environ.get(
        "WAVECAL_OUT_DIR", os.path.join(os.path.expanduser("~"), "Desktop"))
    HEADER = ["像素位置", "标准波长(nm)", "计算波长(nm)", "绝对残差(nm)"]
    COEFF_NAME = ["a1 (一次项)", "a0 (常数项)"]

    def __init__(self):
        self.root = ""
        self.offset = 0
        self.out_name = ""


# ==============================================================
# 工具层：ENVI 数据读取 + 基线校正
# ==============================================================
class SpecTool:
    _dtype_map = {
        1: np.uint8, 2: np.int16, 3: np.int32,
        4: np.float32, 5: np.float64, 12: np.uint16
    }
    _VALID_INTERLEAVE = ("bil", "bsq", "bip")

    @classmethod
    def read_center(cls, hdr_path, window_size=5):
        """读 HDR/RAW，取中心 window_size×window_size 区域的平均光谱"""
        if not os.path.exists(hdr_path):
            raise FileNotFoundError(f"HDR不存在: {hdr_path}")
        params = cls._parse_hdr(hdr_path)
        sx, sy, bands = params['samples'], params['lines'], params['bands']
        dtype = cls._dtype_map.get(params['data type'], np.float32)

        hdr_dir = os.path.dirname(hdr_path)
        hdr_name = os.path.splitext(os.path.basename(hdr_path))[0]
        raw_path = None
        for cand in (hdr_name, hdr_name + '.raw'):
            p = os.path.join(hdr_dir, cand)
            if os.path.exists(p):
                raw_path = p
                break
        if raw_path is None:
            raise FileNotFoundError(f"RAW不存在: {os.path.join(hdr_dir, hdr_name)}")

        data = np.fromfile(raw_path, dtype=dtype)
        total = data.size
        expected = sx * sy * bands

        # HDR 声明的维度与实际文件大小不符时自动修正
        if total != expected:
            if total % (sx * bands) == 0:
                sy = total // (sx * bands)
            elif total % (sy * bands) == 0:
                sx = total // (sy * bands)
            elif total % (sx * sy) == 0:
                bands = total // (sx * sy)
            else:
                raise ValueError(f"维度不匹配：期望{expected}，实际{total}")

        cube = cls._reshape_cube(data, sx, sy, bands, params['interleave'])

        cx, cy = sx // 2, sy // 2
        half = max(1, window_size // 2)
        x0, x1 = max(0, cx - half), min(sx, cx + half + 1)
        y0, y1 = max(0, cy - half), min(sy, cy + half + 1)
        spec = cube[x0:x1, y0:y1, :].mean(axis=(0, 1)).astype(np.float64)

        # 波长轴反向时同步反转光谱，保证像素序与波长序单调递增
        wl = params['wavelength']
        if wl is not None and len(wl) == bands and wl[0] > wl[-1]:
            wl = wl[::-1]
            spec = spec[::-1]
        return spec, wl

    @classmethod
    def remove_baseline(cls, spec, iter=10):
        """多项式迭代下包络基线校正：每次迭代用 2 次多项式包络取低值"""
        x = np.arange(len(spec))
        baseline = spec.copy()
        for _ in range(iter):
            coeff = np.polyfit(x, baseline, 2)
            baseline = np.minimum(baseline, np.polyval(coeff, x))
        out = spec - baseline
        return out - np.min(out), baseline

    @classmethod
    def _parse_hdr(cls, hdr_path):
        with open(hdr_path, 'r') as f:
            text = f.read()

        def _get_val(key, default=None):
            m = re.search(rf'^{key}\s*=\s*([^\n;]+)', text,
                          re.IGNORECASE | re.MULTILINE)
            return m.group(1).strip() if m else default

        interleave = (_get_val('interleave', 'bsq')).lower()
        if interleave not in cls._VALID_INTERLEAVE:
            raise ValueError(f"不支持的存储格式: {interleave}")

        wl_match = re.search(r'^wavelength\s*=\s*\{([^}]+)\}', text,
                             re.IGNORECASE | re.MULTILINE | re.DOTALL)
        if wl_match:
            wl_str = wl_match.group(1).replace('\n', '').replace(' ', '')
            wavelength = np.array([float(v) for v in wl_str.split(',') if v])
        else:
            wavelength = None

        return {
            'samples': int(_get_val('samples')),
            'lines': int(_get_val('lines')),
            'bands': int(_get_val('bands')),
            'data type': int(_get_val('data type')),
            'interleave': interleave,
            'wavelength': wavelength,
        }

    @classmethod
    def _reshape_cube(cls, data, sx, sy, bands, interleave):
        if interleave == 'bsq':
            return data.reshape((bands, sy, sx)).transpose(2, 1, 0)
        if interleave == 'bil':
            return data.reshape((sy, bands, sx)).transpose(2, 0, 1)
        if interleave == 'bip':
            return data.reshape((sy, sx, bands)).transpose(1, 0, 2)
        raise ValueError(f"不支持的存储格式: {interleave}")


# ==============================================================
# 峰匹配：候选检测 → 全局拟合 → RANSAC → 单调匹配 → 局部补全
# ==============================================================
class PeakMatcher:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.SEED)
        self.min_disp = 0.8
        self.max_disp = 1.2

    def run(self, spec, lamp):
        """
        单灯段完整匹配流程。
        返回 (normal_res, overlap_res)：
          normal_res  : {全局标准波长下标: 亚像素位置}  正常匹配 + 局部补全
          overlap_res : {全局标准波长下标: 亚像素位置}  匹配上但被判为重叠的
        """
        cfg = self.cfg
        lamp_idx = cfg.LAMP_IDX[lamp]
        std_wl = cfg.STD_WL[lamp_idx]
        n_bands = len(spec)

        # 色散 = 波长跨度 / 波段数，用于约束 RANSAC 的斜率范围
        wl_range = cfg.STD_WL.max() - cfg.STD_WL.min()
        nominal_disp = wl_range / n_bands
        self.min_disp = nominal_disp * (1 - cfg.DISP_TOL)
        self.max_disp = nominal_disp * (1 + cfg.DISP_TOL)

        # 匹配容差：不超过 RANSAC_TOL，且不超过最小标准线间距的 60%
        tol = self._match_tol(std_wl, cfg)

        # ---- 阶段1：候选峰检测 ----
        candidates = self._detect_peaks(spec)
        if len(candidates) < cfg.RANSAC_MIN:
            print("    候选峰不足，跳过")
            return {}, {}
        coarse_pos = np.array([c["pos"] for c in candidates], dtype=float)
        print(f"    候选峰：{len(candidates)} 个")

        # ---- 阶段2：全局多峰高斯拟合 ----
        precise_pos, overlap_mask = self._global_fit(spec, candidates)

        # ---- 阶段3：RANSAC 求初始色散模型 ----
        best_coeff, _ = self._ransac(coarse_pos, std_wl, tol)
        if best_coeff is None:
            print("    未找到有效匹配模型")
            return {}, {}
        print(f"    RANSAC：a1={best_coeff[0]:.5f}  a0={best_coeff[1]:.2f}")

        # ---- 阶段4：单调贪心匹配 ----
        matches = self._greedy_match(coarse_pos, std_wl,
                                     best_coeff[0], best_coeff[1], tol)
        print(f"    匹配：{len(matches)} 个峰")

        # ---- 阶段5：未匹配标准线诊断 ----
        matched_sj = set(sj for _, sj in matches)
        unmatched = [sj for sj in range(len(std_wl)) if sj not in matched_sj]
        if unmatched:
            a, b = best_coeff
            print(f"    未匹配标准线 {len(unmatched)} 个：")
            for sj in unmatched:
                wl_std = float(std_wl[sj])
                px_expect = (wl_std - b) / a if abs(a) > 1e-9 else -1
                near = np.abs(coarse_pos - px_expect) <= 15
                if np.any(near):
                    idxs = np.where(near)[0]
                    info = "  附近候选峰: " + ", ".join(
                        f"#{int(i)}@{int(coarse_pos[i])}"
                        f"(ov={overlap_mask[i]})" for i in idxs)
                else:
                    info = "  ✗ 附近无候选峰"
                print(f"      · {wl_std:.4f}nm 预期px≈{px_expect:.2f}{info}")

        # ---- 阶段6：局部二次拟合补全未匹配标准线 ----
        extra = self._local_refit_unmatched(
            spec, coarse_pos, overlap_mask, matches, std_wl, best_coeff,
            lamp_idx, unmatched)
        if extra:
            print(f"    局部补全：{len(extra)} 个")

        # ---- 汇总：按索引对位，重叠峰单独归类 ----
        normal_res, overlap_res = {}, {}
        for pi, sj in matches:
            std_idx = lamp_idx[sj]
            if overlap_mask[pi]:
                overlap_res[std_idx] = float(precise_pos[pi])
            else:
                normal_res[std_idx] = float(precise_pos[pi])
        normal_res.update(extra)  # 局部补全并入正常匹配

        print(f"    融合：正常 {len(normal_res)} 个，重叠 {len(overlap_res)} 个")
        return normal_res, overlap_res

    # ----------------------------------------------------------
    # 候选峰检测：find_peaks + 半高宽 + 相邻峰的粗糙重叠标记
    # ----------------------------------------------------------
    def _detect_peaks(self, spec):
        cfg = self.cfg
        smax = float(np.max(spec))
        if smax <= 0:
            return []

        peaks, _ = find_peaks(
            spec,
            height=smax * cfg.PEAK_H_RATIO,
            distance=cfg.PEAK_MIN_DIST,
            prominence=smax * cfg.PEAK_P_RATIO
        )
        if peaks.size == 0:
            return []

        # 向量化半高宽
        widths = peak_widths(spec, peaks, rel_height=0.5)[0]

        res = []
        n = len(peaks)
        for i in range(n):
            p = int(peaks[i])
            h = float(spec[p])
            overlap = False
            # 检查相邻峰：间隔小且山谷高 → 粗糙重叠
            for nb in (i - 1, i + 1):
                if 0 <= nb < n:
                    q = int(peaks[nb])
                    if abs(p - q) < cfg.PEAK_MIN_DIST * 1.5:
                        lo, hi = (p, q) if p < q else (q, p)
                        v = float(np.min(spec[lo:hi + 1]))
                        if v > h * cfg.OVERLAP_RATIO:
                            overlap = True
            res.append({"pos": p, "h": h,
                        "w": float(widths[i]), "overlap": overlap})
        return res

    # ----------------------------------------------------------
    # 全局多峰高斯拟合：一次拟合全部候选峰
    # 参数向量 = [b0, b1] + 每峰 [A, x0, sigma]
    # ----------------------------------------------------------
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
            s0 = float(np.clip(c["w"] / 2.355,
                               cfg.FIT_SIGMA_MIN, cfg.FIT_SIGMA_MAX))
            p0.extend([h, float(p), s0])
            lo.extend([0.0, p - cfg.FIT_CENTER_WIN, cfg.FIT_SIGMA_MIN])
            hi.extend([np.inf, p + cfg.FIT_CENTER_WIN, cfg.FIT_SIGMA_MAX])

        try:
            popt, _ = curve_fit(
                self._gauss_model, x, y, p0=p0, bounds=(lo, hi),
                max_nfev=20000, method='trf'
            )
        except Exception:
            # 拟合失败时退回粗峰位，避免整个流程中断
            precise_pos = np.array([c["pos"] for c in candidates], dtype=float)
            overlap_mask = [c["overlap"] for c in candidates]
            return precise_pos, overlap_mask

        n = len(candidates)
        precise_pos = np.array([float(popt[2 + 3 * i + 1]) for i in range(n)])
        fit_params = [(float(popt[2 + 3 * i]), float(popt[2 + 3 * i + 1]),
                       float(popt[2 + 3 * i + 2])) for i in range(n)]
        overlap_mask = self._mark_overlap(fit_params)
        return precise_pos, overlap_mask

    @staticmethod
    def _gauss_model(x, b0, b1, *peak_params):
        """多峰高斯 + 线性基线的合成模型"""
        y = b0 + b1 * x
        for i in range(0, len(peak_params), 3):
            a, x0, s = peak_params[i], peak_params[i + 1], peak_params[i + 2]
            y += a * np.exp(-(x - x0) ** 2 / (2.0 * s * s))
        return y

    def _mark_overlap(self, fit_params):
        """
        重叠峰判据（核心）：
          距离近 AND 光强相当，二者同时满足才判重叠。
          1) sep = |Δx| / (σi + σj) < OVERLAP_SEP_RATIO
          2) ratio = min(ai,aj) / max(ai,aj) > OVERLAP_INTENSITY_RATIO
        只检查按中心排序后的相邻峰对（O(n log n)）。
        一强一弱即使距离很近，弱峰仍可独立辨位 → 不标重叠。
        """
        cfg = self.cfg
        n = len(fit_params)
        overlap = [False] * n
        if n < 2:
            return overlap

        order = np.argsort([p[1] for p in fit_params])
        for k in range(n - 1):
            i = int(order[k])
            j = int(order[k + 1])
            ai, x0i, si = fit_params[i]
            aj, x0j, sj = fit_params[j]

            sep = abs(x0i - x0j) / (si + sj + 1e-9)
            if sep >= cfg.OVERLAP_SEP_RATIO:
                continue

            ratio = min(ai, aj) / max(ai, aj, 1e-9)
            if ratio < cfg.OVERLAP_INTENSITY_RATIO:
                continue

            overlap[i] = overlap[j] = True
        return overlap

    # ----------------------------------------------------------
    # 二次局部拟合：对未匹配标准线在预期像素位置附近做单峰拟合
    # ----------------------------------------------------------
    def _local_refit_unmatched(self, spec, coarse_pos, overlap_mask, matches,
                               std_wl, coeff, lamp_idx, unmatched):
        cfg = self.cfg
        a, b = coeff
        if abs(a) < 1e-9:
            return {}

        overlap_arr = np.array(overlap_mask, dtype=bool)
        smax = float(np.max(spec))

        result = {}
        for sj in unmatched:
            wl_std = float(std_wl[sj])
            px_expect = (wl_std - b) / a

            # 预期位置附近存在重叠峰 → 跳过（避免污染）
            near = np.abs(coarse_pos - px_expect) <= cfg.LOCAL_FIT_OVERLAP_SKIP
            ov_hit = near & overlap_arr
            if np.any(ov_hit):
                bad = int(np.where(ov_hit)[0][0])
                print(f"      [跳过] {wl_std:.4f}nm 预期px={px_expect:.2f}，"
                      f"邻域存在重叠峰 #{bad}@{int(coarse_pos[bad])}")
                continue

            center, reason = self._local_gauss_fit(spec, px_expect, smax)
            if center is None:
                print(f"      [跳过] {wl_std:.4f}nm 预期px={px_expect:.2f}，"
                      f"局部拟合失败（{reason}）")
                continue

            print(f"      [补全] {wl_std:.4f}nm 预期px={px_expect:.2f}"
                  f" → 拟合px={center:.3f}")
            result[lamp_idx[sj]] = float(center)
        return result

    def _local_gauss_fit(self, spec, px_expect, smax):
        """
        局部窗口内拟合 A·exp(-(x-xc)²/(2σ²)) + b0 + b1·(x-xc)。
        返回 (亚像素中心, 失败原因)；成功时原因为 None。
        """
        cfg = self.cfg
        n = len(spec)
        w = int(cfg.LOCAL_FIT_WINDOW)
        lo = max(0, int(round(px_expect)) - w)
        hi = min(n, int(round(px_expect)) + w + 1)
        if hi - lo < 5:
            return None, "窗口太窄"

        x = np.arange(lo, hi, dtype=float)
        y = spec[lo:hi].astype(float)

        # 用窗口内最大值作为初始峰位置与峰高
        i0 = int(np.argmax(y))
        p0 = [float(y[i0] - np.min(y)), float(x[i0]), 2.0,
              float(np.min(y)), 0.0]

        def model(xx, A, xc, s, b0, b1):
            return (A * np.exp(-(xx - xc) ** 2 / (2 * s * s))
                    + b0 + b1 * (xx - xc))

        try:
            popt, _ = curve_fit(
                model, x, y, p0=p0,
                bounds=(
                    [0.0, px_expect - cfg.LOCAL_FIT_CENTER_WIN,
                     cfg.FIT_SIGMA_MIN, -np.inf, -np.inf],
                    [np.inf, px_expect + cfg.LOCAL_FIT_CENTER_WIN,
                     cfg.FIT_SIGMA_MAX, np.inf, np.inf]
                ),
                max_nfev=3000
            )
        except Exception as e:
            return None, f"curve_fit异常 {type(e).__name__}"

        # SNR 门槛：拟合峰高需显著高于拟合残差
        resid = y - model(x, *popt)
        rms = float(np.sqrt(np.mean(resid ** 2)))
        snr = popt[0] / max(rms, 1e-9)
        if snr < cfg.LOCAL_FIT_MIN_SNR:
            return None, f"SNR={snr:.2f}<{cfg.LOCAL_FIT_MIN_SNR}"

        # 相对峰高门槛：避免把噪声当作峰
        rel = popt[0] / max(smax, 1e-9)
        if rel < cfg.LOCAL_FIT_MIN_REL:
            return None, f"相对峰高={rel:.4f}<{cfg.LOCAL_FIT_MIN_REL}"

        return float(popt[1]), None

    # ----------------------------------------------------------
    # 匹配容差
    # ----------------------------------------------------------
    @staticmethod
    def _match_tol(std_wl, cfg):
        if len(std_wl) > 1:
            min_gap = np.min(np.diff(np.sort(std_wl)))
            tol = min(cfg.RANSAC_TOL, 0.6 * min_gap)
        else:
            tol = cfg.RANSAC_TOL
        return max(tol, 1.0)

    # ----------------------------------------------------------
    # 双指针单调贪心匹配 O(n_p + n_s)
    # 保证像素序与波长序同时单调递增
    # ----------------------------------------------------------
    @staticmethod
    def _greedy_match(peak_pos, std_wl, a, b, tol):
        est = a * np.asarray(peak_pos, dtype=float) + b
        n_s = len(std_wl)
        matches = []
        sj = 0
        for pi in range(len(est)):
            if sj >= n_s:
                break
            seg = np.abs(std_wl[sj:] - est[pi])
            off = int(np.argmin(seg))
            if seg[off] <= tol:
                matches.append((pi, sj + off))
                sj = sj + off + 1
        return matches

    @staticmethod
    def _spacing_ok(peak_pos, std_wl, inliers, a, cfg):
        """相邻匹配对的局部色散应与全局 a 一致"""
        if len(inliers) < 4:
            return True
        p_idx = np.array([p for p, _ in inliers], dtype=int)
        s_idx = np.array([s for _, s in inliers], dtype=int)
        dp = np.diff(peak_pos[p_idx])
        dw = np.diff(std_wl[s_idx])
        if np.any(dp <= 0):
            return False
        return bool(np.all(np.abs(dw / dp - a) <= a * cfg.SPACE_CHECK_TOL))

    def _ransac(self, peak_pos, std_wl, tol):
        """
        RANSAC 求初始色散模型：
          1) 随机取两对 (峰, 标准线)，求解 a、b
          2) 用 a、b 约束在合理色散范围内
          3) 贪心匹配统计内点数，保留内点最多的模型
        """
        cfg = self.cfg
        n_p, n_s = len(peak_pos), len(std_wl)
        max_possible = min(n_p, n_s)
        best_in, best_coeff, best_score = [], None, -1

        for _ in range(cfg.RANSAC_ITER):
            if best_score >= max_possible:
                break
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

            inliers = self._greedy_match(peak_pos, std_wl, a, b, tol)
            n_in = len(inliers)
            if n_in < cfg.RANSAC_MIN:
                continue
            if n_in >= 4 and not self._spacing_ok(peak_pos, std_wl,
                                                  inliers, a, cfg):
                continue
            if n_in > best_score:
                best_score, best_coeff, best_in = n_in, (a, b), inliers
        return best_coeff, best_in


# ==============================================================
# 标定计算：单调校验 → 预过滤 → 迭代 sigma 剔除 → 线性拟合
# ==============================================================
class Calibrator:
    def __init__(self, cfg: Config):
        self.cfg = cfg

    def run(self, bands_std, mask_normal, mask_overlap):
        cfg = self.cfg
        n_std = len(cfg.STD_WL)

        # 合并正常与重叠标记
        used_mask = mask_normal | mask_overlap
        bands_std, used_mask = self._monotonic_check(bands_std, used_mask)
        mask_normal = mask_normal & used_mask
        mask_overlap = mask_overlap & used_mask

        valid_idx = np.where(used_mask)[0]
        if len(valid_idx) < 5:
            raise ValueError(f"有效点仅{len(valid_idx)}个，至少需要5个")

        x, y = bands_std[valid_idx].copy(), cfg.STD_WL[valid_idx].copy()

        # 预过滤：剔除粗差 > 15nm 的点，避免污染初值
        init_coeff = np.polyfit(x, y, 1)
        init_res = np.abs(np.polyval(init_coeff, x) - y)
        x, y = x[init_res < 15.0], y[init_res < 15.0]
        if len(x) < 5:
            raise ValueError("预过滤后有效点不足5个")

        # 迭代 sigma 剔除：每轮残差 > max(sigma*OUTLIER_SIGMA, 1.0nm) 的点被剔除
        inlier = np.ones(len(x), dtype=bool)
        for _ in range(cfg.OUTLIER_ITER):
            if np.sum(inlier) < 5:
                break
            c = np.polyfit(x[inlier], y[inlier], 1)
            res = np.abs(np.polyval(c, x) - y)
            sigma = np.std(res[inlier])
            new_mask = res < max(cfg.OUTLIER_SIGMA * sigma, 1.0)
            if np.array_equal(new_mask, inlier):
                break
            inlier = new_mask
        if np.sum(inlier) < 5:
            raise ValueError("剔除异常后有效点不足5个")

        xf, yf = x[inlier], y[inlier]
        c1 = np.polyfit(xf, yf, 1)

        # 回填预测波长与残差（用于 Excel 显示与异常点判定）
        pred = np.full(n_std, np.nan)
        res = np.full(n_std, np.nan)
        for i in range(n_std):
            if used_mask[i] and not np.isnan(bands_std[i]):
                pred[i] = np.polyval(c1, bands_std[i])
                res[i] = abs(cfg.STD_WL[i] - pred[i])

        rf = np.abs(np.polyval(c1, xf) - yf)
        rmse = float(np.sqrt(np.mean(rf ** 2)))
        ss_tot = np.sum((yf - np.mean(yf)) ** 2)
        r2 = float(1 - np.sum(rf ** 2) / ss_tot) if ss_tot > 0 else 1.0
        passed = rmse < cfg.MAX_RMSE and r2 > cfg.MIN_R2

        return {
            "coeff": c1, "bands": bands_std, "pred": pred, "res": res,
            "mask_normal": mask_normal, "mask_overlap": mask_overlap,
            "rmse": rmse, "r2": r2, "passed": passed,
            "n_valid": int(np.sum(inlier)),
            "n_out": int(len(x) - np.sum(inlier))
        }

    def _monotonic_check(self, bands_std, mask):
        """
        检查像素序与波长序是否一致：
        按标准波长排序后若出现像素逆序，剔除残差较大的一点。
        """
        std_wl = self.cfg.STD_WL
        valid_idx = np.where(mask)[0]
        if len(valid_idx) < 3:
            return bands_std, mask
        sorted_idx = valid_idx[np.argsort(std_wl[valid_idx])]
        sorted_bands = bands_std[sorted_idx]
        i = 0
        while i < len(sorted_bands) - 1:
            if sorted_bands[i] >= sorted_bands[i + 1]:
                wl_i, wl_j = std_wl[sorted_idx[i]], std_wl[sorted_idx[i + 1]]
                a_local = ((wl_j - wl_i) / (sorted_bands[i + 1] - sorted_bands[i])
                           if sorted_bands[i + 1] != sorted_bands[i] else 0)
                b_local = wl_i - a_local * sorted_bands[i]
                res_i = abs(wl_i - (a_local * sorted_bands[i] + b_local))
                res_j = abs(wl_j - (a_local * sorted_bands[i + 1] + b_local))
                bad_idx = sorted_idx[i] if res_i > res_j else sorted_idx[i + 1]
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
# 输出层：xlsxwriter 报表（图表样式参照 400-1000nm 版本）
# ==============================================================
class ExcelWriter:
    """
    Excel 结果输出器。
    - 数字格式统一在 num_format 中指定，避免浮点显示位数漂移
    - 异常点判定：D 列残差 > 色散值 |a1| → 整行红字
    - 图表：标定点黑色细实线连接；趋势线蓝色虚线
    - 文本框：无边框、无填充、黑色 Consolas 字体
    """

    _LAMP_BG = {"AR": "#FF9900", "NE": "#00CCFF", "KR": "#00CC00"}
    _GRAY_BG = "#F2F2F2"

    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.wb = None
        self.fmt = {}

    def start(self, path):
        self.wb = xlsxwriter.Workbook(path, {'nan_inf_to_errors': True})
        self.fmt = self._make_formats()

    def close(self):
        if self.wb is not None:
            self.wb.close()
            self.wb = None
            self.fmt = {}

    def _fmt_base(self, **extra):
        d = {'align': 'center', 'valign': 'vcenter',
             'border': 1, 'border_color': '#BFBFBF'}
        d.update(extra)
        return self.wb.add_format(d)

    def _make_formats(self):
        f = {}
        f['header'] = self._fmt_base(bold=True, bg_color='#D6EAF8')
        f['gray'] = self._fmt_base(bg_color=self._GRAY_BG)
        f['coeff_name'] = self._fmt_base(bg_color=self._GRAY_BG)
        f['coeff_val'] = self._fmt_base(bg_color=self._GRAY_BG,
                                        num_format='0.0000')

        # 4 列数据：黑字/红字各一版，均带 num_format 固定显示位数
        f['pixel'] = self._fmt_base(num_format='0.00')
        f['pixel_red'] = self._fmt_base(num_format='0.00', font_color='#FF0000')
        f['overlap_pixel'] = self._fmt_base(bg_color='#FFF2CC', num_format='0.00')
        f['overlap_pixel_red'] = self._fmt_base(bg_color='#FFF2CC',
                                                num_format='0.00',
                                                font_color='#FF0000')
        f['pred'] = self._fmt_base(num_format='0.00')
        f['pred_red'] = self._fmt_base(num_format='0.00', font_color='#FF0000')
        f['res'] = self._fmt_base(num_format='0.000')
        f['res_red'] = self._fmt_base(num_format='0.000', font_color='#FF0000')

        # B 列灯种底色 + 4 位小数
        for lamp, bg in self._LAMP_BG.items():
            f[f'b_{lamp}'] = self._fmt_base(bg_color=bg, num_format='0.0000')
            f[f'b_{lamp}_red'] = self._fmt_base(bg_color=bg, num_format='0.0000',
                                                 font_color='#FF0000')
        return f

    def write_notice(self, warns):
        """标定失败时输出说明页，避免空工作簿"""
        ws = self.wb.add_worksheet("说明")
        ws.set_column(0, 0, 90)
        ws.write(0, 0, "标定未能完成，告警如下：", self.fmt['header'])
        for i, w in enumerate(warns, start=1):
            ws.write(i, 0, "· " + w)

    def write_sheet(self, res, offset):
        cfg = self.cfg
        n_std = len(cfg.STD_WL)
        sheet_name = f"offset{offset}-1000-2500"
        ws = self.wb.add_worksheet(sheet_name)
        fmt = self.fmt

        ws.set_column(0, 0, 14)
        ws.set_column(1, 1, 16)
        ws.set_column(2, 2, 16)
        ws.set_column(3, 3, 14)
        ws.set_column(4, 4, 16)
        ws.set_column(5, 5, 22)

        idx_to_lamp = {}
        for lamp, indices in cfg.LAMP_IDX.items():
            for idx in indices:
                idx_to_lamp[idx] = lamp

        for col, name in enumerate(cfg.HEADER):
            ws.write(0, col, name, fmt['header'])

        valid_idx = [i for i in range(n_std)
                     if res["mask_normal"][i] or res["mask_overlap"][i]]
        c1 = res["coeff"]
        disp = abs(c1[0]) if c1 is not None else None

        # 主数据区
        for i in range(n_std):
            row = i + 1
            is_normal = res["mask_normal"][i]
            is_overlap = res["mask_overlap"][i]
            is_used = is_normal or is_overlap
            has_pred = not np.isnan(res["pred"][i])
            # 异常点判定：绝对残差 > 色散值（|斜率|）
            is_out = (has_pred and disp is not None and res["res"][i] > disp)

            # A 列：像素位置
            if is_used:
                if is_overlap:
                    pixel_fmt = (fmt['overlap_pixel_red'] if is_out
                                 else fmt['overlap_pixel'])
                else:
                    pixel_fmt = fmt['pixel_red'] if is_out else fmt['pixel']
                ws.write_number(row, 0, float(res["bands"][i] + 1), pixel_fmt)
            else:
                ws.write(row, 0, "不可用", fmt['gray'])

            # B 列：标准波长 + 灯种底色
            lamp = idx_to_lamp.get(i, "AR")
            wl_val = float(cfg.STD_WL[i])
            if not is_used:
                ws.write_number(row, 1, wl_val, fmt['gray'])
            else:
                key = f'b_{lamp}_red' if is_out else f'b_{lamp}'
                ws.write_number(row, 1, wl_val, fmt[key])

            # C 列：计算波长 / D 列：绝对残差
            if has_pred:
                ws.write_number(row, 2, float(res["pred"][i]),
                                fmt['pred_red'] if is_out else fmt['pred'])
                ws.write_number(row, 3, float(res["res"][i]),
                                fmt['res_red'] if is_out else fmt['res'])

        # 系数栏（E/F 列）
        for i, (name, val) in enumerate(zip(cfg.COEFF_NAME, res["coeff"])):
            ws.write(i, 4, name, fmt['coeff_name'])
            ws.write_number(i, 5, float(val), fmt['coeff_val'])

        # 图表 + 文本框
        if len(valid_idx) >= 5 and c1 is not None:
            self._add_chart(ws, res, valid_idx, offset, c1)

    def _add_chart(self, ws, res, valid_idx, offset, c1):
        """绘制标定曲线图：标定点黑线连接 + 蓝色虚线趋势线 + 浮动公式文本框"""
        cfg = self.cfg

        # 隐藏数据表，避免污染用户视图
        data_sheet = f"_chart_data_{offset}"
        ws_data = self.wb.add_worksheet(data_sheet)
        ws_data.hide()

        n_pt = len(valid_idx)
        x_arr = [float(res["bands"][i] + 1) for i in valid_idx]
        y_arr = [float(cfg.STD_WL[i]) for i in valid_idx]

        # 写入标定点
        for k, (xa, ya) in enumerate(zip(x_arr, y_arr)):
            ws_data.write_number(k + 1, 0, xa)
            ws_data.write_number(k + 1, 1, ya)

        # 写入趋势线（100 个采样点）
        x_fit = np.linspace(min(x_arr), max(x_arr), 100)
        y_fit = np.polyval(c1, x_fit)
        for k in range(len(x_fit)):
            ws_data.write_number(k + 1, 3, float(x_fit[k]))
            ws_data.write_number(k + 1, 4, float(y_fit[k]))

        chart = self.wb.add_chart({'type': 'scatter'})

        # 系列1：标定点（圆点 + 黑色细实线，参照 400-1000nm 样式）
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

        chart.set_title({'name': f'offset{offset}-1000-2500波长标定曲线'})
        chart.set_x_axis({'name': 'Pixel Location', 'min': 0})
        chart.set_y_axis({'name': 'Wavelength(nm)'})
        chart.set_size({'width': 800, 'height': 500})
        chart.set_legend({'none': True})
        ws.insert_chart('G2', chart, {'x_offset': 0, 'y_offset': 0})

        # 浮动公式文本框：无边框、无填充、黑色 Consolas 字体
        formula_text = f"y = {c1[0]:.4f}x + {c1[1]:.2f}"
        ws.insert_textbox(8, 10, formula_text, {
            'width': 200, 'height': 32,
            'x_offset': 18, 'y_offset': 10,
            'font': {'name': 'Consolas', 'size': 11, 'color': '#000000'},
            'align': {'horizontal': 'center', 'vertical': 'center'},
            'line': {'none': True},
            'fill': {'none': True},
        })


# ==============================================================
# 主控层
# ==============================================================
class App:
    # 词边界匹配，避免 "dark" 误判为氩灯
    LAMP_RE = {L: re.compile(rf'(?<![A-Z]){L}(?![A-Z])')
               for L in ("AR", "NE", "KR")}

    def __init__(self):
        self.cfg = Config()
        self.spec_tool = SpecTool()
        self.matcher = PeakMatcher(self.cfg)
        self.cal = Calibrator(self.cfg)
        self.writer = ExcelWriter(self.cfg)
        self.warn = []

    def run(self, argv=None):
        batch = self._input(argv)
        subs = self._scan_folders()
        print(f"识别到 {len(subs)} 个灯种子文件夹\n")
        n_std = len(self.cfg.STD_WL)

        out_file = os.path.join(self.cfg.OUT_DIR, f"{self.cfg.out_name}.xlsx")
        self.writer.start(out_file)

        try:
            bands_std = np.full(n_std, np.nan)
            mask_normal = np.zeros(n_std, dtype=bool)
            mask_overlap = np.zeros(n_std, dtype=bool)

            # 按平均波长从短到长依次处理，保持处理顺序与物理一致
            lamp_order = sorted(
                self.cfg.LAMP_IDX.keys(),
                key=lambda l: np.mean(self.cfg.STD_WL[self.cfg.LAMP_IDX[l]])
            )
            for lamp in lamp_order:
                for sf in [s for s in subs if s["lamp"] == lamp]:
                    self._process_lamp(sf, bands_std, mask_normal, mask_overlap)

            n_normal = int(np.sum(mask_normal))
            n_overlap = int(np.sum(mask_overlap))
            print("\n---- 全局标定 ----")
            print(f"正常匹配点：{n_normal} 个，重叠峰：{n_overlap} 个")

            if (n_normal + n_overlap) < 5:
                self.warn.append("有效点不足5个")
                print("× 有效点不足，跳过标定\n")
                self.writer.write_notice(self.warn)
            else:
                try:
                    cal_res = self.cal.run(bands_std, mask_normal, mask_overlap)
                except Exception as e:
                    self.warn.append(f"标定失败：{e}")
                    print(f"× 标定失败：{e}\n")
                    self.writer.write_notice(self.warn)
                else:
                    self._print_summary(cal_res)
                    self.writer.write_sheet(cal_res, self.cfg.offset)
                    print("✅ 工作表已添加\n")
        finally:
            self.writer.close()

        print(f"✅ 所有结果已输出：{out_file}")
        if self.warn:
            print("\n" + "-" * 45)
            print("⚠ 告警汇总：")
            for w in self.warn:
                print(f"  · {w}")

        if not batch and sys.stdin.isatty():
            input("\n按回车键退出")

    @staticmethod
    def _print_summary(cal_res):
        """标定摘要只打印到终端，不写入 Excel"""
        c = cal_res["coeff"]
        print("-" * 45)
        print("[标定摘要]")
        print(f"  有效标定点：{cal_res['n_valid']}")
        print(f"  剔除异常点：{cal_res['n_out']}")
        print(f"  RMSE      ：{cal_res['rmse']:.4f} nm")
        print(f"  R²        ：{cal_res['r2']:.6f}")
        print(f"  拟合公式  ：y = {c[0]:.4f}x + {c[1]:.2f}")
        print(f"  判定结果  ：{'✅ 合格' if cal_res['passed'] else '❌ 不合格'}")
        print("-" * 45)

    def _input(self, argv):
        p = argparse.ArgumentParser(description="高光谱波长标定（1000-2500nm）")
        p.add_argument("--root", help="数据根文件夹路径")
        p.add_argument("--offset", type=int, help="offset 整数值")
        p.add_argument("--out", help="输出文件名（不含扩展名）")
        args = p.parse_args(argv)
        batch = args.root is not None and args.offset is not None

        print("=" * 55)
        print("高光谱波长标定程序（1000-2500nm）")
        print("=" * 55)

        root = args.root
        while not root or not os.path.isdir(root):
            if root:
                print("× 路径不存在，请重新输入")
            root = input("请输入数据根文件夹路径：").strip().strip('"').strip("'")
        self.cfg.root = root

        offset = args.offset
        while offset is None:
            try:
                offset = int(input("请输入offset整数值：").strip())
            except ValueError:
                print("× 请输入整数")
        self.cfg.offset = offset

        default_name = f"{datetime.now():%Y%m%d}_1000-2500波长标定结果"
        if batch and not args.out:
            name = ""
        else:
            name = args.out or input(
                f"请输入输出文件名（不含扩展名，直接回车默认 {default_name}）："
            ).strip()
        name = name.replace(".xlsx", "").replace("\\", "").replace("/", "")
        self.cfg.out_name = name if name else default_name

        print("\n开始处理...\n")
        return batch

    def _scan_folders(self):
        subs = []
        for name in os.listdir(self.cfg.root):
            fp = os.path.join(self.cfg.root, name)
            if not os.path.isdir(fp):
                continue
            upper = name.upper()
            lamp = next((L for L, rx in self.LAMP_RE.items()
                         if rx.search(upper)), None)
            if lamp:
                subs.append({"name": name, "path": fp, "lamp": lamp})
        if not subs:
            raise ValueError("未识别到有效灯种子文件夹，请检查命名格式")
        return subs

    def _process_lamp(self, sf, bands_std, mask_normal, mask_overlap):
        print(f"\n▶ {sf['lamp']} 灯段：{sf['name']}")
        hdr_list = sorted(f for f in os.listdir(sf["path"])
                          if f.lower().endswith(".hdr"))
        if not hdr_list:
            self.warn.append(f"{sf['name']} 无HDR文件")
            print("  × 无HDR文件，跳过")
            return

        # 优先选择同名 RAW 存在的 HDR
        def _raw_ok(f):
            base = os.path.join(sf["path"], os.path.splitext(f)[0])
            return os.path.exists(base) or os.path.exists(base + '.raw')

        with_raw = [f for f in hdr_list if _raw_ok(f)]
        pick = with_raw[0] if with_raw else hdr_list[0]
        if len(with_raw) > 1:
            self.warn.append(f"{sf['name']} 存在多个HDR，已使用 {pick}")

        try:
            hdr_path = os.path.join(sf["path"], pick)
            spec, _ = self.spec_tool.read_center(hdr_path, self.cfg.PROFILE_WINDOW)
            print(f"  光谱尺寸：{len(spec)} 波段")
            spec, _ = self.spec_tool.remove_baseline(spec, self.cfg.BASE_ITER)
        except Exception as e:
            self.warn.append(f"{sf['name']} 读取失败：{e}")
            print(f"  × 读取失败：{e}")
            return

        try:
            normal_res, overlap_res = self.matcher.run(spec, sf["lamp"])
            for k, v in normal_res.items():
                bands_std[k] = v
                mask_normal[k] = True
            for k, v in overlap_res.items():
                bands_std[k] = v
                mask_overlap[k] = True
        except Exception as e:
            self.warn.append(f"{sf['lamp']}匹配失败：{e}")
            print(f"  × 匹配失败：{e}")


if __name__ == "__main__":
    App().run()