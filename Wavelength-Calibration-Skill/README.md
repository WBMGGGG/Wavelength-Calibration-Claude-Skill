\# Wavelength Calibration - Claude Custom Skill

高光谱相机波长标定自动化 Skill，适配 Claude Custom Skill。ENVI 原生 IDL 精准读取 hdr/dat，自动生成标准化多 Sheet 标定 Excel。



\## 功能特性

\- ✅ 单触发词：`波长标定`

\- ✅ 自动识别日期根文件夹下「灯种+bin模式」子文件夹，解析灯种与 bin 模式

\- ✅ ENVI 原生 IDL 脚本读取栅格数据，±2nm 峰归属判定，精度与手动操作一致

\- ✅ 自动生成多 Sheet 标定 Excel，全表居中、B 列按灯段分色、线性拟合图仅显公式

\- ✅ 3 阶多项式标定，线性斜率质量判定

\- ✅ 重叠峰自动标记，缺失谱线 DEBUG 告警



\## 文件结构

| 文件路径 | 说明 |

|---------|------|

| `Wavelength\_Calibration\_Skill.md` | Claude Custom Skill 主文件，直接导入使用 |

| `scripts/ENVI\_WavelengthCalib\_FindPeaks.pro` | ENVI IDL 脚本，精准读取 hdr/dat 并检测峰值 |

| `scripts/Matlab\_WavelengthCalib\_GenerateExcel.m` | Matlab 脚本，回填、拟合、生成带格式的标定 Excel |



\## 使用步骤

\### 1. 导入 Skill

复制 `Wavelength\_Calibration\_Skill.md` 全部内容 → 打开 Claude → 进入 Custom Skills → 新建 Skill → 粘贴保存。



\### 2. 数据准备

\- 新建\*\*日期命名的根文件夹\*\*，例如 `20260907\_calib`

\- 根文件夹内放入多个\*\*「灯种+bin模式」命名的子文件夹\*\*，例如 `HG-1bin`、`NE-1bin`、`AR-1bin`

\- 每个子文件夹内放置对应 `.hdr` + `.dat` 数据文件



\### 3. ENVI 峰检测

1\. 打开 ENVI 软件

2\. 菜单：File → Open → 选择 `scripts/ENVI\_WavelengthCalib\_FindPeaks.pro`

3\. 按提示输入日期根文件夹路径

4\. 运行完成后，复制控制台输出的像素结果



\### 4. Matlab 生成 Excel

1\. 打开 Matlab，运行 `scripts/Matlab\_WavelengthCalib\_GenerateExcel.m`

2\. 按提示输入 offset 整数值、ENVI 输出的像素数据

3\. 运行完成后，桌面生成 `YYYYMMDD.xlsx` 标定文件



\## 环境要求

\- ENVI 5.x 及以上（支持 IDL 脚本）

\- MATLAB R2020b 及以上

\- Windows 系统



\## 注意事项

\- offset 值必须用户提供，脚本无硬编码

\- 重叠峰不参与拟合，但 C/D 列仍计算

\- 有效标定点 ≥ 5 才输出正式标定结果

\- 图表仅显示线性关系公式，不显示 R²



