# qupath-ihc-positive-rate

用 QuPath 0.7 无头模式批量统计免疫组化（IHC / DAB）切片阳性率的可复用流程。
从**已建好的 QuPath 项目**出发，逐张完成：自动组织检测 → 阳性细胞检测 → 每片一行
阳性率 + IOD/MOD（积分与平均光密度）CSV + QC 叠加图（组织轮廓 + 阳性细胞），
全程命令行，不需要开 GUI。

## 项目简介

数字病理 IHC 定量里最容易出错的不是细胞检测，而是**颜色本身**：QuPath 默认色向量假设扫描背景为
纯白 255，而实际扫描仪背景常在 220–230；一旦背景没校正，DAB 光密度（OD）会整体上移，
"阳性率"随之虚高。更进一步，若切片本身 DAB 显色很弱，色解卷积会退化——此时无论阈值怎么调，
算出来的都不是阳性细胞比例。

本流程自动检测扫描背景并校正色向量，同时逐片估计 H–DAB 夹角作为**可信度标注**
（夹角小于 20° 时在 CSV 的 `reliability` 列标为 `degenerate_angle`）。
**流程始终跑完并输出结果，不会因为任何检查不通过而中止**——发现不可靠时的正确做法是把数报出来
并明确标注，而不是不跑。附带的诊断脚本可一次性输出背景众数、色向量夹角、DAB OD 分布与阈值
扫描表，用于判断阳性率是否可用。

## 文件清单

| 文件 | 用途 |
|---|---|
| `SKILL.md` | 流程说明：环境要求、命令行用法、三个必调参数、结果判读规则、Groovy/QuPath API 坑、验收清单 |
| `scripts/qupath_ihc_positive_rate.groovy` | **主脚本**。按 QuPath 项目批量执行：图像类型设定 → 背景校正 + 色向量选择 → 自动组织检测 → 阳性细胞检测 → IOD/MOD 积分 → 写 CSV 与 QC 图 |
| `scripts/ihc_qc_probe.groovy` | **诊断脚本**。输出背景 RGB 众数、分片估计的 H–DAB 夹角、三种腔室的 DAB OD 分位数与「阈值 → 阳性率」扫描表 |

## 运行顺序

两个脚本都是 QuPath 的 Groovy 脚本，以项目为单位运行（`-p` 指定项目，脚本会自动对项目内每张图执行）。

**第 1 步：主脚本**（输入：QuPath 项目 + 已注册的切片；输出：`results/ihc_summary.csv` 与 `results/qc/*.jpg`）

```bash
QP=<QuPath 可执行文件>
S=scripts/qupath_ihc_positive_rate.groovy

"$QP" script -p /abs/项目/project.qpproj -s "$S"
# 可选参数：[组织阈值, 核检测阈值, DAB 阳性阈值, 腔室, 色向量模式, 背景]
"$QP" script -p /abs/项目/project.qpproj -s "$S" --args "[200,0.1,0.1,Cell,default,auto]"
```

- `-p` 必须绝对路径，给相对路径 QuPath 会报 `this.dirBase is null` 中止
- `-s` 会把注释与测量写回 `data/*.qpdata`，之后在 GUI 里改阈值不必重跑检测
- 每张图约 25–40 秒（10 万细胞、10 万 × 8 万像素量级），随细胞数线性增长

**第 2 步：诊断脚本**（输入：第 1 步用 `-s` 存好的项目；输出：`results/qc_probe_<图名>.txt`）

```bash
"$QP" script -p /abs/项目/project.qpproj scripts/ihc_qc_probe.groovy
```

**关于 IOD**：除阳性率外，每片另出 `IOD`（积分光密度，OD·µm²）与 `MOD`（平均光密度）。
这组指标不依赖阳性阈值，改 `PIXEL_OD_THRESHOLD` 不影响 `IOD`/`MOD`，适合直接做组间检验。
⚠️ **IOD 随组织面积缩放**，各片面积不同时要比 `MOD`（已除掉面积），别直接比 `IOD`。

**第 3 步：判读**。先看 `ihc_summary.csv` 的 `estimatedAngleDeg` 与 `reliability` 列
（≥30° 为 `ok`），再看诊断脚本输出的阈值扫描表是否有平台段。两条都通过，阳性率可直接报告；
不通过则照常报数，但必须一并声明结果不可靠及原因。判读细则见 `SKILL.md` 的《结果判读》。

**关于阈值**：`DAB_POS_THRESHOLD` 没有客观默认值，脚本里的 0.1 只是占位。
CSV 的 `thresholdSource` 列会记录该阈值是 `explicit`（显式指定）还是 `default`（用了占位值）。
理想情况下阈值应由阴性对照片定出（阴性片 `Cell: DAB OD mean` 的 p99）。

## 安装为 skill

仓库根目录即是 skill 包，clone 后复制到 skills 目录即可被自动识别：

```bash
git clone git@github.com:baiseraku/qupath-ihc-positive-rate.git
cp -R qupath-ihc-positive-rate ~/.dsh/skills/
```

## 数据说明

- **切片与项目文件不入库**：`.qpproj`、`data/*.qpdata`（单张可达上百 MB）、`classifiers/`
  均在 `.gitignore` 中排除
- **原始切片**（OME-TIFF / SVS 等）需放回本机项目注册时记录的路径，否则项目无法打开
- 复现需要：QuPath 0.7.0（自带 JRE，无需系统装 Java）、已注册切片的 QuPath 项目、
  以及本仓库的两个脚本
- 脚本运行会在**项目目录**下生成 `results/`，与仓库目录无关；仓库内的 `results/` 规则
  仅用于避免误提交本地运行产物
- 版本：QuPath 0.7.0。脚本用到 `EstimateStainVectors`（位于 `qupath-core-processing`），
  该 API 在更早版本中未必存在
