---
name: qupath-ihc-positive-rate
description: >-
  用 QuPath 0.7 无头模式批量统计免疫组化（IHC / DAB）切片阳性率：自动组织检测 →
  阳性细胞检测 → 每片一行 CSV + QC 叠加图（组织轮廓 + 阳性细胞）。内置两道可信度闸门：
  扫描仪背景色校正（默认色向量假设背景 255，实际常为 220–230，不校正会让 DAB OD 整体上移）
  和分片色向量退化报警（H–DAB 夹角小于 20° 说明这张片没有可用的 DAB 信号）。附诊断脚本做阈值扫描。
whenToUse: >-
  用户提出「批量统计组化/免疫组化切片的阳性率」「QuPath 批处理 DAB 阳性细胞」「这批 IHC
  片子阳性率多少」「把组化的阳性率跑出来」，或要在已有 QuPath 项目上批量做组织检测 + 阳性
  细胞计数时。动手前先读本 skill，尤其《结果判读》——阳性率算得出来不等于算得对。
metadata:
  version: 1.0.0
  source: 实操沉淀（IHC 全片批量流程端到端跑通；含背景校正与色向量退化两条结论的发现过程）
disable-model-invocation: false
---

# QuPath 批量 IHC 阳性率

在**已建好的 QuPath 项目**上批量跑：自动组织检测 → 阳性细胞检测 → 每片一行阳性率 CSV，
外加 QC 叠加图。全程无头，不需要开 GUI。

前置：切片已被 QuPath 读得动（OME-TIFF / SVS / NDPI 等；SDPC 需先转，见 `sdpc-to-ometiff` skill），
且已注册进 `.qpproj` 项目。

## 环境事实

- QuPath **0.7.0**，自带 JRE，**不需要系统装 Java**
- 脚本件：`<本 skill 目录>/scripts/qupath_ihc_positive_rate.groovy`（主流程）
- 诊断件：`<本 skill 目录>/scripts/ihc_qc_probe.groovy`（背景 / 夹角 / 阈值扫描）
- QuPath 可执行文件的位置随平台不同：macOS 是
  `.../QuPath-<版本>-<架构>.app/Contents/MacOS/<与应用同名的可执行文件>`，
  Windows 是 `QuPath.exe`，Linux 是解包目录下的 `QuPath`。
  **本 skill 不记录任何机器专属路径**，用前先定位一次再填进下面命令的 `QP=`
- ⚠️ QuPath 0.7 的脚本上下文里**没有** `estimateStainVectors()` 这个便利方法（别照抄网上老脚本）。
  估色向量要用 `qupath.lib.analysis.algorithms.EstimateStainVectors.estimateStains(img, stains, checkColors)`，
  类在 `qupath-core-processing`。同理 `getColorDeconvolutionStains()` 挂在 **ImageData** 上，
  不在 `ImageServerMetadata` 上。

## 一条命令跑完

```bash
QP=<QuPath 可执行文件>                      # 定位方法见上面《环境事实》
S=<本 skill 目录>/scripts/qupath_ihc_positive_rate.groovy

"$QP" script -p /abs/项目/project.qpproj -s "$S"
# 调参：--args "[组织阈值, 核阈值, DAB阈值, 腔室, 色向量模式, 背景]"
"$QP" script -p /abs/项目/project.qpproj -s "$S" --args "[200,0.1,0.1,Cell,default,auto]"
```

- `-p` **必须绝对路径**。给相对路径 QuPath 会抛 `this.dirBase is null` 直接退出
- `-s` 把注释与测量写回 `data/*.qpdata`，之后在 GUI 里改阈值不必重跑。不加 `-s` 就只出 CSV + QC 图
- 脚本按项目逐图执行；每图约 25–40 秒（10 万细胞、10 万×8 万像素量级），随细胞数线性增长
- 产物落在**项目目录**下 `results/`，不在 skill 目录

## 跑之前先定三件事

### 1. 组织阈值 `TISSUE_THRESHOLD`（默认 200）

组织检测插件按灰度阈值分割，取**灰度直方图的谷底**。切片的直方图通常是干净双峰
（组织峰 / 背景峰，中间有明显谷）。先测一次：

```groovy
// 低倍图灰度直方图（跑一次，看谷底落在哪）
double ds = Math.max(1.0, getCurrentServer().getWidth() / 1500.0)
def img = getCurrentServer().readRegion(qupath.lib.regions.RegionRequest.createInstance(
        getCurrentServer().getPath(), ds, 0, 0, getCurrentServer().getWidth(), getCurrentServer().getHeight()))
int[] hist = new int[256]
for (int y = 0; y < img.getHeight(); y++)
    for (int x = 0; x < img.getWidth(); x++) {
        int rgb = img.getRGB(x, y)
        hist[((rgb >> 16 & 0xff) + (rgb >> 8 & 0xff) + (rgb & 0xff)).intdiv(3)]++
    }
for (int i = 0; i < 256; i += 5) {
    int s = 0; for (int k = i; k < Math.min(256, i + 5); k++) s += hist[k]
    if (s > 0) println String.format("%3d-%3d %8d %s", i, i + 4, s, '#' * Math.min(60, (int)(s / 5000)))
}
```

阈值取谷底（双峰之间最小值）。取低了会把淡染组织切掉，取高了会吃进背景和杂质点。
**改完一定要看 QC 图核对**，别只看 CSV。

### 2. 阳性阈值 `DAB_POS_THRESHOLD`（默认 0.1）

**这是全流程唯一没有客观依据、必须由你决定、也最容易出错的地方。**

- 有**阴性对照**（同批、同抗体、省略一抗）：阈值 = 阴性片 `Cell: DAB OD mean` 的 p99
- 没有阴性对照：**不要报绝对阳性率**。改报阈值曲线（诊断脚本可出），或改用连续指标
- 不要在不同批次之间沿用同一个阈值而不做检查——扫描仪、显色时间、抗体批号都会平移整个分布

### 3. 腔室 `COMPARTMENT`（默认 `Cell`）

- 核标记（Ki67、ER、p53…）→ `Nucleus`
- 胞浆 / 膜标记（CD3、CD20、Her2…）→ `Cell` 或 `Cytoplasm`

`Cell` 是「核 + 5 µm 扩张」的整个细胞，均值会被周围空白稀释；胞浆标记用 `Cytoplasm` 更聚焦。

## 结果判读（先看这里，再看数字）

阳性率算得出来 ≠ 算得对。跑完先查 `ihc_summary.csv` 的 **`estimatedAngleDeg`** 列，
再用 `<skill>/scripts/ihc_qc_probe.groovy` 做阈值扫描。

### 闸门一：色向量夹角 `estimatedAngleDeg` ≥ 30°

脚本对每张片独立估一次色向量，把 H 与 DAB 两个方向向量的夹角写进 CSV。
夹角度量的是「棕色在图像里是不是一个独立的颜色方向」：

| 夹角 | 判读 |
|---|---|
| ≥ 30° | 正常，DAB 可以被解卷积分离出来 |
| 20–30° | 偏低，存疑，结合染色形态一起看 |
| **< 20°** | **退化**。估计器只能把苏木素的轴硬劈成两半，DAB OD 基本是串扰+背景偏移，**阳性率不可信** |

夹角退化说明**这批片没有可用的 DAB 显色**，此时调阈值、换腔室、加对照都救不回来，
必须回去查染色（阳性对照 / 抗体浓度 / 抗原修复 / DAB 显色时间）。

### 闸门二：阈值扫描要有平台段

用诊断脚本出「阈值 → 阳性率」曲线。判读：

- **有平台段**（相邻阈值间阳性率变化平缓）→ 该阈值下的阳性率可信
- **一路陡降**（如 0.05→99%、0.10→59%、0.15→8%、0.20→1%）→ 分布是**窄单峰**，
  阈值正好切在分布中央，所谓阳性率只是「哪张片整体染得深一点」的归一化读数

实测案例（3 张片、同一批细胞只换阈值，数值原样保留）：

| 阈值 | 片 1 | 片 2 | 片 3 | 片间差异 |
|---|---|---|---|---|
| 0.05 | 99.6% | 99.5% | 97.3% | 几乎没有 |
| **0.10** | **58.9%** | **56.4%** | **14.2%** | **4 倍** |
| 0.15 | 7.9% | 6.9% | 3.4% | 2 倍 |
| 0.20 | 1.0% | 1.8% | 1.4% | 消失，甚至反向 |
| 0.30 | 0% | 0% | 0% | — |

**「4 倍差异」只存在于 0.10 这一个点上，随阈值左右漂移就消失——这是阈值假象，不是生物学差异。**
同一批切片的分片夹角全部落在 10.9°–15.2°，即全部退化。

判据补充：正常 DAB 阳性细胞的 `Cell: DAB OD mean` 应在 0.3–0.6 量级。
若 p99 只有 0.20 上下、0.30 以上一个细胞都没有，说明**分布右端是空的**，同样是信号缺失的证据。

## Groovy / QuPath API 坑（实测踩过，照抄可避）

| 现象 | 原因与写法 |
|---|---|
| `Cannot invoke "java.io.File.toPath()" because "this.dirBase" is null` | `-p` 给了相对路径。改绝对路径 |
| `clearAnnotations()` / `clearDetections()` 报警告 | 0.7 已废弃，改用 `removeAnnotations()` / `removeDetections()`，再 `fireHierarchyUpdate()` |
| 组织面积大得离谱（占全片 1600%） | `ROI.getArea()` 单位是**全分辨率像素²**，要乘 `pixelWidthMicrons * pixelHeightMicrons` 才是 µm² |
| `No signature of method: getBounds()` | QuPath 的 ROI 没有 `getBounds()` / `getInteriorPoint()`，用 `getBoundsX()/getBoundsY()/getBoundsWidth()/getBoundsHeight()` |
| `No signature of method: getMeasurementValue(String)` | 0.7 的测量表是 `NumericMeasurementList`，用 `.get(name)` 取值、`.getNames()` 取名字 |
| 闭包里的变量莫名报 `MissingPropertyException`，且被 catch 吞掉 | Groovy 闭包**不能前向引用**后面才 `def` 的局部变量。闭包定义在变量声明之前就会失败——把该变量当**显式参数**传进去 |
| `String.format("%d", x)` 抛 `d != java.math.BigDecimal` | Groovy 里 `Integer / Integer` 结果是 BigDecimal。用 `.intdiv()` 或先转 double |
| 无法用 `def` 作 map key / 属性名 | `def` 是 Groovy 关键字，`[(def): x]`、`m.def` 都编译不过。换名字 |
| `$list.size` 在 GString 里变成 `值.size` | `$a.b` 解析为 `${a}.b`。要写 `${list.size()}` |
| `0.0 ?: fallback` 返回了 fallback | Groovy 的 elvis 把 `0.0` 当假值。用 `Double.isNaN(x) ? ... : ...` 显式判断 |
| `double[].count { ... }` 结果恒为 0 | 基本类型数组上的 `count` 闭包不按预期工作。改用显式 for 循环计数 |
| 背景检测静默回退到 255，结果整体偏乐观 | catch 里只 return 不报错。**任何回退都要 logger.error/warn 出声** |
| 图片读出来是 `TYPE_CUSTOM` 之类，`getRaster()` 取通道错位 | 读回后先 `new BufferedImage(w, h, TYPE_INT_RGB)` + `drawImage` 转一次，再用 `raster.getSamples(..., band, null)` 取 R/G/B |

## 输出文件

| 文件 | 内容 |
|---|---|
| `results/ihc_summary.csv` | 每片一行：`image, nCells, tissueAreaMM2, nPositive, positiveRate, dabPosThreshold, compartment, nucleusThreshold, tissueThreshold, backgroundRGB, estimatedAngleDeg, elapsedSec`。同一张图重跑会**覆盖旧行**，不产生重复 |
| `results/qc/<图名>_qc.jpg` | QC 叠加图：绿=自动识别组织轮廓，红=判为阳性的细胞 |
| `results/qc_probe_<图名>.txt` | 诊断脚本产物：背景众数、色向量夹角、DAB OD 分位数、阈值扫描表 |
| `data/*.qpdata` | 只有加 `-s` 才写；含注释、细胞与全部测量值，可在 GUI 里改阈值而不重跑检测 |

## 验收

1. 控制台每张片打一行 `[IHC] …`，无 `ERROR`；退出码 0
2. 夹角 <20° 的片会有 `[IHC-QC]` 警告 —— 有警告就别直接报数
3. **逐张看 QC 图**：绿线是否贴合组织（多块组织、卷绕切片尤其要查）、红点是否落在合理的组织学位置。
   只看 CSV 不看图，组织检测切歪了也发现不了
4. `nCells / tissueAreaMM2`（细胞密度）在同批同组织之间应大致同量级；
   某张片密度差一个数量级，通常是组织检测或细胞检测出问题，不是生物学
5. 报数时说明阈值、腔室、是否做了背景校正和夹角检查
