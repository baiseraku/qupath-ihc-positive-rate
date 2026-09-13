/**
 * IHC 阳性率「诊断三连」：背景 / 色向量夹角 / 阈值扫描
 *
 * 用途：在信任 ihc_summary.csv 的阳性率之前，先判断这个数能不能用。
 * 前提：主脚本用 -s 跑过（细胞测量值已存进 .qpproj 的 data/*.qpdata）。
 *
 * 用法：
 *   QP="<QuPath 可执行文件>"
 *   "$QP" script -p /绝对路径/project.qpproj /绝对路径/ihc_qc_probe.groovy
 *
 * 产出：results/qc_probe_<图名>.txt，并在控制台打印同样内容。
 */
import qupath.lib.regions.RegionRequest
import java.awt.image.BufferedImage

def imageName = getCurrentImageName()
def imageData = getCurrentImageData()
def server = getCurrentServer()
def cal = server.getMetadata().getPixelCalibration()
double pxW = cal.getPixelWidthMicrons(), pxH = cal.getPixelHeightMicrons()
def outDir = buildFilePath(PROJECT_BASE_DIR, 'results')
new File(outDir).mkdirs()
def rep = new StringBuilder()
rep << "###### ${imageName} ######\n"

// ---- 1) 真实背景（扫描仪的白不是 255）----
double ds = Math.max(1.0d, server.getWidth() / 1500.0d)
def sm = server.readRegion(RegionRequest.createInstance(server.getPath(), ds, 0, 0,
        server.getWidth(), server.getHeight()))
def im = new BufferedImage(sm.getWidth(), sm.getHeight(), BufferedImage.TYPE_INT_RGB)
def gg = im.createGraphics(); gg.drawImage(sm, 0, 0, null); gg.dispose()
def rst = im.getRaster()
int w = im.getWidth(), h = im.getHeight()
int[] rr = rst.getSamples(0, 0, w, h, 0, (int[]) null)
int[] g1 = rst.getSamples(0, 0, w, h, 1, (int[]) null)
int[] b1 = rst.getSamples(0, 0, w, h, 2, (int[]) null)
def cnt = new HashMap<Integer, Integer>()
for (int i = 0; i < rr.length; i++) {
    cnt.merge((rr[i] << 16) | (g1[i] << 8) | b1[i], 1, Integer::sum)
}
def top = cnt.entrySet().sort { -it.value }.take(5)
rep << "\n【1】背景众数 Top5（默认色向量假设背景 255，偏差会整体平移 DAB OD）\n"
top.each { e ->
    int k = e.key
    rep << String.format("   RGB(%d,%d,%d) x%d\n", (k >> 16) & 0xff, (k >> 8) & 0xff, k & 0xff, e.value)
}
rep << "   主脚本记录的背景见 ihc_summary.csv 的 backgroundRGB 列\n"

// ---- 2) 分片色向量估计：H-DAB 夹角 ----
def vecOf = { st, int i -> def sv = st.getStain(i); [sv.getRed(), sv.getGreen(), sv.getBlue()] as double[] }
def angleOf = { double[] a, double[] b ->
    double d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    Math.toDegrees(Math.acos(Math.max(-1.0d, Math.min(1.0d, d))))
}
rep << "\n【2】色向量夹角（判据：≥30° 才说明「棕色」是图像里独立的颜色方向）\n"
try {
    double ds2 = Math.max(1.0d, server.getWidth() / 1200.0d)
    def img2 = server.readRegion(RegionRequest.createInstance(server.getPath(), ds2, 0, 0,
            server.getWidth(), server.getHeight()))
    def est = qupath.lib.analysis.algorithms.EstimateStainVectors.estimateStains(
            img2, imageData.getColorDeconvolutionStains(), false)
    def vh = vecOf(est, 1), vd = vecOf(est, 2), vr = vecOf(est, 3)
    double ang = angleOf(vh, vd)
    rep << String.format("   分片估计 H-DAB 夹角 = %.1f°   %s\n", ang,
            ang >= 30.0d ? "→ 正常" : (ang >= 20.0d ? "→ 偏低，存疑" : "→ 退化！几乎没有 DAB 显色"))
    rep << String.format("   估计值 H(%.3f %.3f %.3f) DAB(%.3f %.3f %.3f) Residual(%.3f %.3f %.3f)\n",
            vh[0], vh[1], vh[2], vd[0], vd[1], vd[2], vr[0], vr[1], vr[2])
} catch (Throwable t) {
    rep << "   估计失败: " + t + "\n"
}

// ---- 3) 已存细胞的 DAB OD 分布 + 阈值扫描 ----
def cells = getCellObjects()
rep << "\n【3】已存细胞 " + cells.size() + " 个（需主脚本用 -s 跑过）\n"
if (!cells.isEmpty()) {
    def mnames = cells[0].getMeasurementList().getNames()
    ['Cell: DAB OD mean', 'Nucleus: DAB OD mean', 'Cytoplasm: DAB OD mean'].each { mn ->
        if (!mnames.contains(mn)) { rep << "   " + mn + ": 不存在\n"; return }
        def vals = cells.collect { it.getMeasurementList().get(mn) }.findAll { !Double.isNaN(it) }.sort()
        if (vals.isEmpty()) { rep << "   " + mn + ": 全 NaN\n"; return }
        def q = { double p -> vals[(int) Math.min(vals.size() - 1, Math.max(0, Math.round(p * (vals.size() - 1))))] }
        rep << String.format("\n   == %s ==\n", mn)
        rep << String.format("      中位=%.4f p75=%.4f p90=%.4f p95=%.4f p99=%.4f 最大=%.4f\n",
                q(0.5), q(0.75), q(0.90), q(0.95), q(0.99), vals[-1])
        rep << "      阈值 -> 阳性率:\n"
        [0.02d, 0.05d, 0.10d, 0.15d, 0.20d, 0.30d, 0.40d, 0.50d].each { t ->
            double r = (double) vals.count { double v -> v >= t } / vals.size()
            rep << String.format("        %.2f -> %6.2f%%  %s\n", t, 100 * r, '#' * (int) Math.round(r * 40))
        }
    }
    def ann = getAnnotationObjects().find { it.getROI() != null }
    if (ann != null) {
        rep << String.format("\n   组织注释面积 = %.2f mm²\n", ann.getROI().getArea() * pxW * pxH / 1e6)
    }
}

rep << "\n【判读】\n"
rep << "  · 夹角 <20° → 阳性率不可信，先查染色（阳性对照 / 抗体 / 修复 / 显色时间）\n"
rep << "  · 阈值扫描若在两个相邻阈值间从 ~99% 掉到 ~1% → 分布是窄单峰，\n"
rep << "    阳性率只是「在分布中央切一刀」，组间差异会随阈值漂移，不是阳性细胞比例\n"
rep << "  · 有平台段（相邻阈值阳性率变化平缓）→ 该阈值下的阳性率可信\n"

def f = new File(outDir, 'qc_probe_' + imageName + '.txt')
f.text = rep.toString()
println rep.toString()
