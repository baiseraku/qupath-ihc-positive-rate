/**
 * IHC（组化）批处理：自动组织检测 → 阳性细胞检测 → 每片阳性率
 *
 * 用法：
 *   QP="<QuPath 可执行文件>"
 *   "$QP" script -p /绝对路径/project.qpproj -s /绝对路径/qupath_ihc_positive_rate.groovy
 *   "$QP" script -p /绝对路径/project.qpproj -s /绝对路径/qupath_ihc_positive_rate.groovy \
 *        --args "[200,0.1,0.1,Cell,default,auto]"
 *
 * ⚠️ -p 必须给绝对路径：相对路径会让 QuPath 报 "this.dirBase is null" 而中止。
 * ⚠️ 不给 -s 就不会写回 .qpdata，只出 CSV + QC 图。
 *
 * 参数顺序（都可省，省略用默认）：
 *   [0] TISSUE_THRESHOLD   组织检测灰度阈值，默认 200
 *   [1] NUCLEUS_THRESHOLD  细胞核检测阈值(Hematoxylin OD)，默认 0.1
 *   [2] DAB_POS_THRESHOLD  DAB 阳性阈值(OD)，默认 0.1
 *   [3] COMPARTMENT        阳性判定腔室 Cell|Nucleus|Cytoplasm，默认 Cell
 *   [4] STAIN_MODE         色向量 default|estimate|keep，默认 default
 *   [5] BG_ARG             背景 auto（自动检测）或 "222,222,224"，默认 auto
 *
 * 产出（写在项目目录下 results/）：
 *   ihc_summary.csv   每片一行：细胞数、组织面积、阳性细胞数、阳性率、色向量 QC
 *   qc/*_qc.jpg       QC 叠加图：绿=自动识别组织轮廓，红=判为阳性的细胞
 *
 * 跑完务必看 estimatedAngleDeg 列：<20° 说明这张片没有可用的 DAB 信号，
 * 阳性率不可信（判读规则见 SKILL.md 的《结果判读》）。
 */
import qupath.lib.regions.RegionRequest
import qupath.lib.color.ColorDeconvolutionStains
import qupath.lib.analysis.algorithms.EstimateStainVectors
import javax.imageio.ImageIO
import java.awt.BasicStroke
import java.awt.Color
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage

// ==================== 参数 ====================
def argv = (args == null) ? [] : (args as List)
double TISSUE_THRESHOLD = argv.size() > 0 ? (argv[0] as double) : 200.0
double NUCLEUS_THRESHOLD = argv.size() > 1 ? (argv[1] as double) : 0.1
double DAB_POS_THRESHOLD = argv.size() > 2 ? (argv[2] as double) : 0.1
String COMPARTMENT = argv.size() > 3 ? argv[3].toString() : 'Cell'
String STAIN_MODE = argv.size() > 4 ? argv[4].toString() : 'default'
String BG_ARG = argv.size() > 5 ? argv[5].toString() : 'auto'   // auto 或 "222,222,224"
boolean MAKE_QC = true                 // 输出 QC 叠加图：组织轮廓 + 阳性细胞

// 细胞检测参数（想调就改这里）
double CELL_EXPANSION_UM = 5.0
double DETECT_RESOLUTION_UM = 0.5
double MIN_AREA_UM2 = 10.0
double MAX_AREA_UM2 = 400.0

// QuPath 默认 H-DAB 的染色向量（H-DAB 夹角 36.9°）。
// 注意：默认背景假设 255，但扫描仪实际背景常是 220~230，必须校正，否则 DAB OD 有固定正偏移。
double[] HDAB_H = [0.65103d, 0.70119d, 0.29084d]
double[] HDAB_D = [0.26923d, 0.56754d, 0.77817d]
double[] HDAB_R = [0.71061d, 0.42096d, 0.56473d]

def makeStains = { String nm, double[] vh, double[] vd, double[] vr, int br, int bg, int bb ->
    def j = String.format('{"Name":"%s","Stain 1":"Hematoxylin","Values 1":"%.5f %.5f %.5f ",'
            + '"Stain 2":"DAB","Values 2":"%.5f %.5f %.5f ",'
            + '"Stain 3":"Residual","Values 3":"%.5f %.5f %.5f ",'
            + '"Background":" %d %d %d "}', nm,
            vh[0], vh[1], vh[2], vd[0], vd[1], vd[2], vr[0], vr[1], vr[2], br, bg, bb)
    ColorDeconvolutionStains.parseColorDeconvolutionStainsArg(j)
}

// 取低倍图上出现最多的 RGB 作为背景（扫描仪的白不是 255）
// ⚠️ svr 必须显式传参：Groovy 闭包前向引用后面才 def 的局部变量会失败，
//    而且失败会被 catch 吞掉、静默回退到 255，结果整体偏乐观。
def detectBackground = { svr ->
    try {
        double ds = Math.max(1.0d, svr.getWidth() / 1500.0d)
        def sm = svr.readRegion(RegionRequest.createInstance(svr.getPath(), ds, 0, 0,
                svr.getWidth(), svr.getHeight()))
        def im = new BufferedImage(sm.getWidth(), sm.getHeight(), BufferedImage.TYPE_INT_RGB)
        def g = im.createGraphics(); g.drawImage(sm, 0, 0, null); g.dispose()
        def rst = im.getRaster()
        int w = im.getWidth(), h = im.getHeight()
        int[] rr = rst.getSamples(0, 0, w, h, 0, (int[]) null)
        int[] gg = rst.getSamples(0, 0, w, h, 1, (int[]) null)
        int[] bb = rst.getSamples(0, 0, w, h, 2, (int[]) null)
        def cnt = new HashMap<Integer, Integer>()
        for (int i = 0; i < rr.length; i++) {
            int key = (rr[i] << 16) | (gg[i] << 8) | bb[i]
            cnt.merge(key, 1, Integer::sum)
        }
        def top = cnt.entrySet().sort { -it.value }.take(20)
        // 取前 20 众数的均值，比单个众数稳
        double sr = 0, sg = 0, sb = 0
        top.each { e -> sr += (e.key >> 16) & 0xff; sg += (e.key >> 8) & 0xff; sb += e.key & 0xff }
        return [(int) Math.round(sr / top.size()), (int) Math.round(sg / top.size()), (int) Math.round(sb / top.size())] as int[]
    } catch (Throwable t) {
        logger.error("背景检测失败，回退到 255（结果会偏乐观！）: {}", t.toString())
        return [255, 255, 255] as int[]
    }
}

// ==================== 工具 ====================
def vecOf = { st, int i ->
    def sv = st.getStain(i)
    return [sv.getRed(), sv.getGreen(), sv.getBlue()] as double[]
}
def angleOf = { double[] a, double[] b ->
    double d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    Math.toDegrees(Math.acos(Math.max(-1.0d, Math.min(1.0d, d))))
}
def csvCell = { Object v ->
    if (v == null) return ''
    def s = v.toString()
    (s.contains(',') || s.contains('"')) ? '"' + s.replace('"', '""') + '"' : s
}

// ==================== 主流程 ====================
def imageData = getCurrentImageData()
def server = getCurrentServer()
def imageName = getCurrentImageName()
def cal = server.getMetadata().getPixelCalibration()
double pxW = cal.getPixelWidthMicrons(), pxH = cal.getPixelHeightMicrons()
def resultsDir = buildFilePath(PROJECT_BASE_DIR, 'results')
new File(resultsDir).mkdirs()
def t0 = System.currentTimeMillis()

// ---- 1. 图像类型 ----
setImageType('BRIGHTFIELD_H_DAB')

// ---- 2. 色向量：背景校正 + 按模式选择，并记录分片估计值做 QC ----
def stainsKeep = imageData.getColorDeconvolutionStains()
def bg
if (BG_ARG == 'auto') {
    bg = detectBackground(server)
} else {
    def parts = BG_ARG.split(',')
    bg = [parts[0].trim() as int, parts[1].trim() as int, parts[2].trim() as int] as int[]
}
int bgR = bg[0], bgG = bg[1], bgB = bg[2]
def stainsDefault = makeStains('H-DAB bg-corrected', HDAB_H, HDAB_D, HDAB_R, bgR, bgG, bgB)

// 分片估计色向量，只取 H-DAB 夹角当 QC 指标（不复用它，除非 STAIN_MODE=estimate）
def estimateAngle = { ColorDeconvolutionStains base ->
    try {
        double ds = Math.max(1.0d, server.getWidth() / 1200.0d)
        def img = server.readRegion(RegionRequest.createInstance(server.getPath(), ds, 0, 0,
                server.getWidth(), server.getHeight()))
        def est = EstimateStainVectors.estimateStains(img, base, false)
        def vh = vecOf(est, 1), vd = vecOf(est, 2)
        return [angle: angleOf(vh, vd), stains: est, h: vh, d: vd]
    } catch (Throwable t) {
        return [angle: Double.NaN, stains: null, h: null, d: null]
    }
}

def estInfo = estimateAngle(stainsDefault)
def stainsUse
switch (STAIN_MODE) {
    case 'estimate':
        stainsUse = estInfo.stains ?: stainsDefault
        break
    case 'keep':
        stainsUse = stainsKeep
        break
    default:
        stainsUse = stainsDefault
}
imageData.setColorDeconvolutionStains(stainsUse)

// ---- 3. 清空旧对象，组织检测 ----
resetSelection()
removeAnnotations()
removeDetections()
fireHierarchyUpdate()
runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2',
        '{"threshold": ' + TISSUE_THRESHOLD + ', "requestedPixelSizeMicrons": 20.0,' +
                ' "minAreaMicrons": 10000.0, "maxHoleAreaMicrons": 1000000.0,' +
                ' "darkBackground": false, "smoothImage": true, "medianCleanup": true,' +
                ' "dilateBoundaries": false, "smoothCoordinates": true,' +
                ' "excludeOnBoundary": false, "singleAnnotation": true}')
def anns = getAnnotationObjects()
// ⚠️ ROI.getArea() 的单位是「全分辨率像素²」，不是 µm²，必须乘像素面积换算
double tissueAreaMM2 = anns.sum { it.getROI().getArea() } * pxW * pxH / 1e6

// ---- 4. 阳性细胞检测（只在组织注释内）----
def cells = []
if (anns.isEmpty() || tissueAreaMM2 <= 0) {
    logger.warn("{}: 组织检测为空，跳过细胞检测", imageName)
} else {
    selectAnnotations()
    runPlugin('qupath.imagej.detect.cells.PositiveCellDetection',
            '{"detectionImageBrightfield": "Hematoxylin OD",' +
                    ' "requestedPixelSizeMicrons": ' + DETECT_RESOLUTION_UM + ',' +
                    ' "backgroundRadiusMicrons": 8.0, "backgroundByReconstruction": true,' +
                    ' "medianRadiusMicrons": 0.0, "sigmaMicrons": 1.5,' +
                    ' "minAreaMicrons": ' + MIN_AREA_UM2 + ', "maxAreaMicrons": ' + MAX_AREA_UM2 + ',' +
                    ' "threshold": ' + NUCLEUS_THRESHOLD + ', "maxBackground": 2.0,' +
                    ' "watershedPostProcess": true, "excludeDAB": false,' +
                    ' "cellExpansionMicrons": ' + CELL_EXPANSION_UM + ', "includeNuclei": true,' +
                    ' "smoothBoundaries": true, "makeMeasurements": true,' +
                    ' "thresholdCompartment": "' + COMPARTMENT + ': DAB OD mean",' +
                    ' "thresholdPositive1": ' + DAB_POS_THRESHOLD + ',' +
                    ' "thresholdPositive2": 0.4, "thresholdPositive3": 0.6,' +
                    ' "singleThreshold": true}')
    cells = getCellObjects()
}
resetSelection()
fireHierarchyUpdate()

// ---- 5. 阳性率统计 ----
int nCells = cells.size()
int nPos = cells.count { it.getPathClass() != null && it.getPathClass().toString() == 'Positive' }
double posRate = nCells > 0 ? (double) nPos / nCells : Double.NaN

// ---- 6. 把汇总写到组织注释上（方便在 GUI 里直接看）----
anns.each { a ->
    def ml = a.getMeasurementList()
    ml.put('IHC: nCells', nCells as double)
    ml.put('IHC: tissueAreaMM2', tissueAreaMM2)
    ml.put('IHC: positiveRate', posRate)
}

// ---- 7. 追加/更新汇总 CSV（同一张图重跑会覆盖旧行，不产生重复）----
def header = ['image', 'nCells', 'tissueAreaMM2', 'nPositive', 'positiveRate',
              'dabPosThreshold', 'compartment', 'nucleusThreshold', 'tissueThreshold',
              'backgroundRGB', 'estimatedAngleDeg', 'elapsedSec'] as List
def row = [imageName, nCells, String.format('%.3f', tissueAreaMM2),
           nPos, nCells > 0 ? String.format('%.5f', posRate) : '',
           DAB_POS_THRESHOLD, COMPARTMENT, NUCLEUS_THRESHOLD, TISSUE_THRESHOLD,
           String.format('%d;%d;%d', bgR, bgG, bgB), String.format('%.1f', estInfo.angle),
           String.format('%.1f', (System.currentTimeMillis() - t0) / 1000.0d)]

def summaryFile = new File(resultsDir, 'ihc_summary.csv')
def lines = summaryFile.exists() ? summaryFile.readLines() : []
def kept = lines.findAll { it.trim() && !it.startsWith('image,') }
kept.removeAll { it.startsWith(csvCell(imageName) + ',') }
def outLines = [header.join(',')]
kept.each { outLines << it }
outLines << row.collect(csvCell).join(',')
summaryFile.text = outLines.join('\n') + '\n'

// ---- 8. QC 叠加图：低倍底图 + 组织轮廓(绿) + 阳性细胞(红) ----
if (MAKE_QC) {
    try {
        double qds = Math.max(1.0d, server.getWidth() / 1400.0d)
        def sm = server.readRegion(RegionRequest.createInstance(server.getPath(), qds, 0, 0,
                server.getWidth(), server.getHeight()))
        def canvas = new BufferedImage(sm.getWidth(), sm.getHeight(), BufferedImage.TYPE_INT_RGB)
        def g = canvas.createGraphics()
        g.drawImage(sm, 0, 0, null)
        def at = AffineTransform.getScaleInstance(1.0d / qds, 1.0d / qds)
        g.setStroke(new BasicStroke(3f))
        g.setColor(new Color(0, 170, 0))
        anns.each { a -> g.draw(at.createTransformedShape(a.getROI().getShape())) }
        int dot = Math.max(2, (int) Math.round(3.0d / qds))
        getCellObjects().each { c ->
            if (c.getPathClass() != null && c.getPathClass().toString() == 'Positive') {
                int x = (int) (c.getROI().getCentroidX() / qds)
                int y = (int) (c.getROI().getCentroidY() / qds)
                g.setColor(Color.RED)
                g.fillOval(x - dot, y - dot, dot * 2, dot * 2)
            }
        }
        g.dispose()
        def qcDir = buildFilePath(PROJECT_BASE_DIR, 'results', 'qc')
        new File(qcDir).mkdirs()
        ImageIO.write(canvas, 'JPEG', new File(qcDir, imageName + '_qc.jpg'))
    } catch (Throwable t) {
        logger.warn("QC 图生成失败: {}", t.toString())
    }
}

// ---- 9. 控制台一行摘要 ----
def msg = String.format("[IHC] %s | 组织 %.2f mm² | 细胞 %d | 阳性 %d | 阳性率 %.2f%% (阈 %.2f, %s) | 背景RGB(%d,%d,%d) | 夹角 %.1f° | %.0fs",
        imageName, tissueAreaMM2, nCells, nPos, 100 * posRate, DAB_POS_THRESHOLD, COMPARTMENT,
        bgR, bgG, bgB, estInfo.angle, (System.currentTimeMillis() - t0) / 1000.0d)
println msg
logger.info(msg)

// ---- 10. QC 报警：分片估计的 H-DAB 夹角过小 = 这张片没有可用的 DAB 信号 ----
if (!Double.isNaN(estInfo.angle) && estInfo.angle < 20.0d) {
    def warn = String.format("[IHC-QC] %s 分片估计的 H-DAB 夹角仅 %.1f°（正常应 >30°）：" +
            "这张片很可能几乎没有 DAB 显色，阳性率不可信，请先复查染色。", imageName, estInfo.angle)
    println warn
    logger.warn(warn)
}
