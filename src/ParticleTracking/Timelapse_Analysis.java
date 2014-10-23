package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import UtilClasses.GenUtils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.Straightener;
import ij.plugin.TextReader;
import ij.plugin.filter.GaussianBlur;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import ui.ResultsPreviewInterface;
import ui.UserInterface;

/**
 * Timelapse_Analysis seeks to identify individual particles in a stack of
 * images and plot their trajectories. Particles are first identified in
 * individual colour bands at different time-points by searching for local
 * maxima above a certain threshold, fitting IsoGaussian curves about these
 * maxima, associating detected particles between images and constructing
 * trajectories.
 *
 * @author David J Barry
 * @version 2.0, FEB 2011
 */
public class Timelapse_Analysis implements PlugIn {

//    protected static double hystDiff = 1.25;
    protected double xySigEst; //Initial estimate of standard deviation for IsoGaussian fitting
    protected int xyPartRad; //Radius over which to draw particles in visualisation
    public final int SHOW_RESULTS = -1,
            VERSION = 4;
    protected final double LAMBDA = 650.0, //Wavelength of light
            NUM_AP = 1.4; //Numerical aperture of system
    protected static double colocalThresh = 0.1;
    protected ArrayList<ParticleTrajectory> trajectories = new ArrayList(); //Trajectories of the detected particles
    protected ImagePlus imp; //The active image stack
    protected ImageStack stacks[];
    private long startTime;
    protected DecimalFormat numFormat = new DecimalFormat("0.000");
    protected DecimalFormat intFormat = new DecimalFormat("000");
    String title = "Particle Tracker", ext;
    protected static boolean msdPlot = false, intensPlot = false, trajPlot = false, prevRes = true;
    protected boolean monoChrome;
    private final double SEARCH_SCALE = 1.0;
    private final double TRACK_LENGTH = 5.0;
    private final double TRACK_WIDTH = 4.0;
    public static final float TRACK_OFFSET = 1.0f;
    private static File inputDir = new File("C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences\\TestSequence43\\Input\\C0"),
            c0Dir, c1Dir, outputDir = new File(inputDir.getAbsolutePath()),
            calDir = new File(inputDir.getAbsolutePath());
    private final String delimiter = GenUtils.getDelimiter();
    String parentDir;

//    static {
//        System.loadLibrary("cuda_gauss_tracker"); // Load native library at runtime cudaGaussFitter.dll
//    }
//
//    private native boolean cudaGaussFitter(String folder, String ext, float spatialRes, float sigmaEst, float maxthresh, float fitTol, int startSlice, int endSlice);
//    public static void main(String args[]) {
////        if (imp != null) {
//        Timelapse_Analysis instance = new Timelapse_Analysis();
//        instance.run(null);
////        }
//        System.exit(0);
//    }
    public Timelapse_Analysis(double spatialRes, double timeRes, double trajMaxStep, double chan1MaxThresh, boolean monoChrome, ImagePlus imp, double scale, double minTrajLength) {
        UserVariables.setSpatialRes(spatialRes);
        UserVariables.setTimeRes(timeRes);
        UserVariables.setTrajMaxStep(trajMaxStep);
        UserVariables.setChan1MaxThresh(chan1MaxThresh);
//        Timelapse_Analysis.hystDiff = hystDiff;
        UserVariables.setMinTrajLength(minTrajLength);
        this.monoChrome = monoChrome;
//        this.imp = imp;
//        this.stacks = imp.getStack();
    }

    public Timelapse_Analysis() {
    }

    public Timelapse_Analysis(ImageStack[] stacks) {
        this.stacks = stacks;
    }

    public Timelapse_Analysis(ImagePlus imp, String ext) {
//        this.imp = imp;
//        this.stacks = imp.getImageStack();
        this.ext = ext;
    }

    /**
     * Implements run method from {@link PlugIn}.
     */
    public void run(String arg) {
        Utilities.setLookAndFeel(UserInterface.class);
        title = title + "_v" + VERSION + "." + intFormat.format(Revision.Revision.revisionNumber);
        inputDir = Utilities.getFolder(inputDir, null, true);
        if (inputDir == null) {
            return;
        }
        c0Dir = new File(inputDir.getAbsolutePath() + delimiter + "C0");
        ImagePlus imp = Utils.buildStack(c0Dir);
        stacks = new ImageStack[2];
        stacks[0] = imp.getImageStack();
        this.ext = imp.getTitle();
        c1Dir = new File(inputDir.getAbsolutePath() + delimiter + "C1");
        if (c1Dir.exists()) {
            stacks[1] = (Utils.buildStack(c1Dir)).getImageStack();
        }
//        if (IJ.getInstance() != null) {
//            imp = WindowManager.getCurrentImage();
//        }
        if (showDialog()) {
            analyse();
        }
    }

    public boolean showDialog() {
//        if (imp == null) {
//            Toolkit.getDefaultToolkit().beep();
//            IJ.error("No image stack open.");
//            return false;
//        }
//        stack = imp.getImageStack();
        monoChrome = (stacks[1] == null);
        UserInterface ui = new UserInterface(null, true, title, this);
        ui.setVisible(true);
        return ui.isWasOKed();
    }

    /**
     * Analyses the {@link ImageStack} specified by <code>stack</code>.
     */
    public void analyse() {
        outputDir = Utilities.getFolder(outputDir, "Specify directory for output files...", true);
        parentDir = GenUtils.openResultsDirectory(outputDir + delimiter + title, delimiter);
        if (!monoChrome) {
            calDir = Utilities.getFolder(calDir, "Specify directory containing calibrations...", true);
        }
        if (stacks != null) {
            IJ.register(this.getClass());
            startTime = System.currentTimeMillis();
            //gaussianRadius = 0.139d / spatialRes; // IsoGaussian filter radius set to 139 nm
            calcParticleRadius(UserVariables.getSpatialRes());
            int i, count;
//            int width = stacks[0].getWidth(), height = stacks[0].getHeight();
            if (UserVariables.isGpu()) {
//                cudaFindParticles(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1], monoChrome);
            } else {
                findParticles(SEARCH_SCALE, true, 0, stacks[0].getSize() - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1], monoChrome);
            }
            TextWindow results = new TextWindow(title + " Results", "X\tY\tFrame\tChannel 1 ("
                    + UserVariables.channels[UserVariables.getC1Index()]
                    + ")\tChannel 2 (" + UserVariables.channels[UserVariables.getC2Index()]
                    + ")\tChannel 2 " + '\u03C3' + "x\tChannel 2 " + '\u03C3' + "y\t" + '\u03B8',
                    new String(), 1000, 500);
//            results.append(imp.getTitle() + "\n\n");
            TextWindow resultSummary = new TextWindow(title + " Results Summary",
                    "Particle\tType\t% Colocalisation\tDuration (s)\tDisplacement (" + IJ.micronSymbol
                    + "m)\tVelocity (" + IJ.micronSymbol + "m/s)\tDirectionality\tDiffusion Coefficient ("
                    + IJ.micronSymbol + "m^2/s)" + "\tFractal Dimension"
                    + "\tFluorescence Ratio ("
                    + UserVariables.channels[UserVariables.getC2Index()] + "/"
                    + UserVariables.channels[UserVariables.getC1Index()]
                    + ")\tAngle Spread\tStep Spread\tDC\tCurvature\tC2 Fluor Area\tC2 Fluor Skew",
                    new String(), 1200, 500);
//            resultSummary.append(imp.getTitle() + "\n\n");

            int n = trajectories.size();
            for (i = 0; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (!(traj.getSize() > UserVariables.getMinTrajLength() && ((traj.getType(colocalThresh) == ParticleTrajectory.COLOCAL)
                        || ((traj.getType(colocalThresh) == ParticleTrajectory.NON_COLOCAL) && !UserVariables.isColocal())))) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }
            if (prevRes) {
                ArrayList<Integer> excludeList = previewResults();
                if (excludeList != null) {
                    for (Integer e : excludeList) {
                        trajectories.set(e, null);
                    }
                }
            }
            n = trajectories.size();
            ImageStack maps = mapTrajectories((new RGBStackMerge()).mergeStacks(stacks[0].getWidth(), stacks[0].getHeight(), stacks[0].getSize(), stacks[0], stacks[1], null, true),
                    trajectories, UserVariables.getSpatialRes(), UserVariables.getMinTrajLength(),
                    UserVariables.getTimeRes(), true, 0, trajectories.size() - 1, 1, false);
            for (i = 0, count = 1; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (traj != null) {
                    int type = traj.getType(colocalThresh);
                    if (traj.getSize() > UserVariables.getMinTrajLength() && ((type == ParticleTrajectory.COLOCAL)
                            || ((type == ParticleTrajectory.NON_COLOCAL) && !UserVariables.isColocal()))) {
//                        if (intensPlot) {
//                            plotIntensity(i, count);
//                        }
//                        if (trajPlot) {
//                            plotTrajectory(width, height, i, count);
//                        }
                        printData(i, resultSummary, count);
                        traj.printTrajectory(count, results, numFormat, title);
                        if (type == ParticleTrajectory.COLOCAL) {
                            ImageStack signals[] = extractSignalValues(traj,
                                    (int) Math.round(TRACK_LENGTH / UserVariables.getSpatialRes()),
                                    (int) Math.round(TRACK_WIDTH / UserVariables.getSpatialRes()), TRACK_OFFSET / ((float) UserVariables.getSpatialRes()));
                            if (signals[0].getSize() > 0) {
                                for (int j = 1; j <= signals[0].getSize(); j++) {
                                    IJ.saveAs((new ImagePlus("", signals[0].getProcessor(j))),
                                            "TIF", parentDir + "/C0-" + String.valueOf(count)
                                            + "-" + String.valueOf(j));
                                    IJ.saveAs((new ImagePlus("", signals[1].getProcessor(j))),
                                            "TIF", parentDir + "/C1-" + String.valueOf(count)
                                            + "-" + String.valueOf(j));
                                }
                            }
                        }
                        count++;
                    }
                }
            }
            resultSummary.append("\nAnalysis Time (s): " + numFormat.format((System.currentTimeMillis() - startTime) / 1000.0));
            results.append(toString());
            results.setVisible(true);
            resultSummary.setVisible(true);
            if (maps != null) {
                (new ImagePlus("Trajectory Maps", maps)).show();
                IJ.saveAs((new ImagePlus("", maps)), "TIF", parentDir + "/trajectories.tif");
            }
        }
        printParams(parentDir);
    }

    /**
     * Median filter and IsoGaussian filter the image specified by
     * <code>processor</code>.
     *
     * @param processor the image to be pre-processed.
     */
    public FloatProcessor preProcess(ImageProcessor processor) {
        if (processor == null) {
            return null;
        }
        FloatProcessor fp;
        if (UserVariables.isPreProcess()) {
            TypeConverter tc = new TypeConverter(processor, false);
            fp = (FloatProcessor) tc.convertToFloat(null);
            (new GaussianBlur()).blur(fp, xySigEst);
        } else {
            TypeConverter tc = new TypeConverter(processor, false);
            fp = (FloatProcessor) tc.convertToFloat(null);
        }
        return fp;
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice, double fitTol, ImageStack channel1, ImageStack channel2, boolean monoChrome) {
        if (channel1 == null) {
            return null;
        }
        int i, noOfImages = channel1.getSize(), width = channel1.getWidth(), height = channel1.getHeight(),
                arraySize = endSlice - startSlice + 1;
        int fitRad = (int) Math.ceil(xyPartRad * 4.0 / 3.0);
        int c1X, c1Y, pSize = 2 * fitRad + 1;
        int c2Points[][];
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
//        ImageStack detect_output = new ImageStack(stack.getWidth(), stack.getHeight());
//        ImageStack maxima = new ImageStack(stack.getWidth(), stack.getHeight());
//        ImageStack input_output = new ImageStack(stack.getWidth(), stack.getHeight());
        ProgressDialog progress = new ProgressDialog(null, "Finding Particles...", false, title, false);
        progress.setVisible(true);
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
//            ByteProcessor oslice = new ByteProcessor(detect_output.getWidth(), detect_output.getHeight());
            IJ.freeMemory();
            progress.updateProgress(i - startSlice, arraySize);
            FloatProcessor chan1Proc = preProcess(channel1.getProcessor(i + 1));
            FloatProcessor chan2Proc = !monoChrome ? preProcess(channel2.getProcessor(i + 1)) : null;
            // TODO Maybe express threshold as a percentage? See Ponti et al., 2003
            double c1Threshold = Utils.getPercentileThresh(chan1Proc, UserVariables.getChan1MaxThresh());
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, UserVariables.FOREGROUND, chan1Proc, c1Threshold, true);
//            maxima.addSlice(thisC1Max);
            double c2Threshold = Utils.getPercentileThresh(chan2Proc, UserVariables.getChan2MaxThresh());
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == UserVariables.FOREGROUND) {
                        IsoGaussian c2Gaussian = null;
                        c2Points = Utils.searchNeighbourhood(c1X, c1Y,
                                (int) Math.round(fitRad * searchScale),
                                (int) Math.round(c2Threshold),
                                (ImageProcessor) chan2Proc);
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                        MultiGaussFitter c1Fitter = new MultiGaussFitter(UserVariables.getnMax(), fitRad, pSize);
                        c1Fitter.fit(pixValues, xySigEst);
                        ArrayList<IsoGaussian> c1Fits = c1Fitter.getFits(spatialRes, c1X - fitRad, c1Y - fitRad, c1Threshold, fitTol);

                        if (c2Points != null) {
                            Utils.extractValues(xCoords, yCoords, pixValues,
                                    c2Points[0][0], c2Points[0][1], chan2Proc);
                            MultiGaussFitter c2Fitter = new MultiGaussFitter(1, fitRad, pSize);
                            c2Fitter.fit(pixValues, xySigEst);
                            ArrayList<IsoGaussian> c2Fits = c2Fitter.getFits(spatialRes, c2Points[0][0] - fitRad * searchScale, c2Points[0][1] - fitRad * searchScale, c2Threshold, UserVariables.getCurveFitTol());
                            if (c2Fits != null && c2Fits.size() > 0) {
                                c2Gaussian = c2Fits.get(0);
                            }
                        }

                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        if (c1Fits != null) {
                            for (IsoGaussian c1Fit : c1Fits) {
                                particles.addDetection(i - startSlice, new Particle(i - startSlice, c1Fit, c2Gaussian, null, -1));
//                                Utils.draw2DGaussian(oslice, c1Fit, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false, false);
//                                Utils.draw2DGaussian(chan1Proc, c1Fit, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), false, true);
                            }
                        }
                    }
                }
            }

//            detect_output.addSlice(oslice);
//            input_output.addSlice(chan1Proc.duplicate());
        }
        progress.dispose();
//        IJ.saveAs(new ImagePlus("", detect_output), "TIF", parentDir + "/all_detections.tif");
//        IJ.saveAs(new ImagePlus("", maxima), "TIF", parentDir + "/all_maxima.tif");
//        IJ.saveAs(new ImagePlus("", input_output), "TIF", parentDir + "/input_output.tif");
        if (update) {
            updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(), spatialRes, true);
        }
        return particles;
    }

//    public ParticleArray cudaFindParticles(double searchScale, boolean update, int startSlice, int endSlice, double fitTol, ImageStack channel1, ImageStack channel2, boolean monoChrome) {
//        if (!cudaGaussFitter(c0Dir.getAbsolutePath(), ext, (float) UserVariables.getSpatialRes() * 1000.0f, (float) xySigEst, (float) UserVariables.getChan1MaxThresh(), (float) UserVariables.getCurveFitTol(), startSlice, endSlice)) {
//            IJ.log("CUDA Error");
//            return null;
//        }
//        File cudaFile = new File(c0Dir + delimiter + "cudadata.txt");
//        File fileList[] = {cudaFile};
//        ArrayList<double[]>[] cudaData = GenUtils.readData(5, fileList, delimiter);
//        int arraySize = endSlice - startSlice + 1;
//        int fitRad = (int) Math.ceil(xyPartRad * 4.0 / 3.0);
////        int pSize = 2 * fitRad + 1;
//        int c2Points[][];
////        double[] xCoords = new double[pSize];
////        double[] yCoords = new double[pSize];
////        double[][] pixValues = new double[pSize][pSize];
//        double spatialRes = UserVariables.getSpatialRes();
//        ParticleArray particles = new ParticleArray(arraySize);
//        double c2Thresholds[] = new double[channel2.getSize()];
//        ImageStack procChannel2 = new ImageStack(channel2.getWidth(), channel2.getHeight());
//        for (int i = 0; i < channel2.getSize(); i++) {
//            procChannel2.addSlice(preProcess(channel2.getProcessor(i + 1).duplicate()));
//            c2Thresholds[i] = Utils.getPercentileThresh(procChannel2.getProcessor(i + 1),
//                    UserVariables.getChan2MaxThresh());
//        }
//        for (int f = 0; f < fileList.length; f++) {
//            ProgressDialog progress = new ProgressDialog(null,
//                    "Reading data for file " + f + " of " + fileList.length + "...",
//                    false, true, title);
//            progress.setVisible(true);
//            int size = cudaData[f].size();
//            for (int row = 0; row < size; row++) {
//                progress.updateProgress(row, size);
////                IJ.freeMemory();
//                int t = (int) Math.round(cudaData[f].get(row)[0]);
//                double x = cudaData[f].get(row)[1];
//                double y = cudaData[f].get(row)[2];
//                double mag = cudaData[f].get(row)[3];
//                double fit = cudaData[f].get(row)[4];
////                FloatProcessor chan2Proc = !monoChrome ? preProcess(channel2.getProcessor(t + 1)) : null;
//                IsoGaussian c1Gaussian = new IsoGaussian(x, y, mag, xySigEst, xySigEst, fit);
//                int c1X = (int) Math.round(x / UserVariables.getSpatialRes());
//                int c1Y = (int) Math.round(y / UserVariables.getSpatialRes());
//                IsoGaussian c2Gaussian = null;
//                c2Points = Utils.searchNeighbourhood(c1X, c1Y,
//                        (int) Math.round(fitRad * searchScale),
//                        (int) Math.round(c2Thresholds[t]),
//                        procChannel2.getProcessor(t + 1));
//                if (c2Points != null) {
//                    int xC2 = c2Points[0][0] - (int) Math.round(fitRad * searchScale);
//                    int yC2 = c2Points[0][1] - (int) Math.round(fitRad * searchScale);
//                    c2Gaussian = new IsoGaussian(xC2 * spatialRes, yC2 * spatialRes,
//                            procChannel2.getProcessor(t + 1).getPixelValue(xC2, yC2), xySigEst, xySigEst, 0.0);
////                    Utils.extractValues(xCoords, yCoords, pixValues,
////                            c2Points[0][0], c2Points[0][1], chan2Proc);
////                    MultiGaussFitter c2Fitter = new MultiGaussFitter(1, fitRad, pSize);
////                    c2Fitter.fit(pixValues, xySigEst);
////                    ArrayList<IsoGaussian> c2Fits = c2Fitter.getFits(spatialRes, c2Points[0][0] - fitRad * searchScale, c2Points[0][1] - fitRad * searchScale, c2Threshold, UserVariables.getCurveFitTol());
////                    if (c2Fits != null && c2Fits.size() > 0) {
////                        c2Gaussian = c2Fits.get(0);
////                    }
//                }
//                particles.addDetection(t - startSlice, new Particle(t, c1Gaussian, c2Gaussian, null, -1));
//            }
//            progress.dispose();
//        }
//        if (update) {
//            updateTrajectories(particles, UserVariables.getTimeRes(), UserVariables.getTrajMaxStep(), spatialRes, true);
//        }
//        return particles;
//    }
    public void updateTrajectories(ParticleArray objects, double timeRes, double trajMaxStep, double spatialRes, boolean projectPos) {
        if (objects == null) {
            return;
        }
        int depth = objects.getDepth();
        ParticleTrajectory traj = null;
        double x, y, score, minScore;
        ProgressDialog progress = new ProgressDialog(null, "Building Trajectories...", false, title, false);
        progress.setVisible(true);
        for (int m = 0; m < depth; m++) {
            progress.updateProgress(m, depth);
            for (int k = m; (k < depth) && (((k - m)) < trajMaxStep); k++) {
                int size = trajectories.size();
                ArrayList<Particle> detections = objects.getLevel(k);
                for (int j = 0; j < detections.size(); j++) {
                    Particle currentParticle = detections.get(j);
                    if (currentParticle != null) {
                        IsoGaussian ch1G = currentParticle.getC1Gaussian();
                        /*
                         * If no trajectories have yet been built, start a new
                         * one:
                         */
                        if (k == m) {
                            traj = new ParticleTrajectory(timeRes, spatialRes);
                            /*
                             * Particles need to be cloned as they are set to
                             * null once inserted into trajectories.
                             */
                            traj.addPoint((Particle) currentParticle.clone());
                            trajectories.add(traj);
                            /*
                             * Otherwise, determine whether the current particle
                             * belongs to a pre-existing trajectory:
                             */
                        } else {
                            int i, minScoreIndex;
                            for (minScoreIndex = -1, minScore = Double.MAX_VALUE, i = 0; i < size; i++) {
                                traj = (ParticleTrajectory) trajectories.get(i);
                                Particle last = traj.getEnd();
                                if ((last != null) && (last.getTimePoint() == m) && k != m) {
                                    /*
                                     * Evaluate the probability that the current
                                     * particle belongs to the current
                                     * trajectory, based on the particle's
                                     * distance from last point on the current
                                     * trajectory and the number of frames
                                     * between the current particle, the last
                                     * point of the current trajectory and
                                     * differences in respective intensity
                                     * levels:
                                     */
                                    x = ch1G.getX();
                                    y = ch1G.getY();
                                    double vector1[] = {x, y, currentParticle.getTimePoint(),
                                        ch1G.getMagnitude()};
                                    double vector2[] = {last.getX(), last.getY(),
                                        last.getTimePoint(), last.getC1Intensity()};
                                    score = Utils.calcEuclidDist(vector1, vector2);
                                    if (projectPos) {
                                        double vector3[] = {x, y};
                                        double vector4[] = {last.getX() + traj.getXVelocity(),
                                            last.getY() + traj.getYVelocity()};
                                        score += Utils.calcEuclidDist(vector3, vector4);
                                    }
                                    if (score < minScore) {
                                        minScore = score;
                                        minScoreIndex = i;
                                    }
                                }
                            }
                            /*
                             * If an acceptably low score has been evaluated,
                             * the particle is temporarily assigned to the
                             * "winning" trajectory:
                             */
                            if (traj != null) {
                                if (minScoreIndex > -1) {
                                    traj = (ParticleTrajectory) trajectories.get(minScoreIndex);
                                }
                                if ((minScore < trajMaxStep) && (minScore < traj.getTempScore())) {
                                    traj.addTempPoint((Particle) currentParticle.clone(), minScore, j, k);
                                }
                            }
                        }
                    }
                }
            }
            for (ParticleTrajectory trajectorie : trajectories) {
                traj = (ParticleTrajectory) trajectorie;
                Particle temp = traj.getTemp();
                if (temp != null) {
                    int row = traj.getTempRow();
                    int col = traj.getTempColumn();
                    if (col <= m + 1) {
                        traj.checkDetections(temp, 0.0, 0.0);
                        objects.nullifyDetection(col, row);
                    }
                }
            }
        }
        progress.dispose();
    }

    /**
     * Produces a {@link Plot} of the trajectory specified by
     * <code>particleNumber</code>.
     *
     * @param width the width of the images from which the trajectory was
     * extracted
     * @param height the height of the images from which the trajectory was
     * extracted
     * @param particleNumber the trajectory index
     */
    public boolean plotTrajectory(int width, int height, int particleNumber, int label) {
        if (width * height == 0 || trajectories.size() < 1) {
            return false;
        }
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        if (traj == null) {
            return false;
        }
        Particle current = traj.getEnd();
        int size = traj.getSize();
        double xValues[] = new double[size];
        double yValues[] = new double[size];
        width *= UserVariables.getSpatialRes();
        height *= UserVariables.getSpatialRes();

        for (int i = size - 1; i >= 0; i--) {
            xValues[i] = (double) current.getX() / width;
            /*
             * Y-coordinates are inverted so trajectory is displayed as per the
             * image
             */
            yValues[i] = -(double) current.getY() / height;
            current = current.getLink();
        }

        Plot plot = new Plot("Particle " + label + " Trajectory",
                "X", "Y", xValues, yValues,
                (Plot.X_TICKS + Plot.Y_TICKS + Plot.X_NUMBERS + Plot.Y_NUMBERS));
        plot.setLimits(0, 1.0, -1.0, 0);
        plot.setLineWidth(2);
        plot.setColor(Color.BLUE);
        plot.draw();
        plot.show();

        return true;
    }

    /**
     * Outputs velocity and directionality data on the particle specified by
     * <code>particleNumber</code>. Directionality ( <code>D</code>) is
     * calculated according to: <br> <br>
     * <code>D = 1 / (1 + &lambda<sub>1</sub>
     * &lambda<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where <code>&lambda<sub>1</sub></code> and
     * <code>&lambda<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and      <code>&lambda<sub>1</sub> < &lambda<sub>2</sub></code>.
     *
     */
    public boolean printData(int particleNumber, TextWindow output, int label) {
        if (trajectories.size() < 1) {
            return false;
        }
        DecimalFormat decFormat = new DecimalFormat("0.000");
        DecimalFormat msdFormat = new DecimalFormat("0.000000");
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        if (traj == null) {
            return false;
        }
//        traj.smooth();
        traj.calcMSD(label, -1, msdPlot);
        traj.calcAngleSpread();
        traj.calcStepSpread();
        traj.calcDirectionality();
        double displacement = traj.getDisplacement(traj.getEnd(), traj.getSize());
        double duration = traj.getDuration();
        int type = traj.getType(colocalThresh);
        String trajType = null;
        switch (type) {
            case ParticleTrajectory.COLOCAL:
                trajType = "Colocalised";
                break;
            case ParticleTrajectory.NON_COLOCAL:
                trajType = "Non-Colocalised";
                break;
            case ParticleTrajectory.UNKNOWN:
                trajType = "Unknown";
        }
        output.append(label + "\t" + trajType + "\t"
                + decFormat.format(traj.getDualScore() * 100.0 / traj.getSize()) + "\t"
                + decFormat.format(duration) + "\t"
                + decFormat.format(displacement)
                + "\t" + decFormat.format(displacement / duration) + "\t"
                + decFormat.format(traj.getDirectionality()) + "\t"
                + msdFormat.format(traj.getDiffCoeff()) + "\t"
                + decFormat.format(traj.getBoxCountFD()) + "\t"
                + decFormat.format(traj.getFluorRatio()) + "\t"
                + decFormat.format(traj.getAngleSpread()) + "\t"
                + decFormat.format(traj.getStepSpread()) + "\t"
                + decFormat.format(traj.getLogDC()) + "\t"
                + decFormat.format(traj.getMeanKappa()) + "\t");
        return true;
    }

//    /**
//     * Produces a {@link Plot} of normalised intensity in the red and green
//     * channels, together with a ratio of red:green intensity, for the
//     * trajectory denoted by <code>particleNumber</code>.
//     */
//    public boolean plotIntensity(int particleNumber, int label) {
//        if (trajectories.size() < 1) {
//            return false;
//        }
//        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
//        Particle current = traj.getEnd();
//        int size = traj.getSize();
//        double xvalues[] = new double[size];
//        double redValues[] = new double[size];
//        double greenValues[] = new double[size];
//        double ratios[] = new double[size];
//        double temp, maxVal = -Double.MAX_VALUE, minVal = Double.MAX_VALUE;
//
//        if (traj.getType(colocalThresh) != ParticleTrajectory.COLOCAL) {
//            return false;
//        }
//
//        for (int i = size - 1; i >= 0; i--) {
//            xvalues[i] = current.getTimePoint() * UserVariables.getTimeRes();
//            redValues[i] = current.getC1Intensity() / 255.0d;
//            greenValues[i] = current.getC2Intensity() / 255.0d;
//            ratios[i] = greenValues[i] / redValues[i];
//            if (ratios[i] > maxVal) {
//                maxVal = ratios[i];
//            }
//            temp = Math.min(redValues[i], greenValues[i]);
//            if (temp < minVal) {
//                minVal = temp;
//            }
//            current = current.getLink();
//        }
//
//        Plot plot = new Plot("Particle " + label + " Intensities",
//                "Time (s)", "Normalised Intensity", xvalues, ratios,
//                (Plot.X_TICKS + Plot.Y_TICKS + Plot.X_NUMBERS + Plot.Y_NUMBERS));
//        plot.changeFont(new Font("Serif", Font.BOLD, 14));
//        plot.setLimits(0.0, (stacks.getSize() - 1.0) * UserVariables.getTimeRes(), minVal, maxVal);
//        plot.setLineWidth(2);
//        plot.setColor(Color.BLUE);
//        plot.draw();
//        plot.setColor(Color.RED);
//        plot.addPoints(xvalues, redValues, Plot.LINE);
//        plot.setColor(Color.GREEN);
//        plot.addPoints(xvalues, greenValues, Plot.LINE);
//        plot.show();
//
//        return true;
//    }
    /**
     * Constructed trajectories are drawn onto the original image sequence and
     * displayed as a stack sequence.
     */
    public ImageStack mapTrajectories(ImageStack stack, ArrayList<ParticleTrajectory> trajectories, double spatialRes, double minTrajLength, double timeRes, boolean tracks, int startT, int endT, int index, boolean preview) {
        if (stack == null) {
            return null;
        }
        int radius = xyPartRad;
        int i, j, width = stack.getWidth(), height = stack.getHeight(),
                type, frames = stack.getSize();
        double lastX, lastY;
        ImageStack outputStack = new ImageStack(width, height);
        Particle current;
        ParticleTrajectory traj;
        int length, n = trajectories.size();
        ImageProcessor processor;
//        Rectangle bounds;
//        ByteProcessor trajImage;
        int lastTP;
//        double ptScale = ParticleTrajectory.scale;

        if (n < 1) {
            return null;
        }

        for (i = 0; i < frames; i++) {
            if (monoChrome) {
                processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte());
            } else {
                processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToRGB());
            }
            processor.setInterpolationMethod(ImageProcessor.BICUBIC);
            processor.setInterpolate(true);
            outputStack.addSlice("" + i, processor.resize(width, height));
        }
        Random r = new Random();
        int tLength = (int) Math.round(TRACK_LENGTH / UserVariables.getSpatialRes());
        ProgressDialog progress = new ProgressDialog(null, "Mapping Output...", false, title, false);
        progress.setVisible(true);
        for (i = startT; i <= endT && i < n; i++) {
            progress.updateProgress(i, n);
            traj = (ParticleTrajectory) (trajectories.get(i));
            if (traj != null) {
                Color thiscolor = new Color(r.nextInt(256), r.nextInt(256), r.nextInt(256));
                length = traj.getSize();
                type = traj.getType(colocalThresh);
//                bounds = traj.getBounds();
//                trajImage = new ByteProcessor(bounds.width, bounds.height);
//                trajImage.setColor(Color.white);
//                trajImage.fill();
//                trajImage.setColor(Color.black);
                if (length > minTrajLength && ((type == ParticleTrajectory.COLOCAL)
                        || ((type == ParticleTrajectory.NON_COLOCAL) && !UserVariables.isColocal()))) {
                    current = traj.getEnd();
                    lastX = current.getX();
                    lastY = current.getY();
                    lastTP = current.getTimePoint();
                    current = current.getLink();
                    while (current != null) {
                        for (j = frames - 1; j >= lastTP; j--) {
                            processor = outputStack.getProcessor(j + 1);
                            if (!monoChrome && !preview) {
                                processor.setColor(thiscolor);
                            } else {
                                processor.setColor(Color.white);
                            }
                            if (j - 1 < lastTP) {
                                markParticle(processor, (int) Math.round(lastX / spatialRes) - radius,
                                        (int) Math.round(lastY / spatialRes) - radius, radius, true, "" + index);
                            }
                            if (tracks && j <= lastTP + tLength / timeRes) {
                                int x = (int) (Math.round(current.getX() / spatialRes));
                                int y = (int) (Math.round(current.getY() / spatialRes));
                                processor.drawLine(x, y, (int) Math.round(lastX / spatialRes),
                                        (int) Math.round(lastY / spatialRes));
                            }
                        }
//                        if (tracks) {
//                            trajImage.drawLine((int) Math.round(ptScale * current.getX() - traj.getBounds().x),
//                                    (int) Math.round(ptScale * current.getY() - traj.getBounds().y),
//                                    (int) Math.round(ptScale * lastX - traj.getBounds().x),
//                                    (int) Math.round(ptScale * lastY - traj.getBounds().y));
//                        }
                        lastX = current.getX();
                        lastY = current.getY();
                        lastTP = current.getTimePoint();
                        current = current.getLink();
                    }
                    processor = outputStack.getProcessor(lastTP + 1);
                    if (!monoChrome && !preview) {
                        processor.setColor(thiscolor);
                    } else {
                        processor.setColor(Color.white);
                    }
                    markParticle(processor, (int) Math.round(lastX / spatialRes) - radius,
                            (int) Math.round(lastY / spatialRes) - radius, radius, true, "" + index);
                    index++;
                }
            }
        }
        progress.dispose();
        return outputStack;
    }

    void markParticle(ImageProcessor processor, int x, int y, int radius, boolean string, String label) {
        processor.drawRect(x, y, 2 * radius + 1, 2 * radius + 1);
//        processor.drawOval(x, y, 2 * radius, 2 * radius);
        if (string) {
            processor.drawString(label, x, y);
        }
    }

    public void calcParticleRadius(double spatialRes) {
        double airyRad = 1.22 * LAMBDA / (2.0 * NUM_AP); //Airy radius
//        xyPartRad = (int) Math.ceil(2.0*airyRad / (spatialRes * 1000.0));
        xySigEst = airyRad / (2.0 * spatialRes * 1000.0);
        xyPartRad = (int) Math.ceil(3.0 * xySigEst);
    }

    public ArrayList<Integer> previewResults() {
        if (trajectories.size() < 1) {
            return null;
        }
        ResultsPreviewInterface previewUI = new ResultsPreviewInterface(IJ.getInstance(), true, title, this);
        previewUI.setVisible(true);
        if (previewUI.isWasOKed()) {
            return previewUI.getRemoveList();
        } else {
            return null;
        }
    }

    public Color getDrawColor(int key) {
        switch (key) {
            case UserVariables.RED:
                return Color.red;
            case UserVariables.GREEN:
                return Color.green;
            case UserVariables.BLUE:
                return Color.blue;
            default:
                return Color.white;
        }
    }

    public ArrayList<ParticleTrajectory> getTrajectories() {
        return trajectories;
    }

    public ImageStack[] getStacks() {
        return stacks;
    }

    public boolean isMonoChrome() {
        return monoChrome;
    }

    /**
     * Procedure for obtaining Goshtasby coefficients - see Goshtasby (1988),
     * IEEE Transactions on Geoscience and Remote Sensing, 26:60-64
     *
     * 1. "Track" beads in two channels to provide list of coordinates 2. Input
     * list of coordinates for each channel to MATLAB 3. Run goshtasby.m to
     * obtain translation coefficients for both x and y, mapping 'green'
     * coordinates onto 'red'. 4. Export x-coeffs, y-coeffs and original green
     * channel coords (assuming virus is in red).
     *
     * @param ptraj
     * @param signalLength
     * @param signalWidth
     * @return
     */
    ImageStack[] extractSignalValues(ParticleTrajectory ptraj, int signalLength, int signalWidth, float offset) {
        TextReader reader = new TextReader();
        ImageProcessor xcoeffs = reader.open(calDir + delimiter + "xcoeffs.txt");
        ImageProcessor ycoeffs = reader.open(calDir + delimiter + "ycoeffs.txt");
        ImageProcessor coords = reader.open(calDir + delimiter + "coords.txt");
        Particle sigStartP = ptraj.getEnd();
        if (signalWidth % 2 == 0) {
            signalWidth++;
        }
        int size = (int) Math.round(Math.min(signalLength * 1.2, ptraj.getSize()));
        int iterations = 1 + ptraj.getSize() - size;
        float xSigPoints[] = new float[size + 1];
        float ySigPoints[] = new float[size + 1];
        float xVirPoints[] = new float[size + 1];
        float yVirPoints[] = new float[size + 1];
        ImageProcessor[] sigTemps = new ImageProcessor[iterations];
        ImageProcessor[] virTemps = new ImageProcessor[iterations];
        for (int i = 0; i < iterations; i++) {
            Particle current = sigStartP;
            for (int index = 1; index <= size; index++) {
                double xg = goshtasbyEval(xcoeffs, coords, current.getC1Gaussian().getX(), current.getC1Gaussian().getY());
                double yg = goshtasbyEval(ycoeffs, coords, current.getC1Gaussian().getX(), current.getC1Gaussian().getY());
                xSigPoints[index] = (float) (xg / UserVariables.getSpatialRes());
                ySigPoints[index] = (float) (yg / UserVariables.getSpatialRes());
                xVirPoints[index] = (float) (current.getC1Gaussian().getX() / UserVariables.getSpatialRes());
                yVirPoints[index] = (float) (current.getC1Gaussian().getY() / UserVariables.getSpatialRes());
                current = current.getLink();
            }
            extendSignalArea(xSigPoints, ySigPoints, offset, size + 1, 1);
            extendSignalArea(xVirPoints, yVirPoints, offset, size + 1, 1);
            PolygonRoi sigProi = new PolygonRoi(xSigPoints, ySigPoints, size + 1, Roi.POLYLINE);
            PolygonRoi virProi = new PolygonRoi(xVirPoints, yVirPoints, size + 1, Roi.POLYLINE);
            Straightener straightener = new Straightener();
            ImagePlus sigImp = new ImagePlus("", stacks[1].getProcessor(sigStartP.getTimePoint()));
            ImagePlus virImp = new ImagePlus("", stacks[0].getProcessor(sigStartP.getTimePoint()));
            sigImp.setRoi(sigProi);
            virImp.setRoi(virProi);
            sigTemps[i] = straightener.straighten(sigImp, sigProi, signalWidth);
            virTemps[i] = straightener.straighten(virImp, virProi, signalWidth);
            sigStartP = sigStartP.getLink();
        }
        int xc = (int) Math.ceil(0.75 * offset);
        int yc = (signalWidth - 1) / 2;
        int outputWidth = (int) Math.round(signalLength + offset);
        ImageStack output[] = new ImageStack[2];
        output[0] = new ImageStack(outputWidth, signalWidth);
        output[1] = new ImageStack(outputWidth, signalWidth);

        for (int j = 0; j < iterations; j++) {
            if (sigTemps[j] != null && sigTemps[j].getWidth() >= outputWidth) {
                ImageStack virStack = new ImageStack(virTemps[j].getWidth(), virTemps[j].getHeight());
                virStack.addSlice(virTemps[j]);
                ParticleArray particles = findParticles(0.0, false, 0, 0, UserVariables.getCurveFitTol(), virStack, null, true);
                if (!particles.getLevel(0).isEmpty()) {
                    virTemps[j].setInterpolate(true);
                    virTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    sigTemps[j].setInterpolate(true);
                    sigTemps[j].setInterpolationMethod(ImageProcessor.BICUBIC);
                    double xinc = particles.getLevel(0).get(0).getC1Gaussian().getX() / UserVariables.getSpatialRes() - xc;
                    double yinc = particles.getLevel(0).get(0).getC1Gaussian().getY() / UserVariables.getSpatialRes() - yc;
                    virTemps[j].translate(-xinc, -yinc);
                    sigTemps[j].translate(-xinc, -yinc);
                    FloatProcessor sigSlice = new FloatProcessor(outputWidth, signalWidth);
                    FloatBlitter sigBlitter = new FloatBlitter(sigSlice);
                    sigBlitter.copyBits(sigTemps[j], 0, 0, Blitter.COPY);
                    output[1].addSlice(sigSlice);
                    FloatProcessor virSlice = new FloatProcessor(outputWidth, signalWidth);
                    FloatBlitter virBlitter = new FloatBlitter(virSlice);
                    virBlitter.copyBits(virTemps[j], 0, 0, Blitter.COPY);
                    output[0].addSlice(virSlice);
                }
            }
        }
        return output;
    }

    void extendSignalArea(float[] xpoints, float[] ypoints, float dist, int n, int window) {
        float xdiff = 0.0f, ydiff = 0.0f;
        for (int i = 1; i <= window; i++) {
            xdiff += xpoints[i] - xpoints[i + 1];
            ydiff += ypoints[i] - ypoints[i + 1];
        }
        float ratio = Math.abs(ydiff / xdiff);
        float newX = dist / (float) Math.sqrt(1.0f + (float) Math.pow(ratio, 2.0f));
        float newY = newX * ratio;
        if (xdiff < 0.0) {
            xpoints[0] = xpoints[1] - newX;
        } else {
            xpoints[0] = xpoints[1] + newX;
        }
        if (ydiff < 0.0) {
            ypoints[0] = ypoints[1] - newY;
        } else {
            ypoints[0] = ypoints[1] + newY;
        }
    }

    public int getXyPartRad() {
        return xyPartRad;
    }

    double goshtasbyEval(ImageProcessor coeffs, ImageProcessor coords, double x, double y) {
        int l = coeffs.getHeight();
        double sum = 0.0;
        for (int i = 3; i < l; i++) {
            double r = Math.pow((x - coords.getPixelValue(0, i - 3)), 2.0) + Math.pow((y - coords.getPixelValue(1, i - 3)), 2.0);
            if (r > 0.0) {
                double R = r * Math.log(r);
                sum = sum + coeffs.getPixelValue(0, i) * R;
            }
        }
        return coeffs.getPixelValue(0, 0) + coeffs.getPixelValue(0, 1) * x + coeffs.getPixelValue(0, 2) * y + sum;
    }

    void printParams(String dir) {
        File paramFile;
        PrintWriter paramStream;
        try {
            paramFile = new File(dir + delimiter + "params.csv");
            paramStream = new PrintWriter(new FileOutputStream(paramFile));
        } catch (FileNotFoundException e) {
            System.out.println("Error: Failed to create parameter file.\n");
            System.out.println(e.toString());
            return;
        }
        paramStream.println(title);
        paramStream.println(Utilities.getDate("dd/MM/yyyy HH:mm:ss"));
        paramStream.println();
        paramStream.println(UserInterface.getChannel1LabelText() + "," + UserVariables.getC1Index());
        paramStream.println(UserInterface.getChannel2LabelText() + "," + UserVariables.getC2Index());
        paramStream.println(UserInterface.getSpatResLabelText() + "," + UserVariables.getSpatialRes());
        paramStream.println(UserInterface.getFpsLabelText() + "," + UserVariables.getTimeRes());
        paramStream.println(UserInterface.getMinTrajLengthLabelText() + "," + UserVariables.getMinTrajLength());
        paramStream.println(UserInterface.getMaxLinkDistLabelText() + "," + UserVariables.getTrajMaxStep());
        paramStream.println(UserInterface.getChan1MaxThreshLabelText() + "," + UserVariables.getChan1MaxThresh());
        paramStream.println(UserInterface.getChan2MaxThreshLabelText() + "," + UserVariables.getChan2MaxThresh());
        paramStream.println(UserInterface.getCurveFitTolLabelText() + "," + UserVariables.getCurveFitTol());
        paramStream.println(UserInterface.getnMaxLabelText() + "," + UserVariables.getnMax());
        paramStream.println(UserInterface.getColocalToggleText() + "," + UserVariables.isColocal());
        paramStream.println(UserInterface.getPreprocessToggleText() + "," + UserVariables.isPreProcess());
        paramStream.println(UserInterface.getGpuToggleText() + "," + UserVariables.isGpu());
        paramStream.close();
    }

    public static File getDirectory() {
        return outputDir;
    }

}
