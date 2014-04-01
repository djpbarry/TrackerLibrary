package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.Utils;
import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.*;
import ij.text.TextWindow;
import java.awt.*;
import java.io.File;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;
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

    protected static double timeRes = 2.0d, //Time resolution in s/frame;
            virusDiameter = 350.0,
            chan1MaxThresh = 100.0, //Threshold value for local maxima in channel 1
            chan2MaxThresh = 0.0, //Threshold value for local maxima in channel 2
            //gaussianRadius, //Gaussian filter radius
            trajMaxStep = 2.5, //Tolerance used in evaluating likely trajectories
            //TODO Consider performing a series of analyses with incresing values for trajMaxStep. This may allow discrimination between different modes of motion by identifying slowly-moving particles first (and subsequently removing them from the data set) and then proceeding to more complex trajectories.
            minTrajLength = 0.0, //Minimum trajectory length output in results
            hystDiff = 1.25;
    protected double scale, //Scale factor for visualisation
            xySigEst; //Initial estimate of standard deviation for IsoGaussian fitting
    protected int xyPartRad; //Radius over which to draw particles in visualisation
    public final static int FOREGROUND = 255, //Integer value of foreground pixels
            SHOW_RESULTS = -1,
            VERSION = 3;
    protected final static double VIS_SIZE = 750.0,
            LAMBDA = 650.0, //Wavelength of light
            NUM_AP = 1.4; //Numerical aperture of system
    protected static double c1CurveFitTol = 0.8d, //Tolerance used in determining fit of IsoGaussian curves
            c2CurveFitTol = 0.8d,
            colocalThresh = 0.1;
//    protected static double c1SigmaTol = 3.0, c2SigmaTol = 3.0;
    protected ArrayList<ParticleTrajectory> trajectories = new ArrayList(); //Trajectories of the detected particles
    protected ImagePlus imp; //The active image stack
    protected ImageStack stack;
    private long startTime;
    protected DecimalFormat numFormat = new DecimalFormat("0.000");
    protected DecimalFormat intFormat = new DecimalFormat("000");
    String title = "Particle Tracker";
    protected static boolean colocal = false, msdPlot = false, intensPlot = false,
            preProcess = true, trajPlot = false, prevRes = false;
    protected Co_Localise colocaliser;
    protected boolean monoChrome;
//    private double noiseTol = 0.2;

    public static void main(String args[]) {
        File image = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences"), null);
        ImageStack stack = Utils.buildStack(image);
        ImagePlus imp = new ImagePlus("Stack", stack);
        Timelapse_Analysis instance = new Timelapse_Analysis(imp);
        if (instance.showDialog()) {
            instance.analyse();
        }
    }

    public Timelapse_Analysis(double spatialRes, double timeRes, double trajMaxStep,
            double chan1MaxThresh, double hystDiff, boolean monoChrome, ImagePlus imp, double scale, double minTrajLength) {
        UserVariables.setSpatialRes(spatialRes);
        Timelapse_Analysis.timeRes = timeRes;
        Timelapse_Analysis.trajMaxStep = trajMaxStep;
        Timelapse_Analysis.chan1MaxThresh = chan1MaxThresh;
        Timelapse_Analysis.hystDiff = hystDiff;
        Timelapse_Analysis.minTrajLength = minTrajLength;
        this.monoChrome = monoChrome;
        this.imp = imp;
        this.stack = imp.getStack();
        this.scale = scale;
        colocaliser = new Co_Localise(imp);
        colocaliser.setChannel1(0);
        colocaliser.setChannel2(2);
    }

    public Timelapse_Analysis() {
    }

    public Timelapse_Analysis(ImageStack stack) {
        this.stack = stack;
    }

    public Timelapse_Analysis(ImagePlus imp) {
        this.imp = imp;
        this.stack = imp.getImageStack();
    }

    /**
     * Implements run method from {@link PlugIn}.
     */
    public void run(String arg) {
        title = title + "_v" + VERSION + "." + intFormat.format(Revision.Revision.revisionNumber);
        imp = WindowManager.getCurrentImage();
        if (showDialog()) {
            analyse();
        }
    }

    public boolean showDialog() {
        if (imp == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image stack open.");
            return false;
        }
        colocaliser = new Co_Localise(imp);
        stack = imp.getImageStack();
        monoChrome = !(stack.getProcessor(1).getNChannels() > 1);
        UserInterface ui = new UserInterface(null, true, title, this);
        ui.setVisible(true);
        if (ui.isWasOKed()) {
            return true;
        } else {
            return false;
        }
//        boolean valid = false;
//        while (!valid) {
//            InputDialog dialog = new InputDialog(IJ.getInstance(), true);
//            dialog.setVisible(true);
//            if (dialog.wasOKed()) {
//                valid = dialog.isValidEntries();
//            } else {
//                return false;
//            }
//        }
    }

    /**
     * Analyses the {@link ImageStack} specified by <code>stack</code>.
     */
    public void analyse() {
        if (monoChrome) {
            colocal = false;
        }
        if (stack != null) {
            IJ.register(this.getClass());
            startTime = System.currentTimeMillis();
            //gaussianRadius = 0.139d / spatialRes; // IsoGaussian filter radius set to 139 nm
            calcParticleRadius(UserVariables.getSpatialRes());
            int i, count;
            int width = stack.getWidth(), height = stack.getHeight();
            if (width > VIS_SIZE || height > VIS_SIZE) {
                scale = 1.0;
            } else {
                scale = VIS_SIZE / Math.max(width, height);
            }

            findParticles(1.0, true, 0, stack.getSize() - 1);

            TextWindow results = new TextWindow(title + " Results", "X\tY\tFrame\tChannel 1 ("
                    + Co_Localise.channels[colocaliser.getChannel1()]
                    + ")\tChannel 2 (" + Co_Localise.channels[colocaliser.getChannel2()]
                    + ")\tChannel 2 " + '\u03C3' + "x\tChannel 2 " + '\u03C3' + "y\t" + '\u03B8',
                    null, 1000, 500);
            results.append(imp.getTitle() + "\n\n");
            TextWindow resultSummary = new TextWindow(title + " Results Summary",
                    "Particle\tType\t% Colocalisation\tDuration (s)\tDisplacement (" + IJ.micronSymbol
                    + "m)\tVelocity (" + IJ.micronSymbol + "m/s)\tDirectionality\tDiffusion Coefficient ("
                    + IJ.micronSymbol + "m^2/s)" + "\tFractal Dimension"
                    + "\tFluorescence Ratio ("
                    + Co_Localise.channels[colocaliser.getChannel2()] + "/"
                    + Co_Localise.channels[colocaliser.getChannel1()]
                    + ")\tAngle Spread\tStep Spread\tDC\tCurvature\tC2 Fluor Area\tC2 Fluor Skew",
                    null, 1200, 500);
            resultSummary.append(imp.getTitle() + "\n\n");

            int n = trajectories.size();
            for (i = 0; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (!(traj.getSize() > minTrajLength && ((traj.getType(colocalThresh) == ParticleTrajectory.COLOCAL)
                        || ((traj.getType(colocalThresh) == ParticleTrajectory.NON_COLOCAL) && !colocal)))) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }

            if (prevRes) {
                previewResults();
            }
            n = trajectories.size();
            mapTrajectories(stack, monoChrome, trajectories, scale, UserVariables.getSpatialRes(), xyPartRad, minTrajLength, timeRes, true);
            for (i = 0, count = 1; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (traj.getSize() > minTrajLength && ((traj.getType(colocalThresh) == ParticleTrajectory.COLOCAL)
                        || ((traj.getType(colocalThresh) == ParticleTrajectory.NON_COLOCAL) && !colocal))) {
                    if (intensPlot) {
                        plotIntensity(i, count);
                    }
                    if (trajPlot) {
                        plotTrajectory(width, height, i, count);
                    }
                    printData(i, resultSummary, count);
                    traj.printTrajectory(count, results, numFormat, title);
                    count++;
                }
            }
            resultSummary.append("\nAnalysis Time (s): " + numFormat.format((System.currentTimeMillis() - startTime) / 1000.0));
            results.append(toString());
            results.setVisible(true);
            resultSummary.setVisible(true);
        }
    }

    /**
     * Median filter and IsoGaussian filter the image specified by
     * <code>processor</code>.
     *
     * @param processor the image to be pre-processed.
     */
    public FloatProcessor preProcess(ByteProcessor processor) {
        //TODO Consider Kalman Filtering for noise reduction
        if (processor == null) {
            return null;
        }
        FloatProcessor fp;
        if (preProcess) {
            TypeConverter tc = new TypeConverter(processor, false);
            fp = (FloatProcessor) tc.convertToFloat(null);
            (new GaussianBlur()).blur(fp, xySigEst);
        } else {
            TypeConverter tc = new TypeConverter(processor, false);
            fp = (FloatProcessor) tc.convertToFloat(null);
        }
        return fp;
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice) {
        if (stack == null) {
            return null;
        }
        int i, noOfImages = stack.getSize(), width = stack.getWidth(), height = stack.getHeight(),
                size = width * height, arraySize = endSlice - startSlice + 1;
        byte c1Pix[] = new byte[size], c2Pix[] = new byte[size],
                c3Pix[] = new byte[size];
        int c1X, c1Y, pSize = 2 * xyPartRad + 1;
        int c2Points[][];
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        double spatialRes = UserVariables.getSpatialRes();
        ParticleArray particles = new ParticleArray(arraySize);
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            IJ.showStatus("Finding Particles...");
            IJ.showProgress(i, noOfImages);
            if (!monoChrome) {
                ColorProcessor colProc = new ColorProcessor(width, height);
                colProc.setRGB(colocaliser.getPixels(colocaliser.getChannel1(), i),
                        colocaliser.getPixels(colocaliser.getChannel2(), i),
                        colocaliser.getPixels(-1, i));
                ((ColorProcessor) colProc).getRGB(c1Pix, c2Pix, c3Pix);
            } else {
                c1Pix = (byte[]) (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte().getPixels());
                c2Pix = null;
            }
            FloatProcessor chan1Proc = preProcess(new ByteProcessor(width, height, c1Pix, null));
            FloatProcessor chan2Proc = !monoChrome
                    ? preProcess(new ByteProcessor(width, height, c2Pix, null)) : null;
            // TODO Maybe express threshold as a percentage? See Ponti et al., 2003
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, FOREGROUND, chan1Proc, chan1MaxThresh, true);
            ByteProcessor C2Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, FOREGROUND, chan2Proc, chan2MaxThresh, true);
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == FOREGROUND) {
                        IsoGaussian c1Gaussian;
                        IsoGaussian c2Gaussian = null;
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                        /*
                         * Remove adjacent Gaussians
                         */
//                        removeAdjacentGaussians(xCoords, yCoords, pixValues, chan1Proc, thisC1Max);
                        IsoGaussianFitter c1GF = new IsoGaussianFitter(xCoords, yCoords, pixValues);
                        c1GF.doFit(xySigEst);
                        if (c1GF.getRSquared() > c1CurveFitTol) {
                            c1Gaussian = new IsoGaussian((c1GF.getX0() + c1X - xyPartRad) * spatialRes,
                                    (c1GF.getY0() + c1Y - xyPartRad) * spatialRes, c1GF.getMag(),
                                    c1GF.getXsig(), c1GF.getYsig(), c1GF.getRSquared());
                        } else {
                            c1Gaussian = new IsoGaussian(c1X
                                    * spatialRes, c1Y * spatialRes,
                                    chan1Proc.getPixelValue(c1X, c1Y), xySigEst,
                                    xySigEst, c1GF.getRSquared());
                        }
                        c2Points = Utils.searchNeighbourhood(c1X, c1Y, (int) Math.round(xyPartRad * searchScale), FOREGROUND,
                                C2Max);
                        if (c2Points != null) {
                            Utils.extractValues(xCoords, yCoords, pixValues,
                                    c2Points[0][0], c2Points[0][1], chan2Proc);
                            IsoGaussianFitter c2GF = new IsoGaussianFitter(xCoords, yCoords, pixValues);
                            c2GF.doFit(xySigEst);
                            c2Gaussian = new IsoGaussian((c2GF.getX0() + c2Points[0][0] - xyPartRad) * spatialRes,
                                    (c2GF.getY0() + c2Points[0][1] - xyPartRad) * spatialRes, c2GF.getMag(),
                                    c2GF.getXsig(), c2GF.getYsig(), c2GF.getRSquared());
//                            FloatProcessor image = new FloatProcessor(25, 25);
//                            c2Gaussian.draw(image, spatialRes);
//                            IJ.saveAs((new ImagePlus("", image)), "TIF", "C:/users/barry05/desktop/tail.tif");
                        }
                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        particles.addDetection(i - startSlice,
                                new Particle((i - startSlice), c1Gaussian, c2Gaussian, null, -1));
                    }
                }
            }
        }
        if (update) {
            updateTrajectories(particles, timeRes, trajMaxStep, chan1MaxThresh, hystDiff, spatialRes, true);
        }
        return particles;
    }

    public void removeAdjacentGaussians(double[] xCoords, double[] yCoords,
            double[][] values, ImageProcessor image, ImageProcessor maxima) {
        int radius = (values.length - 1) / 2;
        int currentX = (int) Math.round(xCoords[radius]);
        int currentY = (int) Math.round(yCoords[radius]);
        int pSize = 3 * radius + 1;
        int localMax[][] = Utils.searchNeighbourhood(currentX, currentY, radius,
                FOREGROUND, maxima);
        int m = localMax.length;
        if (m < 2) {
            return;
        }
        for (int i = 0; i < m; i++) {
            int xc = localMax[i][0];
            int yc = localMax[i][1];
            if (!(xc == currentX && yc == currentY)) {
                double cX[] = new double[pSize];
                double cY[] = new double[pSize];
                double cV[][] = new double[pSize][pSize];
                Utils.extractValues(cX, cY, cV, xc, yc, image);
                MultiGaussFitter gf = new MultiGaussFitter(cV, m);
                double maxz = -Double.MAX_VALUE, minz = -maxz;
                ByteProcessor region = new ByteProcessor(values.length, values.length);
                for (int y = 0; y < pSize; y++) {
                    for (int x = 0; x < pSize; x++) {
                        if (cV[x][y] > maxz) {
                            maxz = cV[x][y];
                        }
                        if (cV[x][y] < minz) {
                            minz = cV[x][y];
                        }
                        if (x < values.length && y < values.length) {
                            region.putPixelValue(x, y, values[x][y]);
                        }
                    }
                }
                double estimates[] = {minz, maxz, xc - cX[0], yc - cY[0],
                    2.0, maxz, currentX - cX[0], currentY - cY[0], 2.0};
                gf.doFit(estimates);
                double params[] = gf.getParams();
                IsoGaussian g = new IsoGaussian((params[2] + cX[0]), (params[3] + cY[0]),
                        params[1], params[4], params[4], gf.getRSquared() - c1CurveFitTol);
                for (int y = 0; y < values.length; y++) {
                    for (int x = 0; x < values.length; x++) {
                        values[x][y] -= g.evaluate(xCoords[x], yCoords[y]);
                        region.putPixelValue(x, y, values[x][y]);
                    }
                }
            }
        }
    }

    public void updateTrajectories(ParticleArray objects, double timeRes, double trajMaxStep,
            double chan1MaxThresh, double hystDiff, double spatialRes, boolean projectVel) {
        if (objects == null) {
            return;
        }
        int i, j, k, m, size;
        int depth = objects.getDepth();
        ParticleTrajectory traj = null;
        Particle last;
        IsoGaussian ch1G;
        double x, y, score, minScore;
        int minScoreIndex;
        ArrayList<Particle> detections;

        for (m = 0; m < depth; m++) {
            IJ.showStatus("Building Trajectories...");
            IJ.showProgress(m, depth);
            for (k = m; (k < depth) && (((k - m)) < trajMaxStep); k++) {
                size = trajectories.size();
                detections = objects.getLevel(k);
                for (j = 0; j < detections.size(); j++) {
                    Particle currentParticle = detections.get(j);
                    if (currentParticle != null) {
                        ch1G = currentParticle.getC1Gaussian();
                        /*
                         * If no trajectories have yet been built, start a new
                         * one:
                         */
                        if (k == m && ch1G.getMagnitude() > chan1MaxThresh * hystDiff && ch1G.getFit() > c1CurveFitTol) {
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
                            for (minScoreIndex = -1, minScore = Double.MAX_VALUE, i = 0; i < size; i++) {
                                traj = (ParticleTrajectory) trajectories.get(i);
                                last = traj.getEnd();
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
                                    // TODO Variation in C1 intensity could be interpreted as movement in Z-direction
                                    if (projectVel) {
                                        traj.projectVelocity(x, y);
                                        double vector1[] = {x, y, currentParticle.getTimePoint(),
                                            ch1G.getMagnitude() / 255.0,
                                            traj.getProjectXVel(),
                                            traj.getProjectYVel()};
                                        double vector2[] = {last.getX(), last.getY(),
                                            last.getTimePoint(), last.getC1Intensity() / 255.0,
                                            traj.getXVelocity(),
                                            traj.getYVelocity()};
                                        score = Utils.calcEuclidDist(vector1, vector2);
                                    } else {
                                        double vector1[] = {x, y, currentParticle.getTimePoint(),
                                            ch1G.getMagnitude() / 255.0};
                                        double vector2[] = {last.getX(), last.getY(),
                                            last.getTimePoint(), last.getC1Intensity() / 255.0};
                                        //TODO Useful cost function described in Vallotton et al., 2003
                                        score = Utils.calcEuclidDist(vector1, vector2);
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
                            // TODO This could probably be extended to find a 'globally optimal' solution, so no 'scoreTol' threshold is required.
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
            for (int l = 0; l < trajectories.size(); l++) {
                traj = (ParticleTrajectory) trajectories.get(l);
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
//            if (IJ.getInstance() != null) {
//                IJ.getTextPanel().append("Frame:\t" + m + "\tTotal Count:\t" + trajectories.size());
//            }
        }
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
     * <code>D = 1 / (1 + &lambda<sub>1</sub> &lambda<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where <code>&lambda<sub>1</sub></code> and
     * <code>&lambda<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and      <code>&lambda<sub>1</sub> <
     * &lambda<sub>2</sub></code>.
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
        double displacement = traj.getDisplacement();
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

    /**
     * Produces a {@link Plot} of normalised intensity in the red and green
     * channels, together with a ratio of red:green intensity, for the
     * trajectory denoted by <code>particleNumber</code>.
     */
    public boolean plotIntensity(int particleNumber, int label) {
        if (trajectories.size() < 1) {
            return false;
        }
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        Particle current = traj.getEnd();
        int size = traj.getSize();
        double xvalues[] = new double[size];
        double redValues[] = new double[size];
        double greenValues[] = new double[size];
        double ratios[] = new double[size];
        double temp, maxVal = -Double.MAX_VALUE, minVal = Double.MAX_VALUE;

        if (traj.getType(colocalThresh) != ParticleTrajectory.COLOCAL) {
            return false;
        }

        for (int i = size - 1; i >= 0; i--) {
            xvalues[i] = current.getTimePoint() * timeRes;
            redValues[i] = current.getC1Intensity() / 255.0d;
            greenValues[i] = current.getC2Intensity() / 255.0d;
            ratios[i] = greenValues[i] / redValues[i];
            if (ratios[i] > maxVal) {
                maxVal = ratios[i];
            }
            temp = Math.min(redValues[i], greenValues[i]);
            if (temp < minVal) {
                minVal = temp;
            }
            current = current.getLink();
        }

        Plot plot = new Plot("Particle " + label + " Intensities",
                "Time (s)", "Normalised Intensity", xvalues, ratios,
                (Plot.X_TICKS + Plot.Y_TICKS + Plot.X_NUMBERS + Plot.Y_NUMBERS));
        plot.changeFont(new Font("Serif", Font.BOLD, 14));
        plot.setLimits(0.0, (stack.getSize() - 1.0) * timeRes, minVal, maxVal);
        plot.setLineWidth(2);
        plot.setColor(Color.BLUE);
        plot.draw();
        plot.setColor(Color.RED);
        plot.addPoints(xvalues, redValues, Plot.LINE);
        plot.setColor(Color.GREEN);
        plot.addPoints(xvalues, greenValues, Plot.LINE);
        plot.show();

        return true;
    }

    /**
     * Constructed trajectories are drawn onto the original image sequence and
     * displayed as a stack sequence.
     */
    public ImageStack mapTrajectories(ImageStack stack, boolean monoChrome, ArrayList<ParticleTrajectory> trajectories, double scale, double spatialRes, double xyPartRad, double minTrajLength, double timeRes, boolean tracks) {
        if (stack == null) {
            return null;
        }
        int i, j, width = (int) Math.round(stack.getWidth() * scale), height = (int) Math.round(stack.getHeight() * scale),
                type, frames = stack.getSize(), radius = (int) Math.round(xyPartRad * scale);
        double lastX, lastY;
        ImageStack outputStack = new ImageStack(width, height);
        Particle current;
        ParticleTrajectory traj;
        int length, n = trajectories.size(), count;
        ImageProcessor processor;
        Rectangle bounds;
        ByteProcessor trajImage;
        double scaledSR = scale / spatialRes;
        int lastTP;
        double ptScale = ParticleTrajectory.scale;

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
        for (i = 0, count = 1; i < n; i++) {
            IJ.showStatus("Mapping Output...");
            IJ.showProgress(i, n);
            traj = (ParticleTrajectory) (trajectories.get(i));
            Color thiscolor = new Color(r.nextInt(256), r.nextInt(256), r.nextInt(256));
            length = traj.getSize();
            type = traj.getType(colocalThresh);
            bounds = traj.getBounds();
            trajImage = new ByteProcessor(bounds.width, bounds.height);
            trajImage.setColor(Color.white);
            trajImage.fill();
            trajImage.setColor(Color.black);
            if (length > minTrajLength && ((type == ParticleTrajectory.COLOCAL)
                    || ((type == ParticleTrajectory.NON_COLOCAL) && !colocal))) {
                current = traj.getEnd();
                lastX = current.getX();
                lastY = current.getY();
                lastTP = current.getTimePoint();
                current = current.getLink();
                while (current != null) {
                    for (j = frames - 1; j >= lastTP; j--) {
                        processor = outputStack.getProcessor(j + 1);
                        if (!monoChrome) {
//                            if (type == ParticleTrajectory.NON_COLOCAL) {
                            processor.setColor(thiscolor);
//                            }
                        } else {
                            processor.setValue(255);
                        }
                        if (j - 1 < lastTP) {
                            markParticle(processor, (int) Math.round(lastX * scaledSR) - radius,
                                    (int) Math.round(lastY * scaledSR) - radius, radius, true, "" + count);
                        }
                        if (tracks && j <= lastTP + 20.0 / timeRes) {
                            int x = (int) (Math.round(current.getX() * scaledSR));
                            int y = (int) (Math.round(current.getY() * scaledSR));
                            processor.drawLine(x, y, (int) Math.round(lastX * scaledSR),
                                    (int) Math.round(lastY * scaledSR));
                        }
                    }
                    if (tracks) {
                        trajImage.drawLine((int) Math.round(ptScale * current.getX() - traj.getBounds().x),
                                (int) Math.round(ptScale * current.getY() - traj.getBounds().y),
                                (int) Math.round(ptScale * lastX - traj.getBounds().x),
                                (int) Math.round(ptScale * lastY - traj.getBounds().y));
                    }
                    lastX = current.getX();
                    lastY = current.getY();
                    lastTP = current.getTimePoint();
                    current = current.getLink();
                }
                processor = outputStack.getProcessor(lastTP + 1);
                if (!monoChrome) {
                    processor.setColor(thiscolor);
                } else {
                    processor.setValue(255);
                }
                markParticle(processor, (int) Math.round(lastX * scaledSR) - radius,
                        (int) Math.round(lastY * scaledSR) - radius, radius, true, "" + count);
                count++;
            }
        }
        (new ImagePlus("Trajectories", outputStack)).show();

        return outputStack;
    }

    void markParticle(ImageProcessor processor, int x, int y, int radius, boolean string, String label) {
        processor.drawOval(x, y, 2 * radius, 2 * radius);
        if (string) {
            processor.drawString(label, x, y);
        }
    }

    /*
     * public static void setGaussianRadius(double gaussianRadius) {
     * Timelapse_Analysis.gaussianRadius = gaussianRadius; }
     */
    public void setSigmaEstimate(double sigmaEstimate) {
        this.xySigEst = sigmaEstimate;
    }

    public int getParticleRadius() {
        return xyPartRad;
    }

    public void calcParticleRadius(double spatialRes) {
        double airyRad = 1.22 * LAMBDA / (2.0 * NUM_AP); //Airy radius
        xyPartRad = (int) Math.ceil(airyRad / (spatialRes * 1000.0));
        xySigEst = airyRad / (2.0 * spatialRes * 1000.0);
    }

    public boolean removeTrajectory(int index) {
        if (trajectories.size() < 1) {
            return false;
        }
        trajectories.remove(index);

        return true;
    }

    /*
     * public static void setC1SigmaTol(double c1SigmaTol) {
     * Timelapse_Analysis.c1SigmaTol = c1SigmaTol; }
     *
     * public static void setC2SigmaTol(double c2SigmaTol) {
     * Timelapse_Analysis.c2SigmaTol = c2SigmaTol; }
     */
    public void previewResults() {
        if (trajectories.size() < 1) {
            return;
        }
        int remove = 0;
        boolean done = false;
        stack = imp.getImageStack();
        while (!done) {
            GenericDialog dialog = new GenericDialog("Verify Individual Trajectories", IJ.getInstance());
            if (remove > trajectories.size() - 1) {
                remove--;
            }
            dialog.addSlider("Select Trajectory:", 0, trajectories.size() - 1, remove);
            dialog.setModal(false);
            dialog.setOKLabel("Reject");
            ResultsPreviewer previewer = new ResultsPreviewer();
            dialog.addDialogListener(previewer);
            dialog.showDialog();
            while (!(dialog.wasCanceled() || dialog.wasOKed()));
            previewer.preview.close();
            if (dialog.wasOKed()) {
                remove = (int) dialog.getNextNumber();
                removeTrajectory(remove);
                if (trajectories.size() < 1) {
                    return;
                }
            } else {
                done = true;
            }
        }
    }

    public Color getDrawColor(int key) {
        switch (key) {
            case Co_Localise.RED:
                return Color.red;
            case Co_Localise.GREEN:
                return Color.green;
            case Co_Localise.BLUE:
                return Color.blue;
            default:
                return Color.white;
        }
    }

    public void setColocaliser(Co_Localise colocaliser) {
        this.colocaliser = colocaliser;
    }

    public ArrayList<ParticleTrajectory> getTrajectories() {
        return trajectories;
    }

    public void setTrajectories(ArrayList<ParticleTrajectory> trajectories) {
        this.trajectories = trajectories;
    }

    public ImageStack getStack() {
        return stack;
    }

    public boolean isMonoChrome() {
        return monoChrome;
    }

    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Timelapse_Analysis other = (Timelapse_Analysis) obj;
        if (Double.doubleToLongBits(this.scale) != Double.doubleToLongBits(other.scale)) {
            return false;
        }
        if (this.xyPartRad != other.xyPartRad) {
            return false;
        }
        if (this.trajectories != other.trajectories && (this.trajectories == null || !this.trajectories.equals(other.trajectories))) {
            return false;
        }
        if (this.imp != other.imp && (this.imp == null || !this.imp.equals(other.imp))) {
            return false;
        }
        if (this.stack != other.stack && (this.stack == null || !this.stack.equals(other.stack))) {
            return false;
        }
        if (this.startTime != other.startTime) {
            return false;
        }
        if (this.numFormat != other.numFormat && (this.numFormat == null || !this.numFormat.equals(other.numFormat))) {
            return false;
        }
        return true;
    }
    
    private class ResultsPreviewer implements DialogListener {

        public ImagePlus preview;

        public ResultsPreviewer() {
        }

        public boolean dialogItemChanged(GenericDialog dialog, java.awt.AWTEvent event) {
            if (preview != null) {
                preview.close();
            }
            if (dialog.wasOKed() || dialog.wasCanceled()) {
                return true;
            }
            preview = mapTrajectories((int) ((Scrollbar) (dialog.getSliders().get(0))).getValue(), UserVariables.getSpatialRes());
            preview.show();
            return true;
        }

        public ImagePlus mapTrajectories(int trajNumber, double spatialRes) {
            ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(trajNumber));
            stack = imp.getStack();
            Rectangle bounds = traj.getBounds();
            int visScale = 2, border = 2 * xyPartRad;
            int i, j, width = (bounds.width + 2 * border) * visScale, height = (bounds.height + 2 * border) * visScale,
                    currentX, currentY, lastX, lastY, type, frames = stack.getSize(),
                    xOffset = bounds.x - border, yOffset = bounds.y - border;
            ImageStack outputStack = new ImageStack(width, height);
            Particle current;
            ImageProcessor processor;
            double lastTP;
            Rectangle newBounds = new Rectangle(bounds.x - border, bounds.y - border,
                    bounds.width + 2 * border, bounds.height + 2 * border);

            for (i = 0; i < frames; i++) {
                if (monoChrome) {
                    processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte());
                } else {
                    processor = (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToRGB());
                }
                processor.setRoi(newBounds);
                outputStack.addSlice("" + i, (processor.crop()).resize(width, height));
            }

            type = traj.getType(0.1);
            current = traj.getEnd();
            currentX = lastX = (int) (Math.round(current.getX() / spatialRes - xOffset)) * visScale;
            currentY = lastY = (int) (Math.round(current.getY() / spatialRes - yOffset)) * visScale;
            lastTP = current.getTimePoint();
            current = current.getLink();
            while (current != null) {
                for (j = frames - 1; j >= lastTP; j--) {
                    if (j <= lastTP + 10.0 / timeRes) {
                        processor = outputStack.getProcessor(j + 1);
                        if (!monoChrome) {
//                            if (type == ParticleTrajectory.NON_COLOCAL) {
                            processor.setColor(Color.WHITE);
//                            } else {
//                                processor.setColor(getDrawColor(colocaliser.getChannel2()));
//                            }
                        } else {
                            processor.setColor(255);
                        }
                        currentX = (int) (Math.round(current.getX() / spatialRes - xOffset)) * visScale;
                        currentY = (int) (Math.round(current.getY() / spatialRes - yOffset)) * visScale;
                        processor.drawLine(currentX, currentY, lastX, lastY);
                    }
                }
                processor = outputStack.getProcessor(j + 1);
                if (!monoChrome) {
                    processor.setColor(Color.WHITE);
                } else {
                    processor.setColor(255);
                }
                processor.drawOval(currentX - 2, currentY - 2, 5, 5);
                lastX = currentX;
                lastY = currentY;
                lastTP = current.getTimePoint();
                current = current.getLink();
            }
            return new ImagePlus("Trajectory Number " + (trajNumber + 1), outputStack);
        }
    }
}
