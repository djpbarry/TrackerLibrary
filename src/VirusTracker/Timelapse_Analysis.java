package VirusTracker;

import AnaMorf.Utilities;
import EMSeg.Utils;
import EMSeg.FractalEstimator;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.Toolkit;
import java.io.File;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

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

    protected static double spatialRes = 84.0 / 1000.0, //Spatial resolution in nm/pixel
            timeRes = 1.0d, //Time resolution in s/frame;
            virusDiameter = 350.0,
            chan1MaxThresh = 100.0, //Threshold value for local maxima in channel 1
            chan2MaxThresh = 0.0, //Threshold value for local maxima in channel 2
            //gaussianRadius, //Gaussian filter radius
            trajMaxStep = 1.1, //Tolerance used in evaluating likely trajectories
            //TODO Consider performing a series of analyses with incresing values for trajMaxStep. This may allow discrimination between different modes of motion by identifying slowly-moving particles first (and subsequently removing them from the data set) and then proceeding to more complex trajectories.
            minTrajLength = 0.0 / timeRes, //Minimum trajectory length output in results
            hystDiff = 1.25;
    protected double scale, //Scale factor for visualisation
            xySigEst; //Initial estimate of standard deviation for IsoGaussian fitting
    protected int xyPartRad; //Radius over which to draw particles in visualisation
    public final static int FOREGROUND = 255, //Integer value of foreground pixels
            SHOW_RESULTS = -1;
    protected final static double VIS_SIZE = 750.0,
            LAMBDA = 650.0, //Wavelength of light
            NUM_AP = 1.4; //Numerical aperture of system
    protected final static double curveFitTol = 0.8d; //Tolerance used in determining fit of IsoGaussian curves
    //protected static double c1SigmaTol = 3.0, c2SigmaTol = 3.0;
    protected ArrayList<ParticleTrajectory> trajectories = new ArrayList<ParticleTrajectory>(); //Trajectories of the detected particles
    protected ImagePlus imp; //The active image stack
    protected ImageStack stack;
    private long startTime;
    protected DecimalFormat numFormat = new DecimalFormat("0.000");
    public static String title = "Virus Tracker v3.09";
    protected static boolean colocal = false, msdPlot = false, intensPlot = false,
            preProcess = true, trajPlot = false, prevRes = false;
    protected Co_Localise colocaliser;
    protected boolean monoChrome;
    private double noiseTol = 0.2;

    /*
     * public static void main(String args[]) { File image =
     * Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\Tracking Test
     * Sequences")); ImageStack stack = Utils.buildStack(image); ImagePlus imp =
     * new ImagePlus("Stack", stack); Timelapse_Analysis instance = new
     * Timelapse_Analysis(imp); if (instance.showDialog()) { instance.analyse();
     * } return;
    }
     */
    public Timelapse_Analysis(double spatialRes, double timeRes, double trajMaxStep,
            double chan1MaxThresh, double hystDiff, boolean monoChrome, ImagePlus imp, double scale, double minTrajLength) {
        Timelapse_Analysis.spatialRes = spatialRes;
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
        boolean valid = false;
        while (!valid) {
            valid = true;
            stack = imp.getImageStack();
            if (stack.getProcessor(1).getNChannels() > 1) {
                monoChrome = false;
            } else {
                monoChrome = true;
            }
            InputDialog dialog = new InputDialog(IJ.getInstance(), true);
            dialog.setVisible(true);
            if (dialog.wasOKed()) {
                if (dialog.isValidEntries()) {
                    valid = true;
                } else {
                    valid = false;
                }
            } else {
                return false;
            }
        }
        return true;
    }

    /**
     * Analyses the {@link ImageStack} specified by
     * <code>stack</code>.
     */
    public void analyse() {
        if (stack != null) {
            IJ.register(this.getClass());
            startTime = System.currentTimeMillis();
            //gaussianRadius = 0.139d / spatialRes; // IsoGaussian filter radius set to 139 nm
            calcParticleRadius(); //Virus particles are approximately 350nm in diameter
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
                if (!(traj.getSize() > minTrajLength && ((traj.getType() == ParticleTrajectory.COLOCAL)
                        || ((traj.getType() == ParticleTrajectory.NON_COLOCAL) && !colocal)))) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }

            if (prevRes) {
                previewResults();
            }
            n = trajectories.size();
            mapTrajectories();
            for (i = 0, count = 1; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (traj.getSize() > minTrajLength && ((traj.getType() == ParticleTrajectory.COLOCAL)
                        || ((traj.getType() == ParticleTrajectory.NON_COLOCAL) && !colocal))) {
                    if (intensPlot) {
                        plotIntensity(i, count);
                    }
                    if (trajPlot) {
                        plotTrajectory(width, height, i, count);
                    }
                    printData(i, resultSummary, count);
                    traj.printTrajectory(count, results, numFormat);
                    count++;
                }
            }
            resultSummary.append("\nAnalysis Time (s): " + numFormat.format((System.currentTimeMillis() - startTime) / 1000.0));
            results.append(toString());
            results.setVisible(true);
            resultSummary.setVisible(true);
        }
        return;
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
        ParticleArray particles = new ParticleArray(arraySize);
        ByteProcessor nextC1Max = new ByteProcessor(width, height);
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            IJ.showStatus("Analysing Frame " + i);
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
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, FOREGROUND, chan1Proc, chan1MaxThresh, true); // TODO Maybe express threshold as a percentage? See Ponti et al., 2003
            ByteProcessor C2Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, FOREGROUND, chan2Proc, chan2MaxThresh, true);
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == FOREGROUND) {
                        IsoGaussian c1Gaussian = null;
                        IsoGaussian c2Gaussian = null;
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                        /*
                         * Remove adjacent Gaussians
                         */
                        removeAdjacentGaussians(xCoords, yCoords, pixValues, chan1Proc, thisC1Max);
                        IsoGaussianFitter c1GF = new IsoGaussianFitter(xCoords, yCoords, pixValues);
                        c1GF.doFit();
                        // TODO Estiamtes of intensity need to consider particles moving into/out of focal plane
                        //if (c1GF.getXsig() < (c1SigmaTol * xySigEst)) {
                        /*
                         * if (c1GF.getRSquared() > curveFitTol) {
                         */
                        c1Gaussian = new IsoGaussian((c1GF.getX0() + c1X - xyPartRad) * spatialRes,
                                (c1GF.getY0() + c1Y - xyPartRad) * spatialRes, c1GF.getMag(),
                                c1GF.getXsig(), c1GF.getYsig(), c1GF.getRSquared() - curveFitTol);
                        /*
                         * } else { c1Gaussian = new IsoGaussian(c1X *
                         * spatialRes, c1Y * spatialRes,
                         * chan1Proc.getPixelValue(c1X, c1Y), xySigEst,
                         * xySigEst, c1GF.getRSquared() - curveFitTol);
                        }
                         */
                        /*
                         * c2Points = Utils.searchNeighbourhood(c1X, c1Y, (int)
                         * Math.round(xyPartRad * searchScale), FOREGROUND,
                         * C2Max); if (c2Points != null) {
                         * Utils.extractValues(xCoords, yCoords, pixValues,
                         * c2Points[0][0], c2Points[0][1], chan2Proc);
                         * IsoGaussianFitter c2GF = new
                         * IsoGaussianFitter(xCoords, yCoords, pixValues);
                         * c2GF.doFit(xySigEst, 0.0); //if (c2GF.getXsig() <
                         * (c2SigmaTol * xySigEst)) { c2Gaussian = new
                         * IsoGaussian((c2GF.getX0() + c2Points[0][0] -
                         * xyPartRad) * spatialRes, (c2GF.getY0() +
                         * c2Points[0][1] - xyPartRad) * spatialRes,
                         * c2GF.getMag(), c2GF.getXsig(), c2GF.getYsig(),
                         * c2GF.getRSquared() - curveFitTol); //}
                        }
                         */
                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        if (c1Gaussian != null) {
                            particles.addDetection(i - startSlice, c1Gaussian, c2Gaussian);
                        }
                        //}
                    }
                }
            }
        }
        if (update) {
            updateTrajectories(particles);
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
                /*
                 * TODO try a sum of Gaussians fitter & check work of K. Lidke
                 * or Expectation Maximisation
                 */
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
                        params[1], params[4], params[4], gf.getRSquared() - curveFitTol);
                for (int y = 0; y < values.length; y++) {
                    for (int x = 0; x < values.length; x++) {
                        values[x][y] -= g.evaluate(xCoords[x], yCoords[y]);
                        region.putPixelValue(x, y, values[x][y]);
                    }
                }
            }
        }
        return;
    }

    public void updateTrajectories(ParticleArray objects) {
        //TODO Trajectories should be re-analysed once all frames have been scanned to optimise trajectories
        if (objects == null) {
            return;
        }
        int i, j, k, m, size;
        ParticleTrajectory traj = null;
        Particle last;
        IsoGaussian ch1G, ch2G;
        double x, y, score, greenMag, minScore;
        int minScoreIndex;
        ArrayList detections;

        for (m = 0; m < objects.getDepth(); m++) {
            for (k = m; (k < objects.getDepth()) && (((k - m) * timeRes) < trajMaxStep); k++) {
                size = trajectories.size();
                detections = objects.getLevel(k);
                for (j = 0; j < detections.size(); j++) {
                    ch1G = ((IsoGaussian[]) detections.get(j))[0];
                    if (ch1G != null) {
                        ch2G = ((IsoGaussian[]) detections.get(j))[1];
                        /*
                         * If no trajectories have yet been built, start a new
                         * one:
                         */
                        if (k == m && ch1G.getMagnitude() > chan1MaxThresh * hystDiff && ch1G.getFit() > 0.0) {
                            traj = new ParticleTrajectory(timeRes, spatialRes);
                            traj.addPoint(k * timeRes, ch1G, ch2G);
                            trajectories.add(traj);
                            /*
                             * Otherwise, determine whether the current particle
                             * belongs to a pre-existing trajectory:
                             */
                        } else {
                            for (minScoreIndex = -1, minScore = Double.MAX_VALUE, i = 0; i < size; i++) {
                                traj = (ParticleTrajectory) trajectories.get(i);
                                last = traj.getEnd();
                                if ((last != null) && (last.getTimePoint() == m * timeRes) && k != m) {
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
                                    if (ch2G != null && !(ch2G instanceof RandomDistribution)) {
                                        x = (ch1G.getX() + ch2G.getX()) / 2.0d;
                                        y = (ch1G.getY() + ch2G.getY()) / 2.0d;
                                        greenMag = ch2G.getMagnitude();
                                    } else {
                                        x = ch1G.getX();
                                        y = ch1G.getY();
                                        greenMag = 0.0;
                                    }
                                    traj.projectVelocity(x, y);
                                    double vector1[] = {x, y, k * timeRes,
                                        ch1G.getMagnitude() / 255.0, greenMag / 255.0, traj.getProjectXVel(),
                                        traj.getProjectYVel()};
                                    double vector2[] = {last.getX(), last.getY(),
                                        last.getTimePoint(), last.getC1Intensity() / 255.0,
                                        last.getC2Intensity() / 255.0, traj.getXVelocity(),
                                        traj.getYVelocity()};
                                    //TODO Useful cost function described in Vallotton et al., 2003
                                    score = Utils.calcEuclidDist(vector1, vector2);
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
                            if (minScoreIndex > -1) {
                                traj = (ParticleTrajectory) trajectories.get(minScoreIndex);
                            }
                            if ((minScore < trajMaxStep) && (minScore < traj.getTempScore())) {
                                traj.addTempPoint(k, ch1G, ch2G, minScore, j, k);
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
                        traj.checkDetections(timeRes * temp.getTimePoint(), temp.getC1Gaussian(),
                                temp.getC2Gaussian());
                        objects.nullifyDetection(col, row);
                    }
                }
            }
        }
        return;
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
        width *= spatialRes;
        height *= spatialRes;

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
     * <code>particleNumber</code>. Directionality (
     * <code>D</code>) is calculated according to: <br> <br>
     * <code>D = 1 / (1 + &lambda<sub>1</sub> &lambda<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where
     * <code>&lambda<sub>1</sub></code> and
     * <code>&lambda<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and
     * <code>&lambda<sub>1</sub> <
     * &lambda<sub>2</sub></code>.
     *
     */
    public boolean printData(int particleNumber, TextWindow output, int label) {
        if (trajectories.size() < 1) {
            return false;
        }
        //TODO Add MSD to results
        DecimalFormat decFormat = new DecimalFormat("0.000");
        DecimalFormat msdFormat = new DecimalFormat("0.000000");
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        if (traj == null) {
            return false;
        }
        traj.smooth();
        traj.calcMSD(label, -1, msdPlot);
        traj.calcAngleSpread();
        traj.calcStepSpread();
        traj.calcDirectionality();
        double displacement = traj.getDisplacement();
        double duration = traj.getSize() * timeRes;
        int type = traj.getType();
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
     * trajectory denoted by
     * <code>particleNumber</code>.
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

        if (traj == null || traj.getType() != ParticleTrajectory.COLOCAL) {
            return false;
        }

        for (int i = size - 1; i >= 0; i--) {
            xvalues[i] = current.getTimePoint();
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
    public boolean mapTrajectories() {
        if (imp == null) {
            return false;
        }
        stack = imp.getStack();
        int i, j, width = (int) Math.round(stack.getWidth() * scale), height = (int) Math.round(stack.getHeight() * scale),
                textX, textY, type, frames = stack.getSize(), radius = (int) Math.round(xyPartRad * scale);
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
        FractalEstimator boxCount;

        if (n < 1) {
            return false;
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

        for (i = 0, count = 1; i < n; i++) {
            traj = (ParticleTrajectory) (trajectories.get(i));
            length = traj.getSize();
            type = traj.getType();
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
                textX = (int) (Math.round(current.getX() * scaledSR));
                textY = (int) (Math.round(current.getY() * scaledSR));
                lastTP = (int) Math.round(current.getTimePoint() / timeRes);
                current = current.getLink();
                while (current != null) {
                    for (j = frames - 1; j >= lastTP; j--) {
                        processor = outputStack.getProcessor(j + 1);
                        if (!monoChrome) {
                            if (type == ParticleTrajectory.NON_COLOCAL) {
                                processor.setColor(getDrawColor(colocaliser.getChannel1()));
                            } else {
                                processor.setColor(getDrawColor(colocaliser.getChannel2()));
                            }
                        } else {
                            processor.setValue(255);
                        }
                        int x = (int) (Math.round(current.getX() * scaledSR));
                        int y = (int) (Math.round(current.getY() * scaledSR));
                        processor.drawLine(x, y, (int) Math.round(lastX * scaledSR),
                                (int) Math.round(lastY * scaledSR));
                        if (j - 1 < lastTP) {
                            processor.drawOval(((int) Math.round(lastX * scaledSR) - radius), ((int) Math.round(lastY * scaledSR) - radius), 2 * radius, 2 * radius);
                            System.out.println("j: " + j + " lastTP: " + lastTP);
                        }
                    }
                    trajImage.drawLine((int) Math.round(ptScale * current.getX() - traj.getBounds().x),
                            (int) Math.round(ptScale * current.getY() - traj.getBounds().y),
                            (int) Math.round(ptScale * lastX - traj.getBounds().x),
                            (int) Math.round(ptScale * lastY - traj.getBounds().y));
                    lastX = current.getX();
                    lastY = current.getY();
                    lastTP = (int) Math.round(current.getTimePoint() / timeRes);
                    current = current.getLink();
                }
                /*
                 * for (j = 0; j < frames; j++) { processor =
                 * outputStack.getProcessor(j + 1); if (!monoChrome) { if (type
                 * == ParticleTrajectory.NON_COLOCAL) {
                 * processor.setColor(getDrawColor(colocaliser.getChannel1()));
                 * } else {
                 * processor.setColor(getDrawColor(colocaliser.getChannel2()));
                 * } } else { processor.setValue(255); } processor.drawString(""
                 * + count, textX, textY);
                }
                 */
                count++;
            }
            boxCount = new FractalEstimator();
            double[] fractalDims = boxCount.do2DEstimate(trajImage);
            if (fractalDims != null) {
                traj.setFracDim(fractalDims[0]);
            } else {
                traj.setFracDim(Double.NaN);
            }
        }
        (new ImagePlus("Trajectories", outputStack)).show();

        return true;
    }

    public static void setChan1MaxThresh(double chan1MaxThresh) {
        Timelapse_Analysis.chan1MaxThresh = chan1MaxThresh;
    }

    public static void setChan2MaxThresh(double chan2MaxThresh) {
        Timelapse_Analysis.chan2MaxThresh = chan2MaxThresh;
    }

    public static void setSpatialRes(double spatialRes) {
        Timelapse_Analysis.spatialRes = spatialRes;
    }

    /*
     * public static void setGaussianRadius(double gaussianRadius) {
     * Timelapse_Analysis.gaussianRadius = gaussianRadius;
    }
     */
    public void setSigmaEstimate(double sigmaEstimate) {
        this.xySigEst = sigmaEstimate;
    }

    public int getParticleRadius() {
        return xyPartRad;
    }

    public void calcParticleRadius() {
        double airyRad = 1.22 * LAMBDA / (2.0 * NUM_AP); //Airy radius
        xyPartRad = (int) Math.ceil(airyRad / (spatialRes * 1000.0));
        xySigEst = airyRad / (2.0 * spatialRes * 1000.0);
    }

    public static void setPreProcess(boolean preProcess) {
        Timelapse_Analysis.preProcess = preProcess;
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
     * Timelapse_Analysis.c2SigmaTol = c2SigmaTol;
    }
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
            while (dialog.isVisible());
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
        return;
    }

    public static double getSpatialRes() {
        return spatialRes;
    }

    public static double getChan1MaxThresh() {
        return chan1MaxThresh;
    }

    public static double getChan2MaxThresh() {
        return chan2MaxThresh;
    }

    Color getDrawColor(int key) {
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

    public String toString() {
        Date current = new Date();
        DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
        if (colocaliser != null) {
            return title + "\t" + dateFormat.format(current)
                    + "\n\nChannel 1:\t" + Co_Localise.channels[colocaliser.getChannel1()]
                    + "\nChannel 2:\t" + Co_Localise.channels[colocaliser.getChannel2()]
                    + "\nSpatial Resolution (nm/pixel):\t" + numFormat.format(spatialRes * 1000.0)
                    + "\nFrames per Second (/s):\t" + numFormat.format(1.0 / timeRes)
                    + "\nVirus Diameter (nm):\t" + numFormat.format(virusDiameter)
                    + "\nMinimum Peak Size (" + Co_Localise.channels[colocaliser.getChannel1()]
                    + "):\t" + numFormat.format(chan1MaxThresh)
                    + "\nMinimum Peak Size (" + Co_Localise.channels[colocaliser.getChannel2()]
                    + "):\t" + numFormat.format(chan2MaxThresh)
                    + "\nTrajectory Step Tolerance:\t" + numFormat.format(trajMaxStep)
                    + "\nMinimum Trajectory Length (s):\t" + numFormat.format(minTrajLength * timeRes)
                    + "\nPre-Process Stack:\t" + preProcess;
        } else {
            return title + "\t" + dateFormat.format(current)
                    + "\nSpatial Resolution (nm/pixel):\t" + numFormat.format(spatialRes * 1000.0)
                    + "\nFrames per Second (/s):\t" + numFormat.format(1.0 / timeRes)
                    + "\nVirus Diameter (nm):\t" + numFormat.format(virusDiameter)
                    + "\nMinimum Peak Size (Chan 1):\t" + numFormat.format(chan1MaxThresh)
                    + "\nMinimum Peak Size (Chan 2):\t" + numFormat.format(chan2MaxThresh)
                    + "\nTrajectory Step Tolerance:\t" + numFormat.format(trajMaxStep)
                    + "\nMinimum Trajectory Length (s):\t" + numFormat.format(minTrajLength * timeRes)
                    + "\nPre-Process Stack:\t" + preProcess;
        }
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
        if (this.colocaliser != other.colocaliser && (this.colocaliser == null || !this.colocaliser.equals(other.colocaliser))) {
            return false;
        }
        return true;
    }

    public class InputDialog extends javax.swing.JDialog {

        private ImagePlus preview;
        private int width, height, slice = 0;
        private boolean dialogOKed = false, validEntries = true;

        public InputDialog() {
            super();
        }

        /**
         * Creates new form InputDialog
         */
        public InputDialog(java.awt.Frame parent, boolean modal) {
            super(parent, modal);
            initComponents();
        }

        /**
         * This method is called from within the constructor to initialize the
         * form. WARNING: Do NOT modify this code. The content of this method is
         * always regenerated by the Form Editor.
         */
        private void initComponents() {

            jLabel1 = new javax.swing.JLabel();
            jLabel2 = new javax.swing.JLabel();
            jLabel3 = new javax.swing.JLabel();
            channel1Combo = new javax.swing.JComboBox();
            channel2Combo = new javax.swing.JComboBox();
            jLabel4 = new javax.swing.JLabel();
            virusDiamField = new javax.swing.JTextField();
            jLabel5 = new javax.swing.JLabel();
            jLabel6 = new javax.swing.JLabel();
            spatialResField = new javax.swing.JTextField();
            jLabel7 = new javax.swing.JLabel();
            jLabel8 = new javax.swing.JLabel();
            timeResField = new javax.swing.JTextField();
            jLabel9 = new javax.swing.JLabel();
            jLabel10 = new javax.swing.JLabel();
            minTrajLengthField = new javax.swing.JTextField();
            jLabel11 = new javax.swing.JLabel();
            jLabel12 = new javax.swing.JLabel();
            trajStepTolField = new javax.swing.JTextField();
            jLabel13 = new javax.swing.JLabel();
            jLabel14 = new javax.swing.JLabel();
            ch1MinPeakField = new javax.swing.JTextField();
            ch2MinPeakField = new javax.swing.JTextField();
            prevDetectBox = new javax.swing.JCheckBox();
            prevScroll = new java.awt.Scrollbar();
            jLabel15 = new javax.swing.JLabel();
            prevField = new javax.swing.JTextField();
            preProcessBox = new javax.swing.JCheckBox();
            colocalBox = new javax.swing.JCheckBox();
            intensPlotBox = new javax.swing.JCheckBox();
            trajPlotBox = new javax.swing.JCheckBox();
            prevResBox = new javax.swing.JCheckBox();
            analyseButton = new javax.swing.JButton();

            setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
            setTitle(Timelapse_Analysis.title);

            jLabel1.setText("Channel 2 will be colocalised with Channel 1");

            jLabel2.setText("Channel 1:");

            jLabel3.setText("Channel 2:");

            channel1Combo.setModel(new javax.swing.DefaultComboBoxModel(Co_Localise.channels));
            channel1Combo.setSelectedItem(Co_Localise.channels[Co_Localise.RED]);
            channel1Combo.setMaximumSize(new java.awt.Dimension(72, 20));
            channel1Combo.setMinimumSize(new java.awt.Dimension(72, 20));
            channel1Combo.setPreferredSize(new java.awt.Dimension(72, 20));
            channel1Combo.setEnabled(!monoChrome);

            channel2Combo.setModel(new javax.swing.DefaultComboBoxModel(Co_Localise.channels));
            channel2Combo.setSelectedItem(Co_Localise.channels[Co_Localise.GREEN]);
            channel2Combo.setMaximumSize(new java.awt.Dimension(72, 20));
            channel2Combo.setMinimumSize(new java.awt.Dimension(72, 20));
            channel2Combo.setPreferredSize(new java.awt.Dimension(72, 20));
            channel2Combo.setEnabled(!monoChrome);

            jLabel4.setText("Virus Diameter:");

            virusDiamField.setText(virusDiameter + "");
            virusDiamField.setMaximumSize(new java.awt.Dimension(72, 20));
            virusDiamField.setMinimumSize(new java.awt.Dimension(72, 20));

            jLabel5.setText("nm");

            jLabel6.setText("Spatial Resolution:");

            spatialResField.setText(spatialRes * 1000.0 + "");
            spatialResField.setMaximumSize(new java.awt.Dimension(72, 20));
            spatialResField.setMinimumSize(new java.awt.Dimension(72, 20));

            jLabel7.setText("<html>nm pixel<sup>-1</sup></html>");

            jLabel8.setText("Frames per Second:");

            timeResField.setText(numFormat.format(1.0 / timeRes));
            timeResField.setMaximumSize(new java.awt.Dimension(72, 20));
            timeResField.setMinimumSize(new java.awt.Dimension(72, 20));

            jLabel9.setText("Hz");

            jLabel10.setText("Minimum Trajectory Length:");

            minTrajLengthField.setText(numFormat.format(minTrajLength * timeRes));
            minTrajLengthField.setMaximumSize(new java.awt.Dimension(72, 20));
            minTrajLengthField.setMinimumSize(new java.awt.Dimension(72, 20));

            jLabel11.setText("s");

            jLabel12.setText("Trajectory Step Tolerance:");

            trajStepTolField.setText(trajMaxStep + "");
            trajStepTolField.setMaximumSize(new java.awt.Dimension(72, 20));
            trajStepTolField.setMinimumSize(new java.awt.Dimension(72, 20));

            jLabel13.setText("Minimum Peak Size (Ch 1):");

            jLabel14.setText("Minimum Peak Size (Ch 2):");

            ch1MinPeakField.setText(chan1MaxThresh + "");
            ch1MinPeakField.setMaximumSize(new java.awt.Dimension(72, 20));
            ch1MinPeakField.setMinimumSize(new java.awt.Dimension(72, 20));

            ch2MinPeakField.setText(chan2MaxThresh + "");
            ch2MinPeakField.setMaximumSize(new java.awt.Dimension(72, 20));
            ch2MinPeakField.setMinimumSize(new java.awt.Dimension(72, 20));
            ch2MinPeakField.setEnabled(!monoChrome);

            prevDetectBox.setText("Preview Detections");
            prevDetectBox.addActionListener(new java.awt.event.ActionListener() {

                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    prevDetectBoxActionPerformed(evt);
                }
            });

            prevScroll.setEnabled(prevDetectBox.isSelected());
            prevScroll.setMaximum(stack.getSize());
            prevScroll.setOrientation(java.awt.Scrollbar.HORIZONTAL);
            prevScroll.addAdjustmentListener(new java.awt.event.AdjustmentListener() {

                public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
                    prevScrollAdjustmentValueChanged(evt);
                }
            });

            jLabel15.setText("Preview Frame:");

            prevField.setText(prevScroll.getValue() + "");
            prevField.setEnabled(prevDetectBox.isSelected());
            prevField.setMaximumSize(new java.awt.Dimension(72, 20));
            prevField.setMinimumSize(new java.awt.Dimension(72, 20));
            prevField.setPreferredSize(new java.awt.Dimension(72, 20));
            prevField.addActionListener(new java.awt.event.ActionListener() {

                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    prevScrollAdjustmentValueChanged(null);
                }
            });

            preProcessBox.setText("Pre-Process Stack");
            preProcessBox.setSelected(preProcess);

            colocalBox.setText("Show Co-Localised Trajectories Only");
            colocalBox.setSelected(colocal);
            colocalBox.setEnabled(!monoChrome);

            intensPlotBox.setText("Show Intensity Plots");
            intensPlotBox.setSelected(intensPlot);

            trajPlotBox.setText("Show Trajectory Plots");
            trajPlotBox.setSelected(trajPlot);

            prevResBox.setText("Preview Results");
            prevResBox.setSelected(prevRes);

            analyseButton.setText("Analyse");
            analyseButton.addActionListener(new java.awt.event.ActionListener() {

                public void actionPerformed(java.awt.event.ActionEvent evt) {
                    analyseButtonActionPerformed(evt);
                }
            });

            /*
             * org.jdesktop.layout.GroupLayout layout = new
             * org.jdesktop.layout.GroupLayout(getContentPane());
             * getContentPane().setLayout(layout); layout.setHorizontalGroup(
             * layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().addContainerGap().add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(prevResBox).add(109,
             * 109,
             * 109).add(analyseButton).addContainerGap()).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(trajPlotBox).addContainerGap()).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(intensPlotBox).addContainerGap()).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(colocalBox).addContainerGap()).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(preProcessBox).addContainerGap()).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().add(jLabel15).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 62, Short.MAX_VALUE).add(prevField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(prevScroll,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 0,
             * Short.MAX_VALUE).add(layout.createSequentialGroup().add(jLabel10).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * org.jdesktop.layout.GroupLayout.DEFAULT_SIZE,
             * Short.MAX_VALUE).add(minTrajLengthField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(layout.createSequentialGroup().add(jLabel8).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 41, Short.MAX_VALUE).add(timeResField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(jLabel4).add(layout.createSequentialGroup().add(jLabel6).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 48, Short.MAX_VALUE).add(spatialResField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(layout.createSequentialGroup().add(jLabel3).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 85, Short.MAX_VALUE).add(channel2Combo,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(layout.createSequentialGroup().add(jLabel2).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 85, Short.MAX_VALUE).add(channel1Combo,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(jLabel1).add(org.jdesktop.layout.GroupLayout.TRAILING,
             * virusDiamField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE,
             * 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(org.jdesktop.layout.GroupLayout.TRAILING,
             * layout.createSequentialGroup().add(jLabel12).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 8, Short.MAX_VALUE).add(trajStepTolField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(org.jdesktop.layout.GroupLayout.TRAILING,
             * layout.createSequentialGroup().add(jLabel13).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 12, Short.MAX_VALUE).add(ch1MinPeakField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).add(org.jdesktop.layout.GroupLayout.TRAILING,
             * layout.createSequentialGroup().add(jLabel14).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED,
             * 12, Short.MAX_VALUE).add(ch2MinPeakField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 72,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(jLabel7,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE,
             * org.jdesktop.layout.GroupLayout.DEFAULT_SIZE,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel9).add(jLabel11).add(jLabel5)).add(26,
             * 26,
             * 26)).add(layout.createSequentialGroup().add(prevDetectBox).addContainerGap(175,
             * Short.MAX_VALUE)))))))))); layout.setVerticalGroup(
             * layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING).add(layout.createSequentialGroup().addContainerGap().add(jLabel1).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel2).add(channel1Combo,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel3).add(channel2Combo,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel4).add(virusDiamField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel5)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel6).add(spatialResField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel7,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE,
             * org.jdesktop.layout.GroupLayout.DEFAULT_SIZE,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel8).add(timeResField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel9)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel10).add(minTrajLengthField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel11)).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING).add(layout.createSequentialGroup().addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(trajStepTolField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel12)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(ch1MinPeakField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel13)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(ch2MinPeakField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).add(jLabel14)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(prevDetectBox).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE).add(jLabel15).add(prevField,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 20,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(prevScroll,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE,
             * org.jdesktop.layout.GroupLayout.DEFAULT_SIZE,
             * org.jdesktop.layout.GroupLayout.PREFERRED_SIZE).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(preProcessBox).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(colocalBox).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(intensPlotBox).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(trajPlotBox).addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(prevResBox).addContainerGap(29,
             * Short.MAX_VALUE)).add(layout.createSequentialGroup().addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED).add(analyseButton).addContainerGap()))));
             */
            pack();
        }// </editor-fold>

        private void prevDetectBoxActionPerformed(java.awt.event.ActionEvent evt) {
            prevField.setEnabled(prevDetectBox.isSelected());
            prevScroll.setEnabled(prevDetectBox.isSelected());
            prevScroll.setValue(Integer.parseInt(prevField.getText()));
            if (preview != null) {
                preview.close();
            }
            if (prevDetectBox.isSelected()) {
                preview = new ImagePlus();
                setVariables();
                initialise();
                trajectories.clear();
            }
        }

        /**
         * Analyses the {@link ImageStack} specified by
         * <code>stack</code>.
         */
        public final void initialise() {
            //gaussianRadius = 0.139d / spatialRes;// IsoGaussian filter radius set to 139 nm
            calcParticleRadius(); //Virus particles are approximately 350nm in diameter

            width = stack.getWidth();
            height = stack.getHeight();

            if (width > VIS_SIZE || height > VIS_SIZE) {
                scale = 1.0;
            } else {
                scale = VIS_SIZE / Math.max(width, height);
            }
            ColorProcessor processor = new ColorProcessor(width, height);
            if (!monoChrome) {
                processor.setRGB(colocaliser.getPixels(colocaliser.getChannel1(), slice),
                        colocaliser.getPixels(colocaliser.getChannel2(), slice),
                        colocaliser.getPixels(-1, slice));
            }
            viewDetections(findParticles(1.0, false, slice, slice));

            return;
        }

        /**
         * Creates a visualisation of detected particles in a given frame.
         */
        public void viewDetections(ParticleArray detections) {
            width = (int) Math.round(width * scale);
            height = (int) Math.round(height * scale);
            int radius = (int) Math.round(scale * xyPartRad), i;
            ImageProcessor output;
            if (monoChrome) {
                output = (new TypeConverter((stack.getProcessor(slice + 1)).duplicate(), true)).convertToByte().resize(width, height);
            } else {
                output = (new TypeConverter((stack.getProcessor(slice + 1)).duplicate(), true)).convertToRGB().resize(width, height);
            }
            double scaledSR = scale / spatialRes;
            IsoGaussian c1, c2;
            ArrayList particles = detections.getLevel(0);
            Color c1Color = !monoChrome ? getDrawColor(colocaliser.getChannel1()) : Color.white;
            Color c2Color = !monoChrome ? getDrawColor(colocaliser.getChannel2()) : Color.white;
            for (i = 0; i < particles.size(); i++) {
                c1 = ((IsoGaussian[]) particles.get(i))[0];
                c2 = ((IsoGaussian[]) particles.get(i))[1];
                drawDetections(output, (int) (Math.round(scaledSR * c1.getX())),
                        (int) (Math.round(scaledSR * c1.getY())), radius,
                        (c1.getFit() > 0), c1Color);
                if (c2 != null) {
                    drawDetections(output, (int) (Math.round(scaledSR * c2.getX())),
                            (int) (Math.round(scaledSR * c2.getY())), radius,
                            (c2.getFit() > 0), c2Color);
                }
            }
            preview.setProcessor("Preview of Slice " + slice, output);
            preview.updateAndDraw();
            preview.show();

            return;
        }

        public void drawDetections(ImageProcessor image, int x, int y, int radius,
                boolean drawOval, Color colour) {
            image.setColor(colour);
            if (drawOval) {
                image.drawOval((x - radius), (y - radius), 2 * radius, 2 * radius);
            } else {
                image.drawRect((x - radius), (y - radius), 2 * radius, 2 * radius);
            }
            return;
        }

        private void prevScrollAdjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
            if (evt != null) {
                prevField.setText(prevScroll.getValue() + "");
            }
            prevDetectBoxActionPerformed(null);
        }

        private void analyseButtonActionPerformed(java.awt.event.ActionEvent evt) {
            dialogOKed = true;
            if (preview != null) {
                preview.close();
            }
            setVariables();
            dispose();
        }

        public boolean wasOKed() {
            return dialogOKed;
        }

        public boolean isValidEntries() {
            return validEntries;
        }

        void setVariables() {
            try {
                colocaliser.setChannel1(channel1Combo.getSelectedIndex());
                colocaliser.setChannel2(channel2Combo.getSelectedIndex());
                //virusDiameter = Double.parseDouble(virusDiamField.getText());
                spatialRes = Double.parseDouble(spatialResField.getText()) / 1000.0;
                timeRes = 1.0 / Double.parseDouble(timeResField.getText());
                minTrajLength = Double.parseDouble(minTrajLengthField.getText()) / timeRes;
                trajMaxStep = Double.parseDouble(trajStepTolField.getText());
                chan1MaxThresh = Double.parseDouble(ch1MinPeakField.getText());
                chan2MaxThresh = Double.parseDouble(ch2MinPeakField.getText());
                slice = Integer.parseInt(prevField.getText());
                preProcess = preProcessBox.isSelected();
                colocal = colocalBox.isSelected();
                intensPlot = intensPlotBox.isSelected();
                trajPlot = trajPlotBox.isSelected();
                prevRes = prevResBox.isSelected();
            } catch (NumberFormatException e) {
                Toolkit.getDefaultToolkit().beep();
                IJ.error("Entries in text fields must be numeric.");
                validEntries = false;
            }
        }
        // Variables declaration - do not modify
        private javax.swing.JButton analyseButton;
        private javax.swing.JTextField ch1MinPeakField;
        private javax.swing.JTextField ch2MinPeakField;
        private javax.swing.JComboBox channel1Combo;
        private javax.swing.JComboBox channel2Combo;
        private javax.swing.JCheckBox colocalBox;
        private javax.swing.JCheckBox intensPlotBox;
        private javax.swing.JLabel jLabel1;
        private javax.swing.JLabel jLabel10;
        private javax.swing.JLabel jLabel11;
        private javax.swing.JLabel jLabel12;
        private javax.swing.JLabel jLabel13;
        private javax.swing.JLabel jLabel14;
        private javax.swing.JLabel jLabel15;
        private javax.swing.JLabel jLabel2;
        private javax.swing.JLabel jLabel3;
        private javax.swing.JLabel jLabel4;
        private javax.swing.JLabel jLabel5;
        private javax.swing.JLabel jLabel6;
        private javax.swing.JLabel jLabel7;
        private javax.swing.JLabel jLabel8;
        private javax.swing.JLabel jLabel9;
        private javax.swing.JTextField minTrajLengthField;
        private javax.swing.JCheckBox preProcessBox;
        private javax.swing.JCheckBox prevDetectBox;
        private javax.swing.JTextField prevField;
        private javax.swing.JCheckBox prevResBox;
        private java.awt.Scrollbar prevScroll;
        private javax.swing.JTextField spatialResField;
        private javax.swing.JTextField timeResField;
        private javax.swing.JCheckBox trajPlotBox;
        private javax.swing.JTextField trajStepTolField;
        private javax.swing.JTextField virusDiamField;
        // End of variables declaration
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
            preview = mapTrajectories((int) ((Scrollbar) (dialog.getSliders().get(0))).getValue());
            preview.show();
            return true;
        }

        public ImagePlus mapTrajectories(int trajNumber) {
            ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(trajNumber));
            stack = imp.getStack();
            Rectangle bounds = traj.getBounds();
            int visScale = 20, border = 2 * xyPartRad;
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

            type = traj.getType();
            current = traj.getEnd();
            currentX = lastX = (int) (Math.round(current.getX() / spatialRes - xOffset)) * visScale;
            currentY = lastY = (int) (Math.round(current.getY() / spatialRes - yOffset)) * visScale;
            lastTP = current.getTimePoint() / timeRes;
            current = current.getLink();
            while (current != null) {
                for (j = frames - 1; j >= lastTP; j--) {
                    if (j <= lastTP + 5.0 / timeRes) {
                        processor = outputStack.getProcessor(j + 1);
                        if (!monoChrome) {
                            if (type == ParticleTrajectory.NON_COLOCAL) {
                                processor.setColor(getDrawColor(colocaliser.getChannel1()));
                            } else {
                                processor.setColor(getDrawColor(colocaliser.getChannel2()));
                            }
                        } else {
                            processor.setColor(255);
                        }
                        currentX = (int) (Math.round(current.getX() / spatialRes - xOffset)) * visScale;
                        currentY = (int) (Math.round(current.getY() / spatialRes - yOffset)) * visScale;
                        processor.drawLine(currentX, currentY, lastX, lastY);
                    }
                }
                lastX = currentX;
                lastY = currentY;
                lastTP = current.getTimePoint() / timeRes;
                current = current.getLink();
            }
            return new ImagePlus("Trajectory Number " + (trajNumber + 1), outputStack);
        }
    }
}
