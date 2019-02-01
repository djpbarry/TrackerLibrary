package ParticleTracking;

import Particle.Particle;
import IAClasses.Utils;
import IAClasses.DSPProcessor;
import IAClasses.DataStatistics;
import fiji.plugin.trackmate.Spot;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Rectangle;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Represents a the trajectory followed by a particle through a series of
 * images. It is, in essence, a linked list of individual
 * <code>Particle</code>s.
 *
 * @author David J Barry
 * @version 1.0, JAN 2011
 */
public class ParticleTrajectory {

    protected Particle end = null, temp = null;
    private int size = 0, dualScore = 0, tempRow = -1, tempColumn = -1;
    public final static int NON_COLOCAL = 0, UNKNOWN = 1, COLOCAL = 2; //Flags
    protected double tempScore = Double.MAX_VALUE, xVelocity = 0.0, yVelocity = 0.0, projectXVel,
            projectYVel, diffCoeff, boxCountFD = 0.0, angleSpread = 0.0,
            stepSpread, timeRes, specFD, meanKappa, logDC, directionality,
            peakIntens, peakTime;
    private double[] kappa;
    private double[] smoothXPoints, smoothYPoints;
    private Rectangle bounds;
    public static double scale;
    private static final int segment = 5;
    private double xFluorSpread, yFluorSpread;
    private int startTimeIndex;
    private static Plot msdPlot;
    private static String plotLegend = "";
    private final int MIN_POINTS_TO_AVERAGE = 10;
    private final float D_SCALING = 4.0f;
    private static ArrayList<DescriptiveStatistics> globalMSD = new ArrayList<>();

    public ParticleTrajectory() {
    }

    /**
     * Constructs a new empty trajectory.
     */
    public ParticleTrajectory(double timeRes, double spatRes) {
        this.timeRes = timeRes;
        scale = 1.0 / spatRes;
        peakIntens = 0.0;
    }

    /**
     * Add a new {@link Particle} to this trajectory.
     *
     */
    public boolean addPoint(Particle particle) {
        if (particle == null) {
            return false;
        }
        particle.setLink(end);
        end = particle;
        int newX = (int) Math.round(scale * end.getX());
        int newY = (int) Math.round(scale * end.getY());
        if (size < 1) {
            startTimeIndex = particle.getFrameNumber();
            bounds = new Rectangle(newX, newY, 0, 0);
        }
        if (newX < bounds.x) {
            bounds = new Rectangle(newX, bounds.y, bounds.width + bounds.x - newX, bounds.height);
        } else if (newX > bounds.x + bounds.width) {
            bounds = new Rectangle(bounds.x, bounds.y, newX - bounds.x, bounds.height);
        }
        if (newY < bounds.y) {
            bounds = new Rectangle(bounds.x, newY, bounds.width, bounds.height + bounds.y - newY);
        } else if (newY > bounds.y + bounds.height) {
            bounds = new Rectangle(bounds.x, bounds.y, bounds.width, newY - bounds.y);
        }
        size++;
        temp = null;
        tempScore = Double.MAX_VALUE;
        tempRow = tempColumn = -1;
        if (size > segment) {
            updateVelocity();
        }
        if (particle.getMagnitude() > peakIntens) {
            peakIntens = particle.getMagnitude();
            peakTime = particle.getFrameNumber() * timeRes;
        }
        return true;
    }

    public boolean checkDetections(Particle particle, double tol) {
        if (particle == null) {
            return false;
        }
        boolean c1 = particle.getFeature(Spot.QUALITY) >= tol;
        boolean c2 = particle.getColocalisedParticle() != null;
        if (c1 && c2) {
            dualScore++;
        }
        return addPoint(particle);
    }

    /**
     * Add a new temporary {@link Particle}, which may or may not belong to this
     * trajectory.
     *
     */
    public boolean addTempPoint(Particle particle, double score, int row, int column) {
        particle.setLink(end);
        temp = particle;
        if (temp == null) {
            return false;
        }
        tempScore = score;
        tempRow = row;
        tempColumn = column;

        return true;
    }

    /**
     * Returns the number of <code>Particle</code>s in this trajectory.
     */
    public int getSize() {
        return size;
    }

    /**
     * Returns the total distance travelled on this trajectory.
     */
    public double getDisplacement(Particle start, int steps) {
        double displacement = 0.0;
        Particle current = start;
        if (current == null) {
            return 0.0;
        }
        int s = 0;
        while (current.getLink() != null && s < steps) {
            displacement += Utils.calcDistance(current.getX(), current.getY(),
                    (current.getLink()).getX(), (current.getLink()).getY());
            current = current.getLink();
            s++;
        }
        return displacement;
    }

    public int getNumberOfFrames() {
        Particle current = end;
        if (current == null) {
            return 0;
        }
        return current.getFrameNumber() - getStart().getFrameNumber();
    }

    /**
     * Returns a non-zero score for the temporary particle, or
     * <code>Double.MAX_VALUE</code> if no such particle exists.
     */
    public double getTempScore() {
        return tempScore;
    }

    /**
     * Returns the most recently-added <code>Particle</code> in this trajectory.
     */
    public Particle getEnd() {
        return end;
    }

    /**
     * Returns the most recently-added temporary <code>Particle</code> in this
     * trajectory.
     */
    public Particle getTemp() {
        return temp;
    }

    /**
     * Prints the details of this trajectory to the output display.
     */
    public void printTrajectory(int number, TextWindow output, DecimalFormat formatter, String title, ImagePlus[] inputs) {
        if (output == null || inputs == null) {
            return;
        }
        if (formatter == null) {
            formatter = new DecimalFormat("0.000");
        }
        ImageStack[] stacks = new ImageStack[2];
        stacks[0] = inputs[0].getImageStack();
        if (inputs[1] != null) {
            stacks[1] = inputs[1].getImageStack();
        } else {
            stacks[1] = null;
        }
        Particle current = end;
        while (current != null) {
            int frame = (int) Math.round(current.getFeature(Spot.FRAME));
            double x = current.getX();
            double y = current.getY();
            int xPix = (int) Math.round(x / UserVariables.getSpatialRes());
            int yPix = (int) Math.round(y / UserVariables.getSpatialRes());
            double c1Mag = stacks[0].getProcessor(frame + 1).getPixelValue(xPix, yPix);
            double c2Mag = Double.NaN;
            if (stacks[1] != null) {
                c2Mag = stacks[1].getProcessor(frame + 1).getPixelValue(xPix, yPix);
            }
            output.append(number
                    + "\t" + frame
                    + "\t" + formatter.format(frame / UserVariables.getTimeRes())
                    + "\t" + formatter.format(x)
                    + "\t" + formatter.format(y)
                    + "\t" + formatter.format(c1Mag)
                    + "\t" + formatter.format(c2Mag));
            current = current.getLink();
        }
    }

    public int getType(double thresh) {
        if (size < 1) {
            return UNKNOWN;
        }
        if ((double) dualScore / size > thresh) {
            return COLOCAL;
        } else {
            return NON_COLOCAL;
        }
    }

    void updateVelocity() {
        if (size < segment) {
            xVelocity = yVelocity = 0.0;
            return;
        }
        int i;
        Particle current = end;
        double xDist = 0.0, yDist = 0.0, x = current.getX(), y = current.getY();
        for (i = 0; i < segment && i < size - 1; i++) {
            current = current.getLink();
            xDist += x - current.getX();
            yDist += y - current.getY();
            x = current.getX();
            y = current.getY();
        }
        xVelocity = xDist / i;
        yVelocity = yDist / i;
    }

    void projectVelocity(double testX, double testY) {
        if (size < segment) {
            projectXVel = projectYVel = 0.0;
            return;
        }
        int i;
        Particle current = end;
        double xDist = testX - current.getX();
        double yDist = testY - current.getY();
        double x = current.getX();
        double y = current.getY();
        for (i = 1; i < segment && i < size - 1; i++) {
            try {
                current = current.getLink();
                xDist += x - current.getX();
                yDist += y - current.getY();
            } catch (Exception e) {
                e.toString();
            }
            x = current.getX();
            y = current.getY();
        }
        projectXVel = xDist / i;
        projectYVel = yDist / i;
    }

    public double getFluorRatio() {
        Particle current = end;
        double c1, c2, ratio = 0.0;
        int points = 0;

        while (current != null) {
            if (current.getColocalisedParticle() == null) {
                c2 = 0.0;
            } else {
                c2 = current.getMagnitude();
            }
            c1 = current.getMagnitude();
            ratio += c2 / c1;
            points++;
            current = current.getLink();
        }

        if (points > 0) {
            return ratio / points;
        } else {
            return Double.NaN;
        }
    }

    public double getYVelocity() {
        return yVelocity;
    }

    public double getXVelocity() {
        return xVelocity;
    }

    public double getProjectXVel() {
        return projectXVel;
    }

    public double getProjectYVel() {
        return projectYVel;
    }

    public double getDiffCoeff() {
        return diffCoeff;
    }

    public void setFracDim(double fracDim) {
        this.boxCountFD = fracDim;
    }

    public int getTempColumn() {
        return tempColumn;
    }

    public int getTempRow() {
        return tempRow;
    }

    public double getAngleSpread() {
        return angleSpread;
    }

    public int getDualScore() {
        return dualScore;
    }

    public double getStepSpread() {
        return stepSpread;
    }

    public Rectangle getBounds() {
        return bounds;
    }

    public double[][] getPoints() {
        int length = getNumberOfFrames();
        if (length < 1) {
            return null;
        }
        double points[][] = new double[3][length];
        Particle current = end;
        int lastFrameNumber;
        int interpolate = 1;
        for (int i = length - 1; i >= 0; i--) {
            points[0][i] = current.getX();
            points[1][i] = current.getY();
            points[2][i] = current.getFrameNumber();
            lastFrameNumber = current.getFrameNumber();
            if (current.getLink() != null && lastFrameNumber - current.getLink().getFrameNumber() > interpolate) {
                interpolate++;
            } else {
                interpolate = 1;
                current = current.getLink();
            }
        }
        return points;
    }

    public void smooth() {
        double sigma = 2.5 * timeRes;
        int gLength = (int) Math.round(12.5 * timeRes);
        double input[][] = getPoints();
        if (input == null) {
            return;
        }
        double gaussian[] = Utils.generateGaussian(sigma, gLength);
        double dgaussian[] = Utils.generateGaussianFirstDeriv(sigma, gLength);
        double ddgaussian[] = Utils.generateGaussianSecondDeriv(sigma, gLength);
        double X[] = Utils.convolve(gaussian, input[0]);
        double Y[] = Utils.convolve(gaussian, input[1]);
        if (X.length <= 2 * (gLength - 1)) {
            smoothXPoints = null;
            smoothYPoints = null;
            return;
        }
        double dX[] = Utils.convolve(dgaussian, input[0]);
        double dY[] = Utils.convolve(dgaussian, input[1]);
        double ddX[] = Utils.convolve(ddgaussian, input[0]);
        double ddY[] = Utils.convolve(ddgaussian, input[1]);
        double shrinkX = (1.0 / ddX[0])
                * (1.0 - Math.exp(-sigma * sigma * ddX[0] * ddX[0] / 2.0));
        double shrinkY = (1.0 / ddY[0])
                * (1.0 - Math.exp(-sigma * sigma * ddY[0] * ddY[0] / 2.0));
        smoothXPoints = new double[X.length - 2 * (gLength - 1)];
        smoothYPoints = new double[Y.length - 2 * (gLength - 1)];
        kappa = new double[Y.length - 2 * (gLength - 1)];
        for (int i = gLength - 1; i < X.length - gLength + 1; i++) {
            smoothXPoints[i - gLength + 1] = X[i] - shrinkX;
            smoothYPoints[i - gLength + 1] = Y[i] - shrinkY;
            kappa[i - gLength + 1] = (dX[i] * ddY[i] - dY[i] * ddX[i])
                    / Math.pow(dX[i] * dX[i] + dY[i] * dY[i], 1.5);
        }
        DataStatistics stats = new DataStatistics(95, kappa, kappa.length);
        meanKappa = Math.abs(1000.0 * stats.getMean());
        calcSpec();
    }

    public boolean calcSpec() {
        if (kappa == null) {
            return false;
        }
        int kSize = kappa.length;
        double[] upsampledK = DSPProcessor.upScale(kappa);
        double sampleRate = (1.0 / timeRes) * upsampledK.length / kSize;
        double fSpec[] = DSPProcessor.calcFourierSpec(upsampledK, sampleRate);
        if (fSpec == null) {
            return false;
        }
        logDC = Math.log(fSpec[0]);
        return true;
    }

    public double getBoxCountFD() {
        return boxCountFD;
    }

    public double getLogDC() {
        return logDC;
    }

    public double getMeanKappa() {
        return meanKappa;
    }

    /**
     * Calculates the directionality of the trajectory specified by
     * <code>particleNumber</code>, where directionality ( <code>D</code>) is
     * calculated according to: <br> <br>
     * <code>D = 1 / (1 + &lambda;<sub>1</sub>
     * &lambda;<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where <code>&lambda;<sub>1</sub></code> and
     * <code>&lambda;<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and      <code>&lambda;<sub>1</sub> &lambda;<sub>2</sub></code>.
     *
     * @return the directionality of the specified trajectory.
     */
    public boolean calcDirectionality(double[] xPoints, double[] yPoints) {
        if (xPoints == null) {
            return false;
        }
        int length = xPoints.length;
        double xSum = 0, ySum = 0;
        for (int i = 0; i < length; i++) {
            xSum += xPoints[i];
            ySum += yPoints[i];
        }
        double[] eigenvalues = Utils.calcEigenvalues(
                Utils.covarianceMatrix(xPoints, yPoints, xSum, ySum));
        if (Math.abs(eigenvalues[0]) > Math.abs(eigenvalues[1])) {
            directionality = 1.0d / (1.0d + 1.0d / Math.sqrt(Math.abs(eigenvalues[0] / eigenvalues[1])));
        } else {
            directionality = 1.0d / (1.0d + 1.0d / Math.sqrt(Math.abs(eigenvalues[1] / eigenvalues[0])));
        }
        return true;
    }

    /**
     * Calculates the Mean Square Displacement (MSD) of the given
     * {@link ParticleTrajectory}. The result can be accessed via the
     * <code>getDiffCoeff()</code> method of <code>ParticleTrajectory</code>.
     */
    public boolean calcMSD(int seg, int label) {
        int maxLength;
        double xval, yval;
        double points[][] = getPoints();
        double xPoints[] = points[0], yPoints[] = points[1];
        if (xPoints == null) {
            return false;
        }
        int length = xPoints.length;
        if (seg > 0) {
            maxLength = seg;
        } else {
            maxLength = length;
        }
        ArrayList<Double> timesteps = new ArrayList<>();
        ArrayList<Double> msd = new ArrayList<>();
        for (int i = 0; i < maxLength; i++) {
            DescriptiveStatistics thisMSD = new DescriptiveStatistics();
            for (int j = 0; i + j < maxLength; j++) {
                xval = Math.pow(xPoints[i + j] - xPoints[j], 2.0);
                yval = Math.pow(yPoints[i + j] - yPoints[j], 2.0);
                thisMSD.addValue(xval + yval);
            }
            long N = thisMSD.getN();
            if (N >= MIN_POINTS_TO_AVERAGE) {
                timesteps.add(i / timeRes);
                msd.add(thisMSD.getMean());
            }
        }
        if (!(msd.size() > 0)) {
            return false;
        }
        if (msdPlot == null) {
            msdPlot = new Plot("Mean Square Displacement",
                    "Time (s)", "Mean Square Displacement (" + IJ.micronSymbol + "m^2)");
            msdPlot.setLineWidth(3);
        }
        double[] tsA = new double[timesteps.size()];
        double[] msdA = new double[msd.size()];
        for (int i = 0; i < timesteps.size(); i++) {
            tsA[i] = timesteps.get(i);
        }
        for (int i = 0; i < msd.size(); i++) {
            msdA[i] = msd.get(i);
            if (globalMSD.size() <= i) {
                globalMSD.add(new DescriptiveStatistics());
            }
            globalMSD.get(i).addValue(msd.get(i));
        }
        Random r = new Random();
        msdPlot.setColor(new Color(r.nextFloat(), r.nextFloat(), r.nextFloat()));
        msdPlot.addPoints(timesteps, msd, Plot.CONNECTED_CIRCLES);
        msdPlot.setLimitsToFit(false);
        plotLegend = ((plotLegend.concat("Particle ")).concat(String.valueOf(label))).concat("\n");
        msdPlot.addLegend(plotLegend);
        msdPlot.draw();
        msdPlot.show();
        CurveFitter fitter = new CurveFitter(tsA, msdA);
        fitter.doFit(CurveFitter.STRAIGHT_LINE);
        diffCoeff = (fitter.getParams())[1] / D_SCALING;

        return true;
    }

    public static void drawGlobalMSDPlot() {
        Plot globalMsdPlot = new Plot("Population Mean Square Displacement",
                "Time (s)", "Mean Square Displacement (" + IJ.micronSymbol + "m^2)");
        globalMsdPlot.setLineWidth(3);
        globalMsdPlot.setColor(Color.red);
        int N = globalMSD.size();
        if (!(N > 0)) {
            return;
        }
        double[] tsA = new double[N];
        double[] msdA = new double[N];
        double[] msdErrorA = new double[N];
        for (int i = 0; i < N; i++) {
            DescriptiveStatistics ds = globalMSD.get(i);
            msdErrorA[i] = ds.getStandardDeviation() / Math.sqrt(ds.getN());
            tsA[i] = i / UserVariables.getTimeRes();
            msdA[i] = ds.getMean();
        }
        globalMsdPlot.addPoints(tsA, msdA, Plot.CONNECTED_CIRCLES);
        globalMsdPlot.addErrorBars(msdErrorA);
        globalMsdPlot.setLimitsToFit(false);
        globalMsdPlot.draw();
        globalMsdPlot.show();
    }

    public static Plot getMsdPlot() {
        return msdPlot;
    }

    public static void resetMSDPlot() {
        globalMSD = new ArrayList<>();
        msdPlot = null;
    }

    public boolean calcAngleSpread() {
        if (smoothXPoints == null) {
            return false;
        }
        int length = smoothXPoints.length;
        if (length < 2) {
            return false;
        }
        double angles[] = new double[length - 2];
        double m1, m2;
        for (int i = 1; i < length - 1; i++) {
            m1 = (smoothYPoints[i + 1] - smoothYPoints[i]) / (smoothXPoints[i + 1] - smoothXPoints[i]);
            m2 = (smoothYPoints[i] - smoothYPoints[i - 1]) / (smoothXPoints[i] - smoothXPoints[i - 1]);
            angles[i - 1] = Math.atan(Math.abs((m1 - m2) / (1 + m1 * m2)));
        }
        DataStatistics stats = new DataStatistics(0.0, angles, length - 2);
        angleSpread = stats.getStdDev();
        return true;
    }

    public boolean calcStepSpread() {
        if (smoothXPoints == null || smoothXPoints.length < 3) {
            return false;
        }
        int length = smoothXPoints.length;
        double steps[] = new double[length - 1];
        for (int i = 0; i < length - 1; i++) {
            steps[i] = Utils.calcDistance(smoothXPoints[i], smoothYPoints[i],
                    smoothXPoints[i + 1], smoothYPoints[i + 1]);
        }
        DataStatistics stats = new DataStatistics(0.0, steps, length - 1);
        stepSpread = stats.getStdDev();
        return true;
    }

    public void setDiffCoeff(double diffCoeff) {
        this.diffCoeff = diffCoeff;
    }

    public double getDirectionality() {
        return directionality;
    }

    public double getxFluorSpread() {
        return xFluorSpread;
    }

    public double getyFluorSpread() {
        return yFluorSpread;
    }

    public double getPeakTime() {
        return peakTime;
    }

    public int getStartTimeIndex() {
        return startTimeIndex;
    }

    public Particle getStart() {
        Particle current = end;
        if (current == null) {
            return null;
        }
        while (current.getLink() != null) {
            current = current.getLink();
        }
        return current;
    }

    public void addTrajectory(ParticleTrajectory traj) {
        Particle[] points = new Particle[traj.getSize()];
        Particle current = traj.getEnd();
        int i = points.length - 1;
        while (current != null) {
            points[i--] = current;
            current = current.getLink();
        }
        for (Particle p : points) {
            this.addPoint(p);
        }
    }
}
