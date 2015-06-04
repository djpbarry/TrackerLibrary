package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.Utils;
import IAClasses.DSPProcessor;
import IAClasses.DataStatistics;
import ij.IJ;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.text.TextWindow;
import java.awt.Color;
import java.awt.Rectangle;
import java.text.DecimalFormat;

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
     * @param t the particle's z-position in a stack or image sequence.
     * @param c1Gaussian {@link IsoGaussian} representation of particle in red
     * channel
     * @param c2Gaussian {@link IsoGaussian} representation of particle in green
     * channel
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
            startTimeIndex = particle.getTimePoint();
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
        if (particle.getC1Gaussian().getMagnitude() > peakIntens) {
            peakIntens = particle.getC1Gaussian().getMagnitude();
            peakTime = particle.getTimePoint() * timeRes;
        }
        return true;
    }

    public boolean checkDetections(Particle particle, double c1tol, double c2tol) {
        if (particle == null) {
            return false;
        }
        boolean c1 = (particle.getC1Gaussian().getFit() >= c1tol);
        boolean c2 = (particle.getC2Gaussian() != null) ? (particle.getC2Gaussian().getFit() >= c2tol) : false;
        if (c1 && c2) {
            dualScore++;
        }
        return addPoint(particle);
    }

    /**
     * Add a new temporary {@link Particle}, which may or may not belong to this
     * trajectory.
     *
     * @param t the particle's z-position in a stack or image sequence.
     * @param redGaussian the detected {@link IsoGaussian} in the red channel
     * representing this object.
     * @param c2Gaussian the detected {@link IsoGaussian} in the green channel
     * representing this object.
     * @param score a measure of the likelihood that this particle belongs to
     * the current trajectory. Lower values indicate greater likelihood.
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

    public double getDuration() {
        int duration = 0;
        Particle current = end;
        if (current == null) {
            return 0.0;
        }
        while (current.getLink() != null) {
            duration += (current.getTimePoint() - (current.getLink()).getTimePoint());
            current = current.getLink();
        }
        return duration / timeRes;
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
    public void printTrajectory(int number, TextWindow output, DecimalFormat formatter, String title) {
        if (output == null) {
            output = new TextWindow(title + " Results",
                    "X\tY\tFrame\tChannel 1\tChannel 2\tChannel 2 " + '\u03C3'
                    + "x\tChannel 2 " + '\u03C3' + "y\t" + '\u03B8',new String(), 1000, 500);
            output.setVisible(true);
        }
        if (formatter == null) {
            formatter = new DecimalFormat("0.000");
        }
        output.append("Particle " + number + "\n");
        Particle current = end;
        while (current != null) {
            double xsig, ysig, theta;
            if (current.getC2Gaussian() != null && current.getC2Gaussian().getFit() > UserVariables.getC1CurveFitTol()) {
                xsig = current.getC2Gaussian().getXSigma();
                ysig = current.getC2Gaussian().getYSigma();
                if (current.getC2Gaussian() instanceof NonIsoGaussian) {
                    theta = ((NonIsoGaussian) current.getC2Gaussian()).getTheta();
                } else {
                    theta = Double.NaN;
                }
            } else {
                xsig = ysig = theta = Double.NaN;
            }
            output.append(formatter.format(current.getX()) + "\t" + formatter.format(current.getY())
                    + "\t" + formatter.format(current.getTimePoint() / timeRes) + "\t"
                    + formatter.format(current.getC1Intensity()) + "\t"
                    + formatter.format(current.getC2Intensity()) + "\t"
                    + formatter.format(xsig) + "\t"
                    + formatter.format(ysig) + "\t"
                    + formatter.format(theta));
            current = current.getLink();
        }
        output.append("\n");
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
            if (current.getC2Gaussian() == null) {
                c2 = 0.0;
            } else {
                c2 = (current.getC2Gaussian()).getMagnitude();
            }
            c1 = (current.getC1Gaussian()).getMagnitude();
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
        if (size < 1) {
            return null;
        }
        double points[][] = new double[2][size];
        Particle current = end;
        for (int i = size - 1; i >= 0; i--) {
            points[0][i] = current.getX();
            points[1][i] = current.getY();
            current = current.getLink();
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
     * <code>D = 1 / (1 + &lambda<sub>1</sub> &lambda<sub>2</sub><sup>-1</sup>)</code>
     * <br> <br> where <code>&lambda<sub>1</sub></code> and
     * <code>&lambda<sub>2</sub></code> are the eigenvalues of the trajectory
     * data and      <code>&lambda<sub>1</sub> <
     * &lambda<sub>2</sub></code>.
     *
     * @param particleNumber the trajectory index.
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
     * <code>getDiffCoeff()</code> method of <ParticleTrajectory</code>.
     *
     * @param traj the trajectory to be analysed.
     * @param label the particle number label to displayed on plots.
     * @param seg the number of time-steps the calculation should be limited to.
     * Set to -1 to include all trajectory points.
     * @param showPlot set to true to display a plot of MSD versus time-step,
     * false otherwise.
     */
    public boolean calcMSD(int label, int seg, boolean showPlot, double[] xPoints, double[] yPoints) {
        int i, j, maxLength, maxStepSize;
        double xval, yval;
        if (xPoints == null) {
            return false;
        }
        int length = xPoints.length;
        if (seg > 0) {
            maxLength = seg;
        } else {
            maxLength = length;
        }
        if (maxLength < 2) {
            diffCoeff = 0.0;
            return false;
        }
        maxStepSize = maxLength / 2;
        double timesteps[] = new double[maxStepSize];
        double msd[] = new double[maxStepSize];
        for (i = 0; i < maxStepSize; i++) {
            for (j = 0, msd[i] = 0.0; i + j < maxLength; j++) {
                xval = Math.pow(xPoints[i + j] - xPoints[j], 2.0) / (maxLength + 1.0);
                yval = Math.pow(yPoints[i + j] - yPoints[j], 2.0) / (maxLength + 1.0);
                msd[i] += xval + yval;
            }
            timesteps[i] = i * timeRes;
        }
        if (showPlot) {
            Plot plot = new Plot("Particle " + label + " Mean Square Displacement",
                    "Time (s)", "Mean Square Displacement (" + IJ.micronSymbol + "m^2)", timesteps, msd,
                    (Plot.X_TICKS + Plot.Y_TICKS + Plot.X_NUMBERS + Plot.Y_NUMBERS));
            plot.setLineWidth(2);
            plot.setColor(Color.BLUE);
            plot.draw();
            plot.show();
        }
        CurveFitter fitter = new CurveFitter(timesteps, msd);
        fitter.doFit(CurveFitter.STRAIGHT_LINE);
        diffCoeff = (fitter.getParams())[1] / 4.0;

        return true;
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

    public double getDirectionality() {
        return directionality;
    }

    public void calcFluorSpread() {
        Particle current = end;
        double xspread = 0.0, yspread = 0.0;
        int points = 0;
        while (current != null) {
            if (current.getC2Gaussian() != null) {
                xspread += current.getC2Gaussian().getXSigma();
                yspread += current.getC2Gaussian().getYSigma();
                points++;
            }
            current = current.getLink();
        }
        if (points >= 10) {
            xFluorSpread = xspread / points;
            yFluorSpread = yspread / points;
        } else {
            xFluorSpread = yFluorSpread = Double.NaN;
        }
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
}
