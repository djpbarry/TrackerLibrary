package tracking;

/**
 * 2D Gaussian Curve Fitter based on ImageJ's
 * <code>CurveFitter</code>.
 *
 * TODO Establish a constant reference background level for all fits
 *
 * @author David J Barry
 * @version 1.0, JAN 2011
 */
public class MultiGaussFitter {

    public static final int IterFactor = 500;
    private static final double alpha = -1.0;     // reflection coefficient
    private static final double beta = 0.5;   // contraction coefficient
    private static final double gamma = 2.0;      // expansion coefficient
    private static final double root2 = 1.414214; // square root of 2
    private double[] xData, yData;  // x,y data to fit
    private double[][] zData;  // z data to fit
    private int numPoints;          // number of data points
    private int numParams;          // number of parametres
    private int curves;
    private int numVertices;        // numParams+1 (includes sumLocalResiduaalsSqrd)
    private int worst;          // worst current parametre estimates
    private int nextWorst;      // 2nd worst current parametre estimates
    private int best;           // best current parametre estimates
    private double[][] simp;        // the simplex (the last element of the array at each vertice is the sum of the square of the residuals)
    private double[] next;      // new vertex to be tested
    private int numIter;        // number of iterations so far
    private int maxIter;    // maximum number of iterations per restart
    private int restarts;   // number of times to restart simplex after first soln.
    private static int defaultRestarts = 2;  // default number of restarts
    private int nRestarts;  // the number of restarts that occurred
    private static double maxError = 1e-10;    // maximum error tolerance

    /**
     * Construct a new CurveFitter.
     */
    public MultiGaussFitter(double[][] zVals, int peaks) {
        xData = new double[zVals[0].length];
        yData = new double[xData.length];
        this.zData = zVals;
        this.curves = peaks;
        numPoints = xData.length * yData.length;
        for (int i = 0; i < xData.length; i++) {
            xData[i] = i;
            yData[i] = i;
        }
        numParams = 1 + peaks * 4;
    }

    public boolean doFit(double[] estimates) {
        if (estimates == null || xData == null || yData == null || zData == null) {
            return false;
        }
        initialize(estimates);
        restart(0);

        numIter = 0;
        boolean done = false;
        double[] center = new double[numParams];  // mean of simplex vertices
        while (!done) {
            numIter++;
            for (int i = 0; i < numParams; i++) {
                center[i] = 0.0;
            }
            // get mean "center" of vertices, excluding worst
            for (int i = 0; i < numVertices; i++) {
                if (i != worst) {
                    for (int j = 0; j < numParams; j++) {
                        center[j] += simp[i][j];
                    }
                }
            }
            // Reflect worst vertex through centre
            for (int i = 0; i < numParams; i++) {
                center[i] /= numParams;
                next[i] = center[i] + alpha * (simp[worst][i] - center[i]);
            }
            sumResiduals(next);
            // if it's better than the best...
            if (next[numParams] <= simp[best][numParams]) {
                newVertex();
                // try expanding it
                for (int i = 0; i < numParams; i++) {
                    next[i] = center[i] + gamma * (simp[worst][i] - center[i]);
                }
                sumResiduals(next);
                // if this is even better, keep it
                if (next[numParams] <= simp[worst][numParams]) {
                    newVertex();
                }
            } // else if better than the 2nd worst keep it...
            else if (next[numParams] <= simp[nextWorst][numParams]) {
                newVertex();
            } // else try to make positive contraction of the worst
            else {
                for (int i = 0; i < numParams; i++) {
                    next[i] = center[i] + beta * (simp[worst][i] - center[i]);
                }
                sumResiduals(next);
                // if this is better than the second worst, keep it.
                if (next[numParams] <= simp[nextWorst][numParams]) {
                    newVertex();
                } // if all else fails, contract simplex in on best
                else {
                    for (int i = 0; i < numVertices; i++) {
                        if (i != best) {
                            for (int j = 0; j < numVertices; j++) {
                                simp[i][j] = beta * (simp[i][j] + simp[best][j]);
                            }
                            sumResiduals(simp[i]);
                        }
                    }
                }
            }
            order();

            double rtol = 2 * Math.abs(simp[best][numParams] - simp[worst][numParams])
                    / (Math.abs(simp[best][numParams]) + Math.abs(simp[worst][numParams]) + 0.0000000001);

            if (numIter >= maxIter) {
                done = true;
            } else if (rtol < maxError) {
                restarts--;
                if (restarts < 0) {
                    done = true;
                } else {
                    restart(best);
                }
            }
        }
        return true;
    }

    /**
     * Initialise the simplex
     */
    boolean initialize(double[] estimates) {
        if (estimates == null || xData == null || yData == null || zData == null) {
            return false;
        }
        // Calculate some things that might be useful for predicting parametres
        numVertices = numParams + 1;      // need 1 more vertice than parametres,
        simp = new double[numVertices][numVertices];
        next = new double[numVertices];
        maxIter = IterFactor * numParams * numParams;  // TODO Where does this estimate come from?
        restarts = defaultRestarts;
        nRestarts = 0;
        System.arraycopy(estimates, 0, simp[0], 0, estimates.length);

        return true;
    }

    /**
     * Restart the simplex at the nth vertex
     */
    boolean restart(int n) {
        if (simp == null || n >= simp.length) {
            return false;
        }
        // Copy nth vertice of simplex to first vertice
        System.arraycopy(simp[n], 0, simp[0], 0, numParams);
        sumResiduals(simp[0]);          // Get sum of residuals^2 for first vertex
        double[] step = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            step[i] = simp[0][i] / 2.0;     // Step half the parametre value
            if (step[i] == 0.0) // We can't have them all the same or we're going nowhere
            {
                step[i] = 0.01;
            }
        }
        // Some kind of factor for generating new vertices
        double[] p = new double[numParams];
        double[] q = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            p[i] = step[i] * (Math.sqrt(numVertices) + numParams - 1.0) / (numParams * root2);
            q[i] = step[i] * (Math.sqrt(numVertices) - 1.0) / (numParams * root2);
        }
        // Create the other simplex vertices by modifing previous one.
        for (int i = 1; i < numVertices; i++) {
            for (int j = 0; j < numParams; j++) {
                simp[i][j] = simp[i - 1][j] + q[j];
            }
            simp[i][i - 1] = simp[i][i - 1] + p[i - 1];
            sumResiduals(simp[i]);
        }
        // Initialise current lowest/highest parametre estimates to simplex 1
        best = 0;
        worst = 0;
        nextWorst = 0;
        order();
        nRestarts++;
        return true;
    }

    /**
     * Returns 'fit' formula value for parameters "p" at "x"
     */
    public double evaluate(double[] p, double x, double y) {
        if (p == null) {
            return Double.NaN;
        }
        if (curves == 1) {
            return p[0] + p[1] * Math.exp(-(((x - p[2]) * (x - p[2])) + ((y - p[3]) * (y - p[3]))) / (2 * p[4] * p[4]));
        } else if (curves == 2) {
            return p[0] + p[1] * Math.exp(-(((x - p[2]) * (x - p[2])) + ((y - p[3]) * (y - p[3]))) / (2 * p[4] * p[4]))
                    + p[5] * Math.exp(-(((x - p[6]) * (x - p[6])) + ((y - p[7]) * (y - p[7]))) / (2 * p[8] * p[8]));
        } else {
            return 0.0;
        }
    }

    /**
     * Returns residuals array ie. differences between data and curve.
     */
    public double[] getResiduals() {
        if (!(numPoints > 0)) {
            return null;
        }
        double[] params = getParams();
        double[] residuals = new double[numPoints];
        for (int x = 0; x < xData.length; x++) {
            for (int y = 0; y < yData.length; y++) {
                residuals[x * y] = zData[x][y] - evaluate(params, xData[x], yData[y]);
            }
        }
        return residuals;
    }

    /**
     * Returns R<sup>2</sup>, where 1.0 is best.<br> <br> R<sup>2</sup> = 1.0 -
     * SSE/SSD<br> <br> where SSE is the sum of the squares of the errors and
     * SSD is the sum of the squares of the deviations about the mean.
     */
    public double getRSquared() {
        if (numPoints < 1) {
            return Double.NaN;
        }
        double sumZ = 0.0;
        for (int x = 0; x < xData.length; x++) {
            for (int y = 0; y < yData.length; y++) {
                sumZ += zData[x][y];
            }
        }
        double mean = sumZ / numPoints;
        double sumMeanDiffSqr = 0.0;
        for (int x = 0; x < xData.length; x++) {
            for (int y = 0; y < yData.length; y++) {
                sumMeanDiffSqr += Math.pow(zData[x][y] - mean, 2);
            }
        }
        double rSquared = 0.0;
        if (sumMeanDiffSqr > 0.0) {
            double srs = getSumResidualsSqr();
            rSquared = 1.0 - srs / sumMeanDiffSqr;
        }
        return rSquared;
    }

    /**
     * Adds sum of square of residuals to end of array of parameters
     */
    boolean sumResiduals(double[] x) {
        if (x == null) {
            return false;
        }
        double e;
        x[numParams] = 0.0;
        for (int i = 0; i < xData.length; i++) {
            for (int j = 0; j < yData.length; j++) {
                e = evaluate(x, xData[i], yData[j]) - zData[i][j];
                x[numParams] = x[numParams] + (e * e);
            }
        }
        return true;
    }

    /**
     * Keep the "next" vertex
     */
    boolean newVertex() {
        if (next == null) {
            return false;
        }
        System.arraycopy(next, 0, simp[worst], 0, numVertices);
        return true;
    }

    /**
     * Find the worst, nextWorst and best current set of parameter estimates
     */
    void order() {
        for (int i = 0; i < numVertices; i++) {
            if (simp[i][numParams] < simp[best][numParams]) {
                best = i;
            }
            if (simp[i][numParams] > simp[worst][numParams]) {
                worst = i;
            }
        }
        nextWorst = best;
        for (int i = 0; i < numVertices; i++) {
            if (i != worst) {
                if (simp[i][numParams] > simp[nextWorst][numParams]) {
                    nextWorst = i;
                }
            }
        }
    }

    /**
     * Get the set of parameter values from the best corner of the simplex
     */
    public double[] getParams() {
        order();
        if (simp != null) {
            return simp[best];
        } else {
            return null;
        }
    }

    /*
     * Last "parametre" at each vertex of simplex is sum of residuals for the
     * curve described by that vertex
     */
    public double getSumResidualsSqr() {
        double[] params = getParams();
        if (params != null) {
            return params[getNumParams()];
        } else {
            return Double.NaN;
        }
    }

    public int getNumParams() {
        return numParams;
    }
}
