package ParticleTracking;

import IAClasses.IsoGaussian;
import java.util.ArrayList;

/**
 * 2D Gaussian Curve Fitter based on ImageJ's <code>CurveFitter</code>.
 *
 * TODO Establish a constant reference background level for all fits
 *
 * @author David J Barry
 * @version 1.0, JAN 2011
 */
public class MultiGaussFitter {

    public static final int IterFactor = 500;
    double _2sig2, _maxThresh, xyStepSize, magStepSize, bgStepSize;
    int N_MAX = 1, FIT_RADIUS = 3, FIT_SIZE = 7;
    int STEP_TOL = 2000, ITERATIONS = 100;
    double XY_STEP_SIZE = 0.1, MAG_STEP_SIZE = 10.0, BG_STEP_SIZE = 1.0f;
    private int best;
    double xe[][];
    double ye[][];
    double mag[][];
    double bg[][];
    double r[];
    double sigEst;

    public MultiGaussFitter(int N_MAX, int FIT_RADIUS, int FIT_SIZE) {
        this.N_MAX = N_MAX;
        this.FIT_RADIUS = FIT_RADIUS;
        this.FIT_SIZE = FIT_SIZE;
    }

    void fit(double[][] A, double sigEst) {
        this.sigEst = sigEst;
        _2sig2 = 2.0f * sigEst * sigEst;
        xyStepSize = XY_STEP_SIZE / STEP_TOL;
        magStepSize = MAG_STEP_SIZE / STEP_TOL;
        bgStepSize = BG_STEP_SIZE / STEP_TOL;
        xe = new double[N_MAX][N_MAX];
        ye = new double[N_MAX][N_MAX];
        mag = new double[N_MAX][N_MAX];
        bg = new double[N_MAX][N_MAX];
        r = new double[N_MAX];
        initialiseFitting(A, FIT_RADIUS, xe, ye, mag, bg, r);
    }

    void initialiseFitting(double[][] image, int index, double[][] xe, double[][] ye, double[][] mag, double[][] bg, double[] r) {
        centreOfMass(xe, ye, bg, index, image);
        mag[0][0] = image[FIT_RADIUS + 1][FIT_RADIUS + 1];
        doMultiFit(image, index, 0, xe, ye, mag, bg, r);
        for (int n = 1; n < N_MAX; n++) {
            mag[n][0] = 0.0;
            bg[n][0] = 0.0;
            xe[n][0] = 0.0;
            ye[n][0] = 0.0;
            for (int j = 0; j < FIT_SIZE; j++) {
                for (int i = index - FIT_RADIUS; i < index - FIT_RADIUS + FIT_SIZE; i++) {
                    double residual = image[i][j];
                    for (int m = 0; m < n; m++) {
                        xe[n][m] = xe[n - 1][m];
                        ye[n][m] = ye[n - 1][m];
                        mag[n][m] = mag[n - 1][m];
                        bg[n][m] = bg[n - 1][m];
                        residual -= multiEvaluate(xe[n][m], ye[n][m], mag[n][m], bg[n][m], i, j);
                    }
                    if (residual > mag[n][0]) {
                        mag[n][0] = residual;
                        bg[n][0] = 0.0f;
                        xe[n][0] = i;
                        ye[n][0] = j;
                    }
                }
            }
            doMultiFit(image, index, n, xe, ye, mag, bg, r);
        }
        best = -1;
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < N_MAX; i++) {
            if (r[i] > max) {
                max = r[i];
                best = i;
            }
        }
    }

    void centreOfMass(double[][] x, double[][] y, double[][] bg, int index, double[][] image) {
        float xsum = 0.0f;
        float ysum = 0.0f;
        float sum = 0.0f;
        bg[0][0] = Double.MAX_VALUE;
        for (int j = 0; j < FIT_SIZE; j++) {
            for (int i = index - FIT_RADIUS; i <= index + FIT_RADIUS; i++) {
                xsum += i * image[i][j];
                ysum += j * image[i][j];
                sum += image[i][j];
                if (image[i][j] < bg[0][0]) {
                    bg[0][0] = image[i][j];
                }
            }
        }
        x[0][0] = xsum / sum;
        y[0][0] = ysum / sum;
    }

    void doMultiFit(double[][] M, int x0, int N, double xe[][], double[][] ye, double[][] mag, double[][] bg, double[] r) {
        for (int i = 0; i < ITERATIONS; i++) {
            for (int j = 0; j <= N; j++) {
                float r1 = sumMultiResiduals(x0, xe, ye, mag, bg, M, -XY_STEP_SIZE, 0.0f, 0.0f, 0.0f, j, N);
                float r2 = sumMultiResiduals(x0, xe, ye, mag, bg, M, XY_STEP_SIZE, 0.0f, 0.0f, 0.0f, j, N);
                float r3 = sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, -XY_STEP_SIZE, 0.0f, 0.0f, j, N);
                float r4 = sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, XY_STEP_SIZE, 0.0f, 0.0f, j, N);
                float r5 = sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, 0.0f, -MAG_STEP_SIZE, 0.0f, j, N);
                float r6 = sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, 0.0f, MAG_STEP_SIZE, 0.0f, j, N);
                float r7 = sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, 0.0f, 0.0f, -BG_STEP_SIZE, j, N);
                float r8 = sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, 0.0f, 0.0f, BG_STEP_SIZE, j, N);
                xe[N][j] -= (r2 - r1) * xyStepSize;
                ye[N][j] -= (r4 - r3) * xyStepSize;
                mag[N][j] -= (r6 - r5) * magStepSize;
                bg[N][j] -= (r8 - r7) * bgStepSize;
                if (mag[N][j] < 0.0f) {
                    mag[N][j] = 0.0f;
                }
                if (bg[N][j] < 0.0f) {
                    bg[N][j] = 0.0f;
                }
                if (bg[N][j] > mag[N][j]) {
                    bg[N][j] = mag[N][j];
                }
            }
        }
        r[N] = getRSquared(x0, sumMultiResiduals(x0, xe, ye, mag, bg, M, 0.0f, 0.0f, 0.0f, 0.0f, 0, N), M);
    }

    float sumMultiResiduals(int x0, double[][] xe, double[][] ye, double[][] mag, double[][] bg, double[][] M, double xinc, double yinc, double minc, double bginc, int index, int N) {
        float residuals = 0.0f;
        for (int j = 0; j < FIT_SIZE; j++) {
            for (int i = x0 - FIT_RADIUS; i <= x0 + FIT_RADIUS; i++) {
                float res = 0.0f;
                int k;
                for (k = 0; k < index; k++) {
                    res += multiEvaluate(xe[N][k], ye[N][k], mag[N][k], bg[N][k], i, j);
                }
                res += multiEvaluate(xe[N][k] + xinc, ye[N][k] + yinc, mag[N][k] + minc, bg[N][k] + bginc, i, j);
                for (k = index + 1; k <= N; k++) {
                    res += multiEvaluate(xe[N][k], ye[N][k], mag[N][k], bg[N][k], i, j);
                }
                double e = res - M[i][j];
                residuals += e * e;
            }
        }
        return residuals;
    }

    double getRSquared(int x0, double srs, double[][] M) {
        int y0 = FIT_RADIUS;
        double sumZ = 0.0;
        for (int y = y0 - FIT_RADIUS; y <= y0 + FIT_RADIUS; y++) {
            for (int x = x0 - FIT_RADIUS; x <= x0 + FIT_RADIUS; x++) {
                sumZ += M[x][y];
            }
        }
        double mean = sumZ / FIT_SIZE;
        double sumMeanDiffSqr = 0.0;
        for (int y = y0 - FIT_RADIUS; y <= y0 + FIT_RADIUS; y++) {
            for (int x = x0 - FIT_RADIUS; x <= x0 + FIT_RADIUS; x++) {
                sumMeanDiffSqr += (M[x][y] - mean) * (M[x][y] - mean);
            }
        }
        return 1.0 - (srs / sumMeanDiffSqr);
    }

    double multiEvaluate(double x0, double y0, double mag, double bg, int x, int y) {
        return mag * Math.exp(-(((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (_2sig2))) + bg;
    }

    public ArrayList<IsoGaussian> getFits(double spatialRes, int xoffset, int yoffset, double magThresh) {
        if (best < 0) {
            return null;
        }
        ArrayList<IsoGaussian> fits = new ArrayList();
        for (int i = 0; i <= best; i++) {
            if (mag[best][i] > magThresh) {
                fits.add(new IsoGaussian((xe[best][i] + xoffset) * spatialRes,
                        (ye[best][i] + yoffset) * spatialRes, mag[best][i],
                        sigEst, sigEst, r[best]));
            }
        }
        return fits;
    }
}
