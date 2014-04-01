/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class UserVariables {

    public static final int RED = 0, GREEN = 1, BLUE = 2;
    private static double spatialRes = 0.212;
    private static double timeRes = 1.0;
    private static double trajMaxStep = 2.5;
    private static double minTrajLength = 0.0;
    private static double chan1MaxThresh = 100.0;
    private static double chan2MaxThresh = 0.0;
    private static double curveFitTol = 0.8d;
    private static boolean colocal = false, preProcess = true;
    private static final String[] channels = {"Red", "Green", "Blue"};
    private static int c1Index = RED;
    private static int c2Index = GREEN;

    public static double getSpatialRes() {
        return spatialRes;
    }

    public static void setSpatialRes(double spatialRes) {
        UserVariables.spatialRes = spatialRes;
    }

    public static double getTimeRes() {
        return timeRes;
    }

    public static void setTimeRes(double timeRes) {
        UserVariables.timeRes = timeRes;
    }

    public static double getTrajMaxStep() {
        return trajMaxStep;
    }

    public static void setTrajMaxStep(double trajMaxStep) {
        UserVariables.trajMaxStep = trajMaxStep;
    }

    public static double getMinTrajLength() {
        return minTrajLength;
    }

    public static void setMinTrajLength(double minTrajLength) {
        UserVariables.minTrajLength = minTrajLength;
    }

    public static double getChan1MaxThresh() {
        return chan1MaxThresh;
    }

    public static void setChan1MaxThresh(double chan1MaxThresh) {
        UserVariables.chan1MaxThresh = chan1MaxThresh;
    }

    public static double getChan2MaxThresh() {
        return chan2MaxThresh;
    }

    public static void setChan2MaxThresh(double chan2MaxThresh) {
        UserVariables.chan2MaxThresh = chan2MaxThresh;
    }

    public static boolean isColocal() {
        return colocal;
    }

    public static void setColocal(boolean colocal) {
        UserVariables.colocal = colocal;
    }

    public static boolean isPreProcess() {
        return preProcess;
    }

    public static void setPreProcess(boolean preProcess) {
        UserVariables.preProcess = preProcess;
    }

    public static int getC1Index() {
        return c1Index;
    }

    public static void setC1Index(int c1Index) {
        UserVariables.c1Index = c1Index;
    }

    public static int getC2Index() {
        return c2Index;
    }

    public static void setC2Index(int c2Index) {
        UserVariables.c2Index = c2Index;
    }

    public static double getCurveFitTol() {
        return curveFitTol;
    }

    public static void setCurveFitTol(double curveFitTol) {
        UserVariables.curveFitTol = curveFitTol;
    }

}
