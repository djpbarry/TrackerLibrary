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
    private static double spatialRes = 0.133333;
    private static double timeRes = 1.0;
    private static double trajMaxStep = 2.5;
    private static double minTrajLength = 30.0;
    private static double minTrajDist = 0.0;
    private static double chan1MaxThresh = 0.99;
    private static double chan2MaxThresh = 0.99;
    private static double c1CurveFitTol = 0.5d;
    private static double c2CurveFitTol = 0.0d;
    private static double trackLength = 5.0;
    private static double msdThresh = 0.0;
    private static int nMax = 1;
    private static boolean colocal = true, preProcess = true, gpu = true, prevRes = true,useCals = true, extractsigs = true;
    public static final String[] channels = {"Red", "Green"};
    private static int c1Index = RED;
    private static int c2Index = GREEN;
    public static final int FOREGROUND = 255; //Integer value of foreground pixels

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

    public static double getC1CurveFitTol() {
        return c1CurveFitTol;
    }

    public static void setCurveFitTol(double curveFitTol) {
        UserVariables.c1CurveFitTol = curveFitTol;
    }

    public static double getC2CurveFitTol() {
        return c2CurveFitTol;
    }

    public static void setC2CurveFitTol(double c2CurveFitTol) {
        UserVariables.c2CurveFitTol = c2CurveFitTol;
    }

    public static int getnMax() {
        return nMax;
    }

    public static void setnMax(int nMax) {
        UserVariables.nMax = nMax;
    }

    public static boolean isGpu() {
        return gpu;
    }

    public static void setGpu(boolean gpu) {
        UserVariables.gpu = gpu;
    }

    public static double getMinTrajDist() {
        return minTrajDist;
    }

    public static void setMinTrajDist(double minTrajDist) {
        UserVariables.minTrajDist = minTrajDist;
    }

    public static double getTrackLength() {
        return trackLength;
    }

    public static void setTrackLength(double trackLength) {
        UserVariables.trackLength = trackLength;
    }

    public static boolean isPrevRes() {
        return prevRes;
    }

    public static void setPrevRes(boolean prevRes) {
        UserVariables.prevRes = prevRes;
    }

    public static boolean isUseCals() {
        return useCals;
    }

    public static void setUseCals(boolean useCals) {
        UserVariables.useCals = useCals;
    }

    public static boolean isExtractsigs() {
        return extractsigs;
    }

    public static void setExtractsigs(boolean extractsigs) {
        UserVariables.extractsigs = extractsigs;
    }

    public static double getMsdThresh() {
        return msdThresh;
    }

    public static void setMsdThresh(double msdThresh) {
        UserVariables.msdThresh = msdThresh;
    }

}
