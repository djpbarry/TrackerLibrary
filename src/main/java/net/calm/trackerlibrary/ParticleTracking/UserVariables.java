/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.trackerlibrary.ParticleTracking;

import ij.process.AutoThresholder;

/**
 * Runtime configuration for the deterministic tracking pipeline.
 *
 * <p>This class is an instance holder for tracking settings. The mutable state
 * lives on a single process-wide {@linkplain #getInstance() instance} that can
 * be swapped for a fresh one via {@link #setInstance(UserVariables)} (useful for
 * tests and for isolating one run from another). The static getter/setter
 * methods delegate to that instance purely for backward compatibility with
 * callers that predate the instance-based configuration; new code should obtain
 * (and pass) an explicit {@link UserVariables} instance through the tracking
 * pipeline instead of relying on the process-wide default.</p>
 */
public class UserVariables {

    public static final int RED = 0, GREEN = 1, BLUE = 2;
    public static final int MAXIMA = 3, BLOBS = 4, GAUSS = 5;
    public static final int RANDOM = 6, DIRECTED = 7;
    public static final int FOREGROUND = 255; //Integer value of foreground pixels

    private static UserVariables instance;

    private double spatialRes = 0.1;
    private double timeRes = 1.0;
    private double trajMaxStep = 0.75;
    private double minTrajLength = 10.0;
    private double minTrajDist = 0.0;
    private double curveFitTol = 0.5d;
    private double blobSize = 0.5;
    private double blobThresh = 0.1;
    private double trackLength = 5.0;
    private double msdThresh = 0.0;
    private int nMax = 1;
    private double colocalThresh = 0.25;
    private boolean colocal = true, preProcess = true, gpu = false, useCals = false, extractsigs = false;
    private double sigEstGreen = 0.2;
    private double sigEstRed = 0.3;
    private int minMSDPoints = 10;
    private boolean fitC2 = false, trackRegions = false;
    private int detectionMode = MAXIMA;
    private double filterRadius = 0.133;
    private int motionModel = RANDOM;
    private int maxFrameGap = 3;
    private String c1ThreshMethod = AutoThresholder.Method.Li.toString();
    private String c2ThreshMethod = AutoThresholder.Method.Li.toString();

    public UserVariables() {
    }

    /**
     * Returns the process-wide default instance, creating it lazily.
     */
    public static synchronized UserVariables getInstance() {
        if (instance == null) {
            instance = new UserVariables();
        }
        return instance;
    }

    /**
     * Replaces the process-wide instance. Pass {@code null} to reset it to a
     * fresh default on the next {@link #getInstance()} call.
     */
    public static synchronized void setInstance(UserVariables newInstance) {
        instance = newInstance;
    }

    public double getSpatialRes() {
        return spatialRes;
    }

    public void setSpatialRes(double spatialRes) {
        this.spatialRes = spatialRes;
    }

    public double getTimeRes() {
        return timeRes;
    }

    public void setTimeRes(double timeRes) {
        this.timeRes = timeRes;
    }

    public double getTrajMaxStep() {
        return trajMaxStep;
    }

    public void setTrajMaxStep(double trajMaxStep) {
        this.trajMaxStep = trajMaxStep;
    }

    public double getMinTrajLength() {
        return minTrajLength;
    }

    public void setMinTrajLength(double minTrajLength) {
        this.minTrajLength = minTrajLength;
    }

    public String getC1ThreshMethod() {
        return c1ThreshMethod;
    }

    public void setC1ThreshMethod(String c1ThreshMethod) {
        this.c1ThreshMethod = c1ThreshMethod;
    }

    public String getC2ThreshMethod() {
        return c2ThreshMethod;
    }

    public void setC2ThreshMethod(String c2ThreshMethod) {
        this.c2ThreshMethod = c2ThreshMethod;
    }

    public boolean isColocal() {
        return colocal;
    }

    public void setColocal(boolean colocal) {
        this.colocal = colocal;
    }

    public boolean isPreProcess() {
        return preProcess;
    }

    public void setPreProcess(boolean preProcess) {
        this.preProcess = preProcess;
    }

    public double getCurveFitTol() {
        return curveFitTol;
    }

    public void setCurveFitTol(double curveFitTol) {
        this.curveFitTol = curveFitTol;
    }

    public int getnMax() {
        return nMax;
    }

    public void setnMax(int nMax) {
        this.nMax = nMax;
    }

    public boolean isGpu() {
        return gpu;
    }

    public void setGpu(boolean gpu) {
        this.gpu = gpu;
    }

    public double getMinTrajDist() {
        return minTrajDist;
    }

    public void setMinTrajDist(double minTrajDist) {
        this.minTrajDist = minTrajDist;
    }

    public double getTrackLength() {
        return trackLength;
    }

    public void setTrackLength(double trackLength) {
        this.trackLength = trackLength;
    }

    public boolean isUseCals() {
        return useCals;
    }

    public void setUseCals(boolean useCals) {
        this.useCals = useCals;
    }

    public boolean isExtractsigs() {
        return extractsigs;
    }

    public void setExtractsigs(boolean extractsigs) {
        this.extractsigs = extractsigs;
    }

    public double getMsdThresh() {
        return msdThresh;
    }

    public void setMsdThresh(double msdThresh) {
        this.msdThresh = msdThresh;
    }

    public double getColocalThresh() {
        return colocalThresh;
    }

    public void setColocalThresh(double colocalThresh) {
        this.colocalThresh = colocalThresh;
    }

    public double getSigEstGreen() {
        return sigEstGreen;
    }

    public void setSigEstGreen(double sigEstGreen) {
        this.sigEstGreen = sigEstGreen;
    }

    public double getSigEstRed() {
        return sigEstRed;
    }

    public void setSigEstRed(double sigEstRed) {
        this.sigEstRed = sigEstRed;
    }

    public int getMinMSDPoints() {
        return minMSDPoints;
    }

    public void setMinMSDPoints(int minMSDPoints) {
        this.minMSDPoints = minMSDPoints;
    }

    public int getDetectionMode() {
        return detectionMode;
    }

    public void setDetectionMode(int detectionMode) {
        this.detectionMode = detectionMode;
    }

    public boolean isFitC2() {
        return fitC2;
    }

    public void setFitC2(boolean fitC2) {
        this.fitC2 = fitC2;
    }

    public boolean isTrackRegions() {
        return trackRegions;
    }

    public void setTrackRegions(boolean trackRegions) {
        this.trackRegions = trackRegions;
    }

    public double getBlobSize() {
        return blobSize;
    }

    public void setBlobSize(double blobSize) {
        this.blobSize = blobSize;
    }

    public double getFilterRadius() {
        return filterRadius;
    }

    public void setFilterRadius(double filterRadius) {
        this.filterRadius = filterRadius;
    }

    public int getMotionModel() {
        return motionModel;
    }

    public void setMotionModel(int motionModel) {
        this.motionModel = motionModel;
    }

    public int getMaxFrameGap() {
        return maxFrameGap;
    }

    public void setMaxFrameGap(int maxFrameGap) {
        this.maxFrameGap = maxFrameGap;
    }

    public double getBlobThresh() {
        return blobThresh;
    }

    public void setBlobThresh(double blobThresh) {
        this.blobThresh = blobThresh;
    }
}
