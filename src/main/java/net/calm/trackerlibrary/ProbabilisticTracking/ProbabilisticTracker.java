/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.trackerlibrary.ProbabilisticTracking;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Vector;

import net.calm.iaclasslibrary.IAClasses.Utils;

/**
 *
 * @author barry05
 */
public class ProbabilisticTracker extends PFTracking3D {

    private float mBackground = 1;
    private float[] mSigmaOfDynamics = {20, 20, 20, 1};
    private boolean mDoPrecisionCorrection = true;
    int mDimOfState = 7;

    @Override
    public boolean getMDoPrecisionOptimization() {
        return mDoPrecisionCorrection;
    }

    @Override
    public String[] getMDimensionsDescription() {
        String[] vS = {"x[nm]", "y[nm]", "z[nm]", "Intensity"};
        return vS;
    }

    @Override
    public float[] getMSigmaOfRandomWalk() {
        return new float[]{1f, 1f, 1f, 0, 0, 0, 1f};
    }

    @Override
    protected void drawFromProposalDistribution(float[] particle,
            float pxWidthInNm, float pxDepthInNm) {
        particle[3] = particle[3] + (float) mRandomGenerator.nextGaussian() * (mSigmaOfDynamics[0] / pxWidthInNm);
        particle[4] = particle[4] + (float) mRandomGenerator.nextGaussian() * (mSigmaOfDynamics[1] / pxWidthInNm);
        particle[5] = particle[5] + (float) mRandomGenerator.nextGaussian() * (mSigmaOfDynamics[2] / pxDepthInNm);
        particle[0] = particle[0] + particle[3];
        particle[1] = particle[1] + particle[4];
        particle[2] = particle[2] + particle[5];
        particle[6] = particle[6] + (float) mRandomGenerator.nextGaussian() * mSigmaOfDynamics[3];
        if (particle[6] < mBackground + 1) {
            particle[6] = mBackground + 1;
        }
    }

    @Override
    protected void paintOnCanvas(Graphics aG, double magnification,
            int activeFrame) {
        if (mStateVectorsMemory.elementAt(activeFrame - 1) == null) {
            return;
        }
        for (float[] vState : mStateVectorsMemory.elementAt(activeFrame - 1)) {
            int vX = (int) Math.round(vState[0] * magnification);
            int vY = (int) Math.round(vState[1] * magnification);
            aG.setColor(Color.yellow);
            aG.drawLine(vX - 5, vY, vX + 5, vY);
            aG.drawLine(vX, vY - 5, vX, vY + 5);
        }
    }

    @Override
    protected float[][][] generateIdealImage_3D(int aw, int ah, int as,
            float[] particle, int background, float pxWidthInNm,
            float pxDepthInNm) {
        float[][][] vIdealImage = new float[as][ah][aw];
        addBackgroundToImage(vIdealImage, mBackground);
        addFeaturePointTo3DImage(vIdealImage, new Point3D(particle[0], particle[1], particle[2]), particle[6], aw, ah, as, pxWidthInNm, pxDepthInNm, null);
        return vIdealImage;
    }

    @Override
    protected boolean showParameterDialog() {
        GenericDialog vGenericDialog = new GenericDialog("Enter search radius parameters", IJ.getInstance());
        for (int vD = 0; vD < mSigmaOfDynamics.length; vD++) {
            vGenericDialog.addNumericField(mDimensionsDescription[vD], mSigmaOfDynamics[vD], 2);
        }
        vGenericDialog.showDialog();
        for (int vD = 0; vD < mSigmaOfDynamics.length; vD++) {
            mSigmaOfDynamics[vD] = (float) vGenericDialog.getNextNumber();
        }
        if (vGenericDialog.wasCanceled()) {
            return false;
        }
        return true;
    }

    @Override
    protected void calcFromHereButtonPressed() {
        ImagePlus imp = getMOriginalImagePlus();
        ImageProcessor proc = imp.getImageStack().getProcessor(1);
        ArrayList<int[]> maxima = Utils.findLocalMaxima(3, proc, 200.0, true, true);
        int size = maxima.size();
        for (int i = 0; i < size; i++) {
            int[] thismax = maxima.get(i);
            float[] firstState = new float[]{thismax[0], thismax[1], 0.0f, 0.0f,
                0.0f, 0.0f, proc.getPixelValue(thismax[0], thismax[1])};
            mStateVectors.add(firstState);
        }
        setMStateOfFilter(STATE_OF_FILTER.READY_TO_RUN);
        //relaunch the filter from the current slice on with the current initialization if there is one
        if (mStateOfFilter != STATE_OF_FILTER.VISUALIZING && mStateOfFilter != STATE_OF_FILTER.READY_TO_RUN) {
            IJ.showMessage("No initialization done yet.");
            return;
        }
        if (mStateOfFilter == STATE_OF_FILTER.VISUALIZING) { //restart from corrected frame
            if (mStateVectorsMemory.elementAt(mZProjectedImagePlus.getCurrentSlice() - 1) != null) {
                //first setup the state vector then run
                mStateVectors = copyStateVector(mStateVectorsMemory.elementAt(mZProjectedImagePlus.getCurrentSlice() - 1));
                mFrameOfInitialization = mZProjectedImagePlus.getCurrentSlice();
                ProbabilisticTracker.this.run(new FloatProcessor(1, 1));
                mStateOfFilter = STATE_OF_FILTER.RUNNING;
            } else {
                IJ.showMessage("No initialization in this frame.");
            }
        }
        if (mStateOfFilter == STATE_OF_FILTER.READY_TO_RUN) {
            //we're sure that there is a correct initialization at a certain frame
            mStateOfFilter = STATE_OF_FILTER.RUNNING;
            ProbabilisticTracker.this.run(new FloatProcessor(1, 1));
            mStateOfFilter = STATE_OF_FILTER.VISUALIZING;
        }
    }

    @Override
    protected boolean[][][] generateParticlesIntensityBitmap_3D(Vector< float[]> setOfParticles, int aW, int aH, int aS) {
        boolean[][][] vBitmap = new boolean[aS][aH][aW];
// convert to pixel distance and multiply with 3: 4
        float vMaxDistancexy = 3 * mSigmaPSFxy / getPixelWidthInNm();
        float vMaxDistancez = 3 * mSigmaPSFz / getPixelDepthInNm();
// get a bounding box around the each feature point
        int vXStart, vXEnd, vYStart, vYEnd, vZStart, vZEnd;
        for (float[] vParticle : setOfParticles) {
            if (vParticle[0] - vMaxDistancexy < 0) {
                vXStart = 0;
            } else {
                vXStart = (int) Math.round(vParticle[0] - vMaxDistancexy);
            }
            if (vParticle[0] + vMaxDistancexy >= aW) {
                vXEnd = aW - 1;
            } else {
                vXEnd = (int) Math.round(vParticle[0] + vMaxDistancexy);
            }
            if (vParticle[1] - vMaxDistancexy < 0) {
                vYStart = 0;
            } else {
                vYStart = (int) Math.round(vParticle[1] - vMaxDistancexy);
            }
            if (vParticle[1] + vMaxDistancexy >= aH) {
                vYEnd = aH - 1;
            } else {
                vYEnd = (int) Math.round(vParticle[1] + vMaxDistancexy);
            }
            if (vParticle[2] - vMaxDistancez < 0) {
                vZStart = 0;
            } else {
                vZStart = (int) Math.round(vParticle[2] - vMaxDistancez);
            }
            if (vParticle[2] + vMaxDistancez >= aS) {
                vZEnd = aS - 1;
            } else {
                vZEnd = (int) Math.round(vParticle[2] + vMaxDistancez);
            }
            for (int vZ = vZStart; vZ <= vZEnd && vZ < vBitmap.length; vZ++) {
                for (int vY = vYStart; vY <= vYEnd && vY < vBitmap[vZ].length; vY++) {
                    for (int vX = vXStart; vX <= vXEnd && vX < vBitmap[vZ][vY].length; vX++) {
                        vBitmap[vZ][ vY][ vX] = true;
                    }
                }
            }
        }
        return vBitmap;
    }
}
