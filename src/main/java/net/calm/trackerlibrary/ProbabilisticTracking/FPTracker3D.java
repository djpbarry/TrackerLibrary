package net.calm.trackerlibrary.ProbabilisticTracking;

import ij.IJ;
import ij.ImageStack;
import ij.gui.GenericDialog;
import java.awt.Color;
import java.awt.Graphics;
import java.util.Vector;

/**
 * This is a derived class of the Particle Filter Tracking base class
 * (PFTracking3D). It implements a feature point tracker without any further
 * dynamic model (random walk).
 *
 * @author Janick Cardinale, ETHZ 2008
 * @version 1.0
 * @see PFTracking3D
 *
 * <p><b>Disclaimer</b> <br>IN NO EVENT SHALL THE ETH BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
 * EVEN IF THE ETH HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE ETH
 * SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE ETH HAS NO
 * OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.<p>
 *
 */
public class FPTracker3D extends PFTracking3D {

    private float mBackground = 1;
    private float[] mSigmaOfDynamics = {100, 100, 100, 1};
    private boolean mDoPrecisionCorrection = true;

    @Override
    protected float[][][] generateIdealImage_3D(int aw, int ah, int as,
            float[] particle, int background, float pxWidthInNm,
            float pxDepthInNm) {
        float[][][] vIdealImage = new float[as][ah][aw];
        addBackgroundToImage(vIdealImage, mBackground);
        addFeaturePointTo3DImage(vIdealImage, new Point3D(particle[0], particle[1], particle[2]), particle[3], aw, ah, as, pxWidthInNm, pxDepthInNm, null);
        return vIdealImage;
    }

    @Override
    protected void mouseReleased(int ax, int ay) {
        if (getMStateOfFilter() == PFTracking3D.STATE_OF_FILTER.INIT) {
            setMStateOfFilter(STATE_OF_FILTER.READY_TO_RUN);
            ImageStack vInitFrame = getAFrameCopy(getMOriginalImagePlus(), getMZProjectedImagePlus().getCurrentSlice());
            float[] vZInformation = calculateExpectedZPositionAt(ax, ay, vInitFrame);
            float[] vFirstState = new float[]{ax, ay, vZInformation[0], vZInformation[1]};
            mStateVectors.add(vFirstState);
        }
    }

    /**
     * Calculates the expected mean of a gaussian fitted to a Ray trough the
     * imagestack.
     *
     * @param aX The x position of the ray
     * @param aY The y position of the ray
     * @param aIS The imageStack where the intensities are read out.
     * @return the expected z position of a Gaussian in [1; aIS.getSize<code></code>]
     * at position 0 and the maximal intensity at position (ax,ay) in this
     * stack.
     */
    private float[] calculateExpectedZPositionAt(int aX, int aY, ImageStack aIS) {
        float vMaxInt = 0;
        int vMaxSlice = 0;
        for (int vZ = 0; vZ < mNSlices; vZ++) {
            float vThisInt;
            if ((vThisInt = aIS.getProcessor(vZ + 1).getf(aX, aY)) > vMaxInt) {
                vMaxInt = vThisInt;
                vMaxSlice = vZ;
            }

        }
        float vSumOfIntensities = 0f;
        float vRes = 0f;
        int vStartSlice = Math.max(0, vMaxSlice - 2);
        int vStopSlice = Math.min(mNSlices - 1, vMaxSlice + 2);
        for (int vZ = vStartSlice; vZ <= vStopSlice; vZ++) {
            vSumOfIntensities += aIS.getProcessor(vZ + 1).getf(aX, aY);
            vRes += (vZ + 1) * aIS.getProcessor(vZ + 1).getf(aX, aY);
        }
        return new float[]{vRes / vSumOfIntensities, vMaxInt};
    }

    private void addBackgroundToImage(float[][][] aImage, float aBackground) {
        for (float[][] vSlice : aImage) {
            for (float[] vRow : vSlice) {
                for (int vI = 0; vI < vRow.length; vI++) {
                    vRow[vI] += aBackground;
                }
            }
        }
    }

    @Override
    protected void paintParticleOnCanvas(Graphics ag, float[] particle,
            double magnification) {
        ag.setColor(Color.green);
        ag.drawRect((int) (magnification * particle[0] + .5f), (int) (magnification * particle[1] + .5f),
                1, 1);
    }

    private void addFeaturePointTo3DImage(float[][][] aImage, Point3D aPoint, float aIntensity, int aW, int aH, int aS, float aPxWidthInNm, float aPxDepthInNm, float aGhostImage[][]) {
        float vVarianceXYinPx = mSigmaPSFxy * mSigmaPSFxy / (aPxWidthInNm * aPxWidthInNm);
        float vVarianceZinPx = mSigmaPSFz * mSigmaPSFz / (aPxDepthInNm * aPxDepthInNm);
        float vMaxDistancexy = 3 * mSigmaPSFxy / aPxWidthInNm;
        float vMaxDistancez = 3 * mSigmaPSFz / aPxDepthInNm; //in pixel!

        int vXStart, vXEnd, vYStart, vYEnd, vZStart, vZEnd;//defines a bounding box around the tip
        if (aPoint.mX + .5f - (vMaxDistancexy + .5f) < 0) {
            vXStart = 0;
        } else {
            vXStart = (int) (aPoint.mX + .5f) - (int) (vMaxDistancexy + .5f);
        }
        if (aPoint.mY + .5f - (vMaxDistancexy + .5f) < 0) {
            vYStart = 0;
        } else {
            vYStart = (int) (aPoint.mY + .5f) - (int) (vMaxDistancexy + .5f);
        }
        if (aPoint.mZ + .5f - (vMaxDistancez + .5f) < 0) {
            vZStart = 0;
        } else {
            vZStart = (int) (aPoint.mZ + .5f) - (int) (vMaxDistancez + .5f);
        }
        if (aPoint.mX + .5f + (vMaxDistancexy + .5f) >= aW) {
            vXEnd = aW - 1;
        } else {
            vXEnd = (int) (aPoint.mX + .5f) + (int) (vMaxDistancexy + .5f);
        }
        if (aPoint.mY + .5f + (vMaxDistancexy + .5f) >= aH) {
            vYEnd = aH - 1;
        } else {
            vYEnd = (int) (aPoint.mY + .5f) + (int) (vMaxDistancexy + .5f);
        }
        if (aPoint.mZ + .5f + (vMaxDistancez + .5f) >= aS) {
            vZEnd = aS - 1;
        } else {
            vZEnd = (int) (aPoint.mZ + .5f) + (int) (vMaxDistancez + .5f);
        }

        for (int vZ = vZStart; vZ <= vZEnd; vZ++) {
            for (int vY = vYStart; vY <= vYEnd; vY++) {
                for (int vX = vXStart; vX <= vXEnd; vX++) {
                    aImage[vZ][vY][vX] += (float) (aIntensity
                            * Math.pow(Math.E, -(Math.pow(vX - aPoint.mX + .5f, 2) + Math.pow(vY - aPoint.mY + .5f, 2)) / (2 * vVarianceXYinPx))
                            * Math.pow(Math.E, -Math.pow(vZ - aPoint.mZ + .5f, 2) / (2 * vVarianceZinPx)));
                }
            }
        }
    }

    @Override
    protected void drawFromProposalDistribution(float[] particle,
            float pxWidthInNm, float pxDepthInNm) {
        particle[0] = particle[0] + (float) mRandomGenerator.nextGaussian() * (mSigmaOfDynamics[0] / pxWidthInNm);
        particle[1] = particle[1] + (float) mRandomGenerator.nextGaussian() * (mSigmaOfDynamics[1] / pxWidthInNm);
        particle[2] = particle[2] + (float) mRandomGenerator.nextGaussian() * (mSigmaOfDynamics[2] / pxDepthInNm);
        particle[3] = particle[3] + (float) mRandomGenerator.nextGaussian() * mSigmaOfDynamics[3];
        if (particle[3] < mBackground) {
            particle[3] = mBackground + 1;
        }
    }

    @Override
    protected boolean[][][] generateParticlesIntensityBitmap_3D(
            Vector<float[]> setOfParticles, int aW, int aH, int aS) {
        boolean[][][] vBitmap = new boolean[aS][aH][aW];
        float vMaxDistancexy = 3 * mSigmaPSFxy / getPixelWidthInNm();
        float vMaxDistancez = 3 * mSigmaPSFz / getPixelDepthInNm(); //in pixel!

        int vXStart, vXEnd, vYStart, vYEnd, vZStart, vZEnd;//defines a bounding box around the tip
        for (float[] vParticle : setOfParticles) {
            if (vParticle[0] + .5f - (vMaxDistancexy + .5f) < 0) {
                vXStart = 0;
            } else {
                vXStart = (int) (vParticle[0] + .5f) - (int) (vMaxDistancexy + .5f);
            }
            if (vParticle[1] + .5f - (vMaxDistancexy + .5f) < 0) {
                vYStart = 0;
            } else {
                vYStart = (int) (vParticle[1] + .5f) - (int) (vMaxDistancexy + .5f);
            }
            if (vParticle[2] + .5f - (vMaxDistancez + .5f) < 0) {
                vZStart = 0;
            } else {
                vZStart = (int) (vParticle[2] + .5f) - (int) (vMaxDistancez + .5f);
            }
            if (vParticle[0] + .5f + (vMaxDistancexy + .5f) >= aW) {
                vXEnd = aW - 1;
            } else {
                vXEnd = (int) (vParticle[0] + .5f) + (int) (vMaxDistancexy + .5f);
            }
            if (vParticle[1] + .5f + (vMaxDistancexy + .5f) >= aH) {
                vYEnd = aH - 1;
            } else {
                vYEnd = (int) (vParticle[1] + .5f) + (int) (vMaxDistancexy + .5f);
            }
            if (vParticle[2] + .5f + (vMaxDistancez + .5f) >= aS) {
                vZEnd = aS - 1;
            } else {
                vZEnd = (int) (vParticle[2] + .5f) + (int) (vMaxDistancez + .5f);
            }

            for (int vZ = vZStart; vZ <= vZEnd; vZ++) {
                for (int vY = vYStart; vY <= vYEnd; vY++) {
                    for (int vX = vXStart; vX <= vXEnd; vX++) {
                        vBitmap[vZ][vY][vX] = true;
                    }
                }
            }
        }
        return vBitmap;
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
    public float[] getMSigmaOfRandomWalk() {
        return new float[]{1, 1, 1, 1};
    }

    @Override
    public void setMBackground(float background) {
        mBackground = background;

    }

    @Override
    public String[] getMDimensionsDescription() {
        String[] vS = {"x[nm]", "y[nm]", "z[nm]", "Intensity"};
        return vS;
    }

    @Override
    public boolean getMDoPrecisionOptimization() {
        return mDoPrecisionCorrection;
    }

    /**
     * Shows the dialog to enter the search radii 'sigma' of each paramter in
     * the state vector.
     *
     * @return false if cancelled, else true.
     */
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
}
