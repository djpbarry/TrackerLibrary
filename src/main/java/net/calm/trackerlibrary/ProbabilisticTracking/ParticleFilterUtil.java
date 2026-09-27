package net.calm.trackerlibrary.ProbabilisticTracking;

import ij.ImagePlus;
import ij.ImageStack;
import java.util.Vector;

/**
 * Stateless helper operations used by the particle-filter trackers. Extracted
 * from {@link PFTracking3D} so the base class focuses on the filter state
 * machine rather than low-level array/stack utilities.
 *
 * @author Janick Cardinale, ETH Zurich
 */
public final class ParticleFilterUtil {

    private ParticleFilterUtil() {
    }

    /**
     * Copies a <code>Vector&lt;float[]&gt;</code> data structure.
     *
     * @param aOrig the vector to copy.
     * @return the (deep) copy.
     */
    public static Vector<float[]> copyStateVector(Vector<float[]> aOrig) {
        Vector<float[]> vResVector = new Vector<float[]>(aOrig.size());
        for (float[] vA : aOrig) {
            float[] vResA = new float[vA.length];
            for (int vI = 0; vI < vA.length; vI++) {
                vResA[vI] = vA[vI];
            }
            vResVector.add(vResA);
        }
        return vResVector;
    }

    /**
     * Copies a <code>Vector&lt;Vector&lt;float[]&gt;&gt;</code> data structure.
     * Used to copy the particle vector.
     *
     * @param aOrig the vector to copy.
     * @return the (deep) copy.
     */
    public static Vector<Vector<float[]>> copyParticleVector(Vector<Vector<float[]>> aOrig) {
        Vector<Vector<float[]>> vResVector = new Vector<Vector<float[]>>(aOrig.size());
        for (Vector<float[]> vP : aOrig) {
            vResVector.add(copyStateVector(vP));
        }
        return vResVector;
    }

    /**
     * Recursively searches the brightest voxel in the neighbourhood. Might be
     * used for the initialization.
     *
     * @param aStartX start x coordinate.
     * @param aStartY start y coordinate.
     * @param aStartZ start z coordinate.
     * @param aImageStack the image stack.
     * @return an int array with 3 entries: x, y and z coordinate.
     */
    public static int[] searchLocalMaximumIntensityWithSteepestAscent(int aStartX, int aStartY, int aStartZ, ImageStack aImageStack) {
        int[] vRes = new int[]{aStartX, aStartY, aStartZ};
        float vMaxValue = aImageStack.getProcessor(aStartZ).getPixelValue(aStartX, aStartY);
        for (int vZi = -1; vZi < 2; vZi++) {
            if (aStartZ + vZi > 0 && aStartZ + vZi <= aImageStack.getSize()) {
                for (int vXi = -1; vXi < 2; vXi++) {
                    for (int vYi = -1; vYi < 2; vYi++) {
                        if (aImageStack.getProcessor(aStartZ + vZi).getPixelValue(aStartX + vXi, aStartY + vYi) > vMaxValue) {
                            vMaxValue = aImageStack.getProcessor(aStartZ + vZi).getPixelValue(aStartX + vXi, aStartY + vYi);
                            vRes[0] = aStartX + vXi;
                            vRes[1] = aStartY + vYi;
                            vRes[2] = aStartZ + vZi;
                        }
                    }
                }
            }
        }
        if (vMaxValue > aImageStack.getProcessor(aStartZ).getPixelValue(aStartX, aStartY)) {
            return searchLocalMaximumIntensityWithSteepestAscent(vRes[0], vRes[1], vRes[2], aImageStack);
        }
        return vRes;
    }

    /**
     * Adds the intensities of two 2D arrays.
     *
     * @param aResult the first image, accumulated in place.
     * @param aImageToAdd the image added to {@code aResult}.
     */
    public static void addImage(float[][] aResult, float[][] aImageToAdd) {
        int vIMax = Math.min(aResult.length, aImageToAdd.length);
        int vJMax = Math.min(aResult[0].length, aImageToAdd[0].length);
        for (int vI = 0; vI < vIMax; vI++) {
            for (int vJ = 0; vJ < vJMax; vJ++) {
                aResult[vI][vJ] += aImageToAdd[vI][vJ];
            }
        }
    }

    /**
     * Sets all values in the array to {@code aValue}.
     *
     * @param aArray the array to fill.
     * @param aValue the value to assign.
     */
    public static void initArrayToValue(float[][] aArray, float aValue) {
        for (int vI = 0; vI < aArray.length; vI++) {
            for (int vJ = 0; vJ < aArray[0].length; vJ++) {
                aArray[vI][vJ] = aValue;
            }
        }
    }

    /**
     * Returns a copy of a single frame. Note that the properties of the
     * ImagePlus have to be correct.
     *
     * @param aMovie the movie.
     * @param aFrameNumber the (1-based) frame number.
     * @return the frame copy.
     */
    public static ImageStack getAFrameCopy(ImagePlus aMovie, int aFrameNumber) {
        if (aFrameNumber > aMovie.getNFrames() || aFrameNumber < 1) {
            throw new IllegalArgumentException();
        }
        int vS = aMovie.getNSlices();
        return getSubStackFloatCopy(aMovie.getStack(), (aFrameNumber - 1) * vS + 1, aFrameNumber * vS);
    }

    /**
     * Returns a copy of a substack (i.e. frames).
     *
     * @param aImageStack the stack to crop.
     * @param aStartPos 1 &le; aStartPos &le; aImageStack.size().
     * @param aEndPos 1 &le; aStartPos &le; aEndPos &le; aImageStack.size().
     * @return a copy of the substack.
     */
    public static ImageStack getSubStackFloatCopy(ImageStack aImageStack, int aStartPos, int aEndPos) {
        ImageStack res = new ImageStack(aImageStack.getWidth(), aImageStack.getHeight());
        if (!(aStartPos < 1 || aEndPos < 0)) {
            for (int vI = aStartPos; vI <= aEndPos; vI++) {
                res.addSlice(aImageStack.getSliceLabel(vI), aImageStack.getProcessor(vI).convertToFloat().duplicate());
            }
        }
        return res;
    }

    /**
     * Returns a substack (i.e. frames) as float, without duplicating the pixel
     * data.
     *
     * @param aImageStack the stack to crop.
     * @param aStartPos 1 &le; aStartPos &le; aImageStack.size().
     * @param aEndPos 1 &le; aStartPos &le; aEndPos &le; aImageStack.size().
     * @return the substack.
     */
    public static ImageStack getSubStackFloat(ImageStack aImageStack, int aStartPos, int aEndPos) {
        ImageStack res = new ImageStack(aImageStack.getWidth(), aImageStack.getHeight());
        if (!(aStartPos < 1 || aEndPos < 0)) {
            for (int vI = aStartPos; vI <= aEndPos; vI++) {
                res.addSlice(aImageStack.getSliceLabel(vI), aImageStack.getProcessor(vI).convertToFloat());
            }
        }
        return res;
    }
}
