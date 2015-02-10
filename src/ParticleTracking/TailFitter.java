/*
 * Copyright (C) 2014 David Barry <david.barry at cancer.org.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package ParticleTracking;

import UtilClasses.Utilities;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ZProjector;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Rectangle;
import java.io.File;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.MathArrays;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class TailFitter extends IsoGaussianFitter {

    private static double spatialRes = 0.133333;
    private static double sigmaEst = 0.145;
    double sqrt2 = Math.pow(2.0, 0.5);

//    public static void main(String args[]) {
////        Random r = new Random();
//        File directory = Utilities.getFolder(null, "Select input folder", true);
//        File files[] = directory.listFiles();
//        int dirSize = files.length;
//        ImagePlus temp = IJ.openImage(files[0].getAbsolutePath());
//        ImageProcessor tempIP = temp.getProcessor();
//        int stackwidth = tempIP.getWidth();
//        int stackheight = tempIP.getHeight();
//        temp.close();
//        System.out.println(directory);
////        for (int i = 0; i < 100; i++) {
////            System.out.print(i);
//        ImageStack stack = new ImageStack(stackwidth, stackheight);
//        for (int j = 0; j < dirSize; j++) {
////                int fileindex = r.nextInt(dirSize);
////                ImagePlus imp = IJ.openImage(files[fileindex].getAbsolutePath());
//            ImagePlus imp = IJ.openImage(files[j].getAbsolutePath());
//            stack.addSlice(imp.getProcessor().duplicate());
//            imp.close();
//        }
//        ZProjector zproj = new ZProjector(new ImagePlus("", stack));
//        zproj.setMethod(ZProjector.AVG_METHOD);
//        zproj.doProjection();
//        ImageProcessor stackAverage = zproj.getProjection().getProcessor();
//        Rectangle cropRoi = new Rectangle(0, 6, stackAverage.getWidth() - 2, stackAverage.getHeight() - 12);
//        stackAverage.setRoi(cropRoi);
//        stackAverage = stackAverage.crop();
//        ImageStatistics stats = stackAverage.getStatistics();
//        double max = stats.max;
//        double min = stats.min;
//        stackAverage.subtract(min);
//        stackAverage.multiply(1.0 / (max - min));
//        int width = stackAverage.getWidth();
//        int height = stackAverage.getHeight();
//        double xVals[] = new double[width];
//        double yVals[] = new double[height];
//        double zVals[][] = new double[width][height];
//        for (int y = 0; y < height; y++) {
//            yVals[y] = y * spatialRes;
//            for (int x = 0; x < width; x++) {
//                xVals[x] = x * spatialRes;
//                zVals[x][y] = stackAverage.getPixelValue(x, y);
//            }
//        }
//        TailFitter tf = new TailFitter(xVals, yVals, zVals);
//        tf.doFit(TailFitter.sigmaEst);
//        tf.printParams();
////        IJ.saveAs((new ImagePlus("", stackAverage)), "TIFF", "C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails\\Average.tif");
////        }
//        tf.printImage();
//        System.exit(0);
//    }

    public TailFitter(double[] xVals, double[] yVals, double[][] zVals) {
        super();
        numParams = 7;
        this.xData = xVals;
        this.yData = yVals;
        this.zData = new double[xData.length * yData.length];
        for (int x = 0; x < xData.length; x++) {
            for (int y = 0; y < yData.length; y++) {
                this.zData[y * xData.length + x] = zVals[x][y];
            }
        }
        if (xData != null && yData != null) {
            numPoints = xVals.length * yVals.length;
            for (int i = xVals.length - 1; i >= 0; i--) {
                xData[i] -= xData[0];
            }
            for (int j = yVals.length - 1; j >= 0; j--) {
                yData[j] -= yData[0];
            }
        } else {
            numPoints = 0;
        }
    }

    boolean initialize(double sigmaEst) {
        if (sigmaEst <= 0.0 || xData == null || yData == null || zData == null) {
            return false;
        }
        // Calculate some things that might be useful for predicting parametres
        numVertices = numParams + 1;      // need 1 more vertice than parametres,
        simp = new double[numVertices][numVertices];
        next = new double[numVertices];
        maxIter = IterFactor * numParams * numParams;  // Where does this estimate come from?
        restarts = defaultRestarts;
        nRestarts = 0;
        simp[0][0] = 1.0; //lambda
        simp[0][1] = 0.14 * xData.length * spatialRes; //mu
        simp[0][2] = 0.3; //sigma
        simp[0][3] = 0.5; //nu
        simp[0][4] = 0.0; //noise
        simp[0][5] = yData.length * spatialRes / 2.0;
        simp[0][6] = sigmaEst;

        return true;
    }

    public double evaluate1D(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double p22 = p[2] * p[2];
        double a = 0.5 * p[0] * (2.0 * p[1] + p[0] * p22 - 2.0 * x);
        double b = (p[1] + p[0] * p22 - x) / (sqrt2 * p[2]);

        return p[3] * p[0] * Math.exp(a) * Erf.erfc(b);
    }

    public double evaluate2D(double[] p, double xVal, double y) {
        if (p == null) {
            return Double.NaN;
        }
        double v = Math.pow((y - p[5]) / (p[6] * sqrt2), 2.0);

        return xVal * Math.exp(-0.5 * v) + p[4];
    }

    protected boolean sumResiduals(double[] x) {
        if (x == null) {
            return false;
        }
        /*
         * x[numParams] = sumResiduals(x, xData, yData, zData); return true;
         */
        double e;
        x[numParams] = 0.0;
        double tail1d[] = buildTail(x);
        for (int i = 0; i < xData.length; i++) {
            for (int j = 0; j < yData.length; j++) {
                e = evaluate2D(x, tail1d[i + xData.length / 2], yData[j]) - zData[j * xData.length + i];
                x[numParams] = x[numParams] + (e * e);
            }
        }
        return true;
    }

    double[] buildTail(double[] p) {
        double[] emg = new double[xData.length];
        Gaussian gauss = new Gaussian((emg.length - 1.0) / 2.0, sigmaEst / spatialRes);
        double[] gaussian = new double[emg.length];
        for (int i = 0; i < emg.length; i++) {
            gaussian[i] = gauss.value(i);
            emg[i] = evaluate1D(p, xData[i]);
        }
        return MathArrays.convolve(emg, gaussian);
    }

    public void printParams() {
        double params[] = getParams();

        for (int i = 0; i < numParams; i++) {
            System.out.print(" p[" + String.valueOf(i) + "]: " + params[i]);
        }
        System.out.println();
    }

    void printImage() {
        FloatProcessor deconvolved = new FloatProcessor(xData.length, yData.length);
        FloatProcessor convolved = new FloatProcessor(xData.length, yData.length);
        double tail1d[] = buildTail(simp[best]);
        for (int y = 0; y < yData.length; y++) {
            for (int x = 0; x < xData.length; x++) {
                double xVal = tail1d[x + xData.length / 2];
                convolved.putPixelValue(x, y, xVal);
                deconvolved.putPixelValue(x, y, evaluate2D(simp[best], xVal, yData[y]));
            }
        }
        IJ.saveAs((new ImagePlus("", convolved)), "text image", "C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails\\Convolved.txt");
        IJ.saveAs((new ImagePlus("", deconvolved)), "text image", "C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails\\Deconvolved.txt");
    }
}
