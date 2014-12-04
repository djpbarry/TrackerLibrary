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
import java.util.Random;
import org.apache.commons.math3.special.Erf;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class TailFitter extends IsoGaussianFitter {

    private static double spatialRes = 0.133333;
    private static double sigmaEst = 1.06;
    double sqrt2 = Math.pow(2.0, 0.5);

    public static void main(String args[]) {
//        Random r = new Random();
        File directory = Utilities.getFolder(null, "Select input folder", true);
        File files[] = directory.listFiles();
        int dirSize = files.length;
        ImagePlus temp = IJ.openImage(files[0].getAbsolutePath());
        ImageProcessor tempIP = temp.getProcessor();
        int stackwidth = tempIP.getWidth();
        int stackheight = tempIP.getHeight();
        temp.close();
        System.out.println(directory);
//        for (int i = 0; i < 100; i++) {
//            System.out.print(i);
            ImageStack stack = new ImageStack(stackwidth, stackheight);
            for (int j = 0; j < dirSize; j++) {
//                int fileindex = r.nextInt(dirSize);
//                ImagePlus imp = IJ.openImage(files[fileindex].getAbsolutePath());
            ImagePlus imp = IJ.openImage(files[j].getAbsolutePath());
                stack.addSlice(imp.getProcessor().duplicate());
                imp.close();
            }
            ZProjector zproj = new ZProjector(new ImagePlus("", stack));
            zproj.setMethod(ZProjector.AVG_METHOD);
            zproj.doProjection();
            ImageProcessor stackAverage = zproj.getProjection().getProcessor();
            ImageStatistics stats = stackAverage.getStatistics();
            double max = stats.max;
            stackAverage.multiply(1.0 / max);
            Rectangle cropRoi = new Rectangle(0, 6, stackAverage.getWidth() - 2, stackAverage.getHeight() - 12);
            stackAverage.setRoi(cropRoi);
            stackAverage = stackAverage.crop();
            int width = stackAverage.getWidth();
            int height = stackAverage.getHeight();
            double xVals[] = new double[width];
            double yVals[] = new double[height];
            double zVals[][] = new double[width][height];
            for (int y = 0; y < height; y++) {
                yVals[y] = y * spatialRes;
                for (int x = 0; x < width; x++) {
                    xVals[x] = x * spatialRes;
                    zVals[x][y] = stackAverage.getPixelValue(x, y);
                }
            }
            TailFitter tf = new TailFitter(xVals, yVals, zVals);
            tf.doFit(TailFitter.sigmaEst);
            tf.findPeak(Timelapse_Analysis.TRACK_OFFSET, 2.0 - 6.0 * spatialRes);
//            IJ.saveAs((new ImagePlus("", stackAverage)), "TIFF", "C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails\\Averages\\" + i);
//        }
        tf.printImage();
        System.exit(0);
    }

    public TailFitter(double[] xVals, double[] yVals, double[][] zVals) {
        super();
        numParams = 8;
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
        simp[0][0] = 0.15 * xData.length * spatialRes;
        simp[0][1] = 0.25;
        simp[0][2] = 0.4 * xData.length * spatialRes;
        simp[0][3] = 1.1;
        simp[0][4] = 0.5 * yData.length * spatialRes;
        simp[0][5] = 0.2;
        simp[0][6] = 0.4;
        simp[0][7] = 0.35;

        return true;
    }

    public double evaluate(double[] p, double x, double y) {
        if (p == null) {
            return Double.NaN;
        }
        double z = (x - p[0]) / (p[1] * sqrt2);
        double w = (x - p[2]) / (p[3] * sqrt2);
        double v = Math.pow((y - p[4]) / (p[5] * sqrt2), 2.0);

        return p[6] + p[7] * (Erf.erf(z) - Erf.erf(w)) * Math.exp(-0.5 * v);
    }

    public void findPeak(double xOffset, double yOffset) {
        double params[] = getParams();

        double x1 = params[0];
        double x12 = x1 * x1;
        double s1 = params[1];
        double s12 = s1 * s1;
        double x2 = params[2];
        double x22 = x2 * x2;
        double s2 = params[3];
        double s22 = s2 * s2;
        double yCentre = params[4] - yOffset;

        double c = 2.0 * s12 * s22 * (Math.log(s2) - Math.log(s1));
        double d = c + s12 * x22 - s22 * x12;
        double b = 2.0 * (s22 * x1 - s12 * x2);
        double a = s12 - s22;

        double root1 = (-b + Math.sqrt(b * b - 4.0 * a * d)) / (2.0 * a);
        double root2 = (-b - Math.sqrt(b * b - 4.0 * a * d)) / (2.0 * a);
        double xCentre;
        if (root1 < xOffset) {
            xCentre = root2 - xOffset;
        } else {
            xCentre = root1 - xOffset;
        }

        System.out.print(" p[0]: " + params[0]
                + " p[1]: " + params[1]
                + " p[2]: " + params[2]
                + " p[3]: " + params[3]
                + " p[4]: " + params[4]
                + " p[5]: " + params[5]
                + " p[6]: " + params[6]
                + " p[7]: " + params[7]);
        System.out.print(" x0: " + xCentre + " y0: " + yCentre + " R^2: " + params[numParams]);
        System.out.println();
    }

    void printImage() {
        FloatProcessor output = new FloatProcessor(xData.length, yData.length);
        for (int y = 0; y < yData.length; y++) {
            for (int x = 0; x < xData.length; x++) {
                output.putPixelValue(x, y, evaluate(simp[best], xData[x], yData[y]));
            }
        }
        IJ.saveAs((new ImagePlus("", output)), "text image", "C:\\Users\\barry05\\Desktop\\Average.txt");
    }
}
