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
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.MathArrays;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class TailFitter extends IsoGaussianFitter {

    private static double spatialRes = 0.133333;
    private static double sigmaEst = 0.158;
    double sqrt2 = Math.pow(2.0, 0.5);

//    public static void main(String args[]) {
//        File parentDir = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails"), "Select input folder", true);
//        File subDirs[] = parentDir.listFiles();
//        ImagePlus temp = IJ.openImage(subDirs[0].listFiles()[0].getAbsolutePath());
//        ImageProcessor tempIP = temp.getProcessor();
//        int stackwidth = tempIP.getWidth();
//        int stackheight = tempIP.getHeight();
//        temp.close();
//        System.out.println(parentDir);
//        for (int i = 0; i < 100; i++) {
//            System.out.print(i + ",");
//            TailFitter tf = new TailFitter();
//            ImageProcessor stackAverage = tf.buildStackAverage(stackwidth, stackheight, subDirs);
//            Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 2, 1);
//            stackAverage.setRoi(cropRoi);
//            stackAverage = stackAverage.crop();
//            ImageStatistics stats = stackAverage.getStatistics();
//            double max = stats.max;
//            double min = stats.min;
//            stackAverage.subtract(min);
//            stackAverage.multiply(1.0 / (max - min));
//            int width = stackAverage.getWidth();
//            int height = stackAverage.getHeight();
//            double xVals[] = new double[width];
//            double yVals[] = new double[height];
//            double zVals[][] = new double[width][height];
//            for (int y = 0; y < height; y++) {
//                yVals[y] = y * spatialRes;
//                for (int x = 0; x < width; x++) {
//                    xVals[x] = x * spatialRes;
//                    zVals[x][y] = stackAverage.getPixelValue(x, y);
//                }
//            }
//            tf.loadData(xVals, yVals, zVals);
//            tf.doFit(TailFitter.sigmaEst);
//            tf.printParams();
//        }
//        System.exit(0);
//    }
    ImageProcessor buildStackAverage(int stackwidth, int stackheight, File[] subDirs) {
        ImageStack overallStack = new ImageStack(stackwidth, stackheight);
        int nDirs = subDirs.length;
        Random r = new Random();
        for (int j = 0; j < nDirs; j++) {
            File files[] = subDirs[j].listFiles();
            ArrayList<ArrayList> sortedFiles = sortFiles(files);
            int n = sortedFiles.size();
            ImageStack cellStack = new ImageStack(stackwidth, stackheight);
            for (int i = 0; i < n; i++) {
                ArrayList theseFiles = sortedFiles.get(i);
                int m = theseFiles.size();
                ImageStack tailStack = new ImageStack(stackwidth, stackheight);
                for (int l = 0; l < m; l++) {
                    int fileindex = r.nextInt(m);
                    ImagePlus imp = IJ.openImage(files[fileindex].getAbsolutePath());
                    tailStack.addSlice(imp.getProcessor());
                    imp.close();
                }
                if (tailStack.getSize() > 0) {
                    cellStack.addSlice(projectStack(tailStack));
                }
            }
            if (cellStack.getSize() > 0) {
                overallStack.addSlice(projectStack(cellStack));
            }
        }
        if (overallStack.getSize() > 0) {
            return projectStack(overallStack);
        } else {
            return null;
        }
    }

    ImageProcessor projectStack(ImageStack stack) {
        ZProjector zproj = new ZProjector(new ImagePlus("", stack));
        zproj.setMethod(ZProjector.AVG_METHOD);
        zproj.doProjection();
        return zproj.getProjection().getProcessor();
    }

    ArrayList<ArrayList> sortFiles(File[] unsorted) {
        ArrayList<ArrayList> files = new ArrayList();
        int n = unsorted.length;
        for (int i = 0; i < n; i++) {
            String name = unsorted[i].getName();
            Scanner scanner = new Scanner(name).useDelimiter("-");
            String channel = scanner.next();
            int index = Integer.parseInt(scanner.next());
            while (files.size() < index) {
                files.add(new ArrayList<File>());
            }
            files.get(index - 1).add(unsorted[i]);
        }
        return files;
    }

//    public static void main(String args[]) {
//        File directory = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\SuperRes Actin Tails"), "Select input folder", true);
//        File files[] = directory.listFiles();
//        int dirSize = files.length;
//        ImagePlus temp = IJ.openImage(files[0].getAbsolutePath());
//        ImageProcessor tempIP = temp.getProcessor();
//        int stackwidth = tempIP.getWidth();
//        int stackheight = tempIP.getHeight();
//        temp.close();
//        System.out.println(directory);
//        ImageStack stack = new ImageStack(stackwidth, stackheight);
//        for (int j = 0; j < dirSize; j++) {
//            ImagePlus imp = IJ.openImage(files[j].getAbsolutePath());
//            ImageProcessor slice = imp.getProcessor().duplicate();
//            slice.setInterpolationMethod(ImageProcessor.BICUBIC);
//            stack.addSlice(slice.resize(stackwidth, stackheight));
//            imp.close();
//        }
//        ZProjector zproj = new ZProjector(new ImagePlus("", stack));
//        zproj.setMethod(ZProjector.AVG_METHOD);
//        zproj.doProjection();
//        ImageProcessor stackAverage = zproj.getProjection().getProcessor();
//        IJ.saveAs((new ImagePlus("", stackAverage)), "text image", directory.getParent() + "/" + directory.getName() + "_Average.txt");
//        zproj.setMethod(ZProjector.SD_METHOD);
//        zproj.doProjection();
//        IJ.saveAs((new ImagePlus("", zproj.getProjection().getProcessor())), "text image", directory.getParent() + "/" + directory.getName() + "_SD.txt");
//        Rectangle cropRoi = new Rectangle(0, 15, stackAverage.getWidth() - 2, 1);
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
//        TailFitter tf = new TailFitter();
//        tf.loadData(xVals, yVals, zVals);
//        tf.doFit(TailFitter.sigmaEst);
//        tf.printParams();
//        tf.printImage(directory);
//        IJ.saveAs((new ImagePlus("", stackAverage)), "text image", directory.getParent() + "/" + directory.getName() + "_NormAverage.txt");
//        System.exit(0);
//    }

    public TailFitter() {
        super();
        numParams = 5;
    }

    public void loadData(double[] xVals, double[] yVals, double[][] zVals) {
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
        Random r = new Random();
        double noise = 0.1;
        simp[0][0] = 1.0 + r.nextDouble() * noise; //lambda
        simp[0][1] = 0.14 * xData.length * spatialRes + r.nextDouble() * noise; //mu
        simp[0][2] = 0.3 + r.nextDouble() * noise; //sigma
        simp[0][3] = 0.5 + r.nextDouble() * noise; //nu
        simp[0][4] = 0.0 + r.nextDouble() * noise; //noise
//        simp[0][5] = yData.length * spatialRes / 2.0;
//        simp[0][6] = sigmaEst;
//        simp[0][0] = 1.0; //A
//        simp[0][1] = 0.95; //mu
//        simp[0][2] = 0.1; //sigma
//        simp[0][3] = 0.0; //noise
//        simp[0][4] = 0.05;//nu

        return true;
    }

    public double evaluate1DEMG(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double p22 = p[2] * p[2];
        double a = 0.5 * p[0] * (2.0 * p[1] + p[0] * p22 - 2.0 * x);
        double b = (p[1] + p[0] * p22 - x) / (sqrt2 * p[2]);

        return p[3] * p[0] * Math.exp(a) * Erf.erfc(b) + p[4];
    }

    public double evaluate1DEMGFirstDerivative(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double p22 = p[2] * p[2];
        double a = (p[0] / 2.0) * (2.0 * p[1] + p[0] * p22 - 2.0 * x);
        double b = (p[1] + p[0] * p22 - x) / (sqrt2 * p[2]);
        double b2 = b * b;
        double c = p[0] * p[3] * sqrt2 / (p[2] * Math.sqrt(Math.PI));
        double d = p[0] * p[0] * p[3];

        return Math.exp(a) * (c * Math.exp(-b2) - d * Erf.erfc(b));
    }

    double findEMGRoot(int Nmax, double a, double b, double tol, double[] p) {
        int N = 1;
        while (N < Nmax) {
            double fa = evaluate1DEMGFirstDerivative(p, a);
//            double fb = evaluate1DEMGFirstDerivative(p, b);
            double c = (a + b) / 2;
            double fc = evaluate1DEMGFirstDerivative(p, c);
            if (fc == 0.0 || (b - a) / 2 < tol) {
                return c;
            }
            N++;
            if (fa * fc > 0.0) {
                a = c;
            } else {
                b = c;
            }
        }
        return Double.NaN;
    }

    public double evaluate1DGaussianPlusEMG(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }
        double lambda = 0.6648;
        double mu = 0.2377 + .75;
        double sigma = 0.2722;
        double a = 0.5 * lambda * (2.0 * mu + lambda * sigma * sigma - 2.0 * x);
        double b = (mu + lambda * sigma * sigma - x) / (Math.sqrt(2.0) * sigma);

        return p[0] * Math.exp(-0.5 * (Math.pow((x - p[1]) / p[2], 2.0))) + p[4] * lambda * Math.exp(a) * Erf.erfc(b) + p[3];
    }

    public double evaluate1DGaussian(double[] p, double x) {
        if (p == null) {
            return Double.NaN;
        }

        return p[0] * Math.exp(-0.5 * (Math.pow((x - p[1]) / p[2], 2.0))) + p[3];
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
                e = tail1d[i + xData.length / 2] - zData[j * xData.length + i];
//                e = tail1d[i] - zData[j * xData.length + i];
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
            emg[i] = evaluate1DEMG(p, xData[i]);
        }
        return MathArrays.convolve(emg, gaussian);
//        return emg;
    }

    public void printParams() {
        double params[] = getParams();

        for (int i = 0; i < numParams; i++) {
            System.out.print("p[" + String.valueOf(i) + "]:," + params[i] + ",");
        }
        System.out.print("Peak:,x=," + findEMGRoot(10000, params[1], params[1] + 2.0 * params[2], 1.0E-10, params) + ",");
        System.out.println();
    }

    void printImage(File directory) {
        FloatProcessor deconvolved = new FloatProcessor(xData.length, yData.length);
        FloatProcessor convolved = new FloatProcessor(xData.length, yData.length);
        FloatProcessor derivative = new FloatProcessor(xData.length, yData.length);
        double tail1d[] = buildTail(simp[best]);
        for (int y = 0; y < yData.length; y++) {
            for (int x = 0; x < xData.length; x++) {
                convolved.putPixelValue(x, y, tail1d[x + xData.length / 2]);
                deconvolved.putPixelValue(x, y, evaluate1DEMG(simp[best], xData[x]));
                derivative.putPixelValue(x, y, evaluate1DEMGFirstDerivative(simp[best], xData[x]));
            }
        }
        IJ.saveAs((new ImagePlus("", convolved)), "text image", directory.getParent() + "/" + directory.getName() + "_Convolved.txt");
        IJ.saveAs((new ImagePlus("", deconvolved)), "text image", directory.getParent() + "/" + directory.getName() + "_Deconvolved.txt");
        IJ.saveAs((new ImagePlus("", derivative)), "text image", directory.getParent() + "/" + directory.getName() + "_Derivative.txt");
    }
}
