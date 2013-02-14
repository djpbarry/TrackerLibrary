/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package VirusTracking;

import IAClasses.IsoGaussian;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.text.DecimalFormat;
import java.util.Random;

/**
 *
 * @author barry05
 */
public class TestGenerator {

    DecimalFormat indFormat = new DecimalFormat("000");
    private double noise = 3.0;
    private Random rand = new Random();
    private double numAp = 1.4;
    private double lambda = 602.0;
    private double res = 0.083;
    private double sigmaEstPix = 0.305 * lambda / (numAp * res * 1000.0);

    public static void main(String args[]) {
        TestGenerator tg = new TestGenerator();
        tg.generateFluorophoreCircle(150, 1000, 1000, 200, 83.0, 0.01,
                "C:/Users/barry05/Desktop/Test_Data_Sets/Test_Generator_Output");
    }

    public TestGenerator() {
    }

    public void generate() {
        DecimalFormat indFormat = new DecimalFormat("000");
        int width = 640, height = 480;
        Co_Localise cl = new Co_Localise();
        double res = Timelapse_Analysis.getSpatialRes();
        for (int i = 0; i < 50; i++) {
            ByteProcessor image = new ByteProcessor(width, height);
            IsoGaussian g1 = new IsoGaussian((i * 0.2 + 320) * res, 240.0 * res, 100.0, 2.0 * res, 2.0 * res, 0.1);
            IsoGaussian g2 = new IsoGaussian((320 - i * 0.2) * res, 240.0 * res, 100.0, 2.0 * res, 2.0 * res, 0.1);
            cl.draw2DGaussian(image, g1, 0.0, Timelapse_Analysis.spatialRes);
            cl.draw2DGaussian(image, g2, 0.0, Timelapse_Analysis.spatialRes);
            IJ.saveAs(new ImagePlus("", image.duplicate()), "PNG",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\"
                    + indFormat.format(i));
        }
    }

    public void generateMulti(int n, int width, int height, int length) {
        int totalcount = n;
        Co_Localise cl = new Co_Localise();
        MotileGaussian particles[] = new MotileGaussian[n];
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            particles[i] = new MotileGaussian(width * res * r.nextDouble(), height * res * r.nextDouble(),
                    255.0, sigmaEstPix, sigmaEstPix, 0.1, 0.02, true, false);
        }
        for (int i = 0; i < length; i++) {
            ByteProcessor image = new ByteProcessor(width, height);
            image.setColor(255);
            for (int j = 0; j < n; j++) {
                if (particles[j] != null) {
                    cl.draw2DGaussian(image, particles[j], 0.0, res);
                    particles[j].updateVelocity();
                    particles[j].updatePosition();
                    double x = particles[j].getX() / res;
                    double y = particles[j].getY() / res;
                    if (x < -2.0 * particles[j].getXSigma()
                            || x > width + 2.0 * particles[j].getXSigma()
                            || y < -2.0 * particles[j].getYSigma()
                            || y > height + 2.0 * particles[j].getYSigma()) {
//                        particles[j] = null;
                        particles[j] = new MotileGaussian(width * res * r.nextDouble(),
                                height * res * r.nextDouble(), 255.0, sigmaEstPix, sigmaEstPix,
                                0.1, 0.02, true, false);
                        totalcount++;

                    }
                }
            }
            IJ.saveAs(new ImagePlus("", image.duplicate()), "PNG",
                    "C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences\\Simulation\\"
                    + indFormat.format(i));
            System.out.println("Frame:\t" + i + "\tTotal Count:\t" + totalcount);
        }
    }

    public void generateFluorophoreCircle(int radius, int width, int height, int length,
            double finalRes, double thresh, String outputDir) {
        int circum = (int)Math.ceil(2.0 * Math.PI * radius);
        Fluorophore dots[] = new Fluorophore[circum];
        double sigma = 0.305f * 602.0 / 1.4;
        double theta;
        int cx = width / 2;
        int cy = height / 2;
        for (int i = 0; i < circum; i++) {
            theta = i * 2.0 * Math.PI / circum;
            double x = cx + radius * Math.cos(theta);
            double y = cy + radius * Math.sin(theta);
            dots[i] = new Fluorophore(x, y, 255.0, thresh);
//            System.out.println("x: " + x + " y: " + y + " theta: " + theta);
        }
        runGenerator(length, width, height, dots, sigma, finalRes, outputDir);
    }

//    public void generateFilledFluorophoreCircle(int n, int width, int height,
//            int length, double finalRes) {
//        DecayingFluorophore dots[] = new DecayingFluorophore[n];
//        Random r = new Random();
//        double radius = width
//                / 15.0;
//        double res = 125.0 / radius;
//        double sigma = (0.305f * 602.0 / 1.4)
//                / res;
//        double maxNoise = width / 600.0;
//        int cx = width / 2;
//        int cy =
//                height / 2;
//        GaussianBlur gb = new GaussianBlur();
//        for (int i = 0; i < n;
//                i++) {
//            double theta = 2.0 * Math.PI * r.nextDouble();
//            double nradius =
//                    radius * (1.0 - Math.pow(r.nextDouble(), 2.0));
//            double x = cx + nradius
//                    * Math.cos(theta);
//            double y = cy + nradius * Math.sin(theta);
//            dots[i] = new DecayingFluorophore(x, y, 255.0, 0.05);
//        }
//        for (int i = 0; i < length;
//                i++) {
//            FloatProcessor image = new FloatProcessor(width, height);
//            image.setValue(0.0);
//            image.fill();
//            for (int j = 0; j < n; j++) {
//                dots[j].updateMag();
//                int x = (int) Math.round(dots[j].getX());
//                int y =
//                        (int) Math.round(dots[j].getY());
//                double mag = dots[j].getCurrentMag()
//                        + image.getPixelValue(x, y);
//                image.putPixelValue(x, y, mag);
//            }
//            IJ.saveAs(new ImagePlus("", image), "TIF", "C:\\Users\\barry05\\Desktop\\Test Data Sets\\Tracking Test      Sequences\\Simulation\\Original_" + indFormat.format(i));
//            gb.blurGaussian(image, sigma, sigma, 0.001);
//            image.setInterpolationMethod(ImageProcessor.BICUBIC);
//            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width * res / finalRes))),
//                    "TIF", "C:\\Users\\barry05\\Desktop\\Test Data Sets\\Tracking Test      Sequences\\Simulation\\BlurredAndScaled_" + indFormat.format(i));
//        }
//    }
    public void generateFilledFluorophoreCircle(int n, int width, int height, int length, double finalRes) {
        DecayingFluorophore dots[] = new DecayingFluorophore[n];
        Random r = new Random();
        double radius = 150.0;
        double radius2 = radius * radius;
        double res = 1.0;
        double sigma = (0.305f * 602.0 / 1.4) / res;
        double scope = 2 * radius;
        int xc = width / 2;
        int yc = height / 2;
        GaussianBlur gb = new GaussianBlur();
        for (int i = 0; i < n;) {
            double x = xc - radius + r.nextDouble() * scope;
            double y = yc - radius + r.nextDouble() * scope;
            if (Math.pow(x - xc, 2.0) + Math.pow(y - yc, 2.0) <= radius2) {
                dots[i] = new DecayingFluorophore(x, y, 255.0, 0.05);
                i++;
            }
        }
        for (int i = 0; i < length; i++) {
            FloatProcessor image = new FloatProcessor(width, height);
            image.setValue(0.0);
            image.fill();
            for (int j = 0; j < n; j++) {
                dots[j].updateMag();
                int x = (int) Math.round(dots[j].getX());
                int y = (int) Math.round(dots[j].getY());
                double mag = dots[j].getCurrentMag() + image.getPixelValue(x, y);
                image.putPixelValue(x, y, mag);
            }
            IJ.saveAs(new ImagePlus("", image), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Test Data Sets\\Tracking Test Sequences\\Simulation\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, sigma, sigma, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width * res / finalRes))), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Test Data Sets\\Tracking Test Sequences\\Simulation\\BlurredAndScaled_"
                    + indFormat.format(i));
        }
    }

    public void generateFilledFluorophoreSquare(int dw, int dh, int width, int height,
            int length, double finalRes, String outputDir, double thresh) {
        Fluorophore dots[] = new Fluorophore[dw * dh];
        double sigma = 0.305f * lambda / numAp;
        int x0 = (width - dw) / 2;
        int y0 = (height - dh) / 2;
        for (int i = x0; i < x0 + dw; i++) {
            for (int j = y0; j < y0 + dh; j++) {
                dots[(j - y0) * dw + (i - x0)] = new Fluorophore(i, j, 255.0, thresh);
            }
        }
        runGenerator(length, width, height, dots, sigma, finalRes, outputDir);
    }

    public void generateFilament(int n, int m, int width, int height, int length, double finalRes) {
        DecayingFluorophore dots[] = new DecayingFluorophore[n];
        Random r = new Random();
        double radius = 250.0;
        double res = 1.0;
        double sigma = (0.305f * 602.0 / 1.4) / res;
        GaussianBlur gb = new GaussianBlur();
        double yc = height / 2.0 - radius;
        double yincs[] = new double[m];
        for (int k = 0; k < m; k++) {
            yincs[k] = k * 2.0 * radius / (m - 1);
        }
        for (int i = 0; i < n; i++) {
            double x = r.nextDouble() * width + noise * rand.nextGaussian();
            double y = yc + yincs[r.nextInt(m)] + noise * rand.nextGaussian();
            dots[i] = new DecayingFluorophore(x, y, 0.0, 0.0);
        }
        for (int i = 0; i < length; i++) {
            FloatProcessor image = new FloatProcessor(width, height);
            image.setValue(0.0);
            image.fill();
            for (int j = 0; j < n; j++) {
                dots[j].updateXMag(i);
                int x = (int) Math.round(dots[j].getX());
                int y = (int) Math.round(dots[j].getY());
                double mag = dots[j].getCurrentMag() + image.getPixelValue(x, y);
                image.putPixelValue(x, y, mag);
            }
            IJ.saveAs(new ImagePlus("", image), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Test Data Sets\\Tracking Test Sequences\\Simulation\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, sigma, sigma, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width * res / finalRes))), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Test Data Sets\\Tracking Test Sequences\\Simulation\\BlurredAndScaled_"
                    + indFormat.format(i));
        }
    }

    public void generateDeterministicGrid(int width, int height, int length, int sep, int border) {
        int m = (int) Math.round((width - 2.0 * border) / sep) + 1;
        int n = (int) Math.round((height - 2.0 * border) / sep) + 1;
        int size = m * n;
        DecayingFluorophore dots[] = new DecayingFluorophore[size];
        GaussianBlur gb = new GaussianBlur();
        double initialSig = 0.305 * lambda / numAp;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                dots[i + j * n] = new DecayingFluorophore(border + i * sep, border + j * sep, 255.0, 0.05);
            }
        }
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < size; j++) {
                FloatProcessor image = new FloatProcessor(width, height);
                image.setValue(0.0);
                image.fill();
                dots[j].updateMag(dots[j].getInitialMag() * (1.0 - (double) i / length));
                for (int k = 0; k < size; k++) {
                    int x = (int) Math.round(dots[k].getX());
                    int y = (int) Math.round(dots[k].getY());
                    double mag = dots[k].getCurrentMag() + image.getPixelValue(x, y);
                    image.putPixelValue(x, y, mag);
//                    System.out.println("x: " + x + " y: " + y + " mag: " + mag);
                }
                ImageProcessor imageCopy = image.duplicate();
                IJ.saveAs(new ImagePlus("", imageCopy), "TIF",
                        "C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences\\Simulation\\Original_"
                        + indFormat.format(i * size + j));
                gb.blurGaussian(imageCopy, initialSig, initialSig, 0.001);
                imageCopy.setInterpolationMethod(ImageProcessor.BICUBIC);
                IJ.saveAs(new ImagePlus("", imageCopy.resize((int) Math.round(width / (res * 1000.0)))), "TIF",
                        "C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences\\Simulation\\BlurredAndScaled_"
                        + indFormat.format(i * size + j));
            }
        }
    }

    public void generateStochasticGrid(int width, int height, int length, int sep, int border) {
        int m = (int) Math.round((width - 2.0 * border) / sep) + 1;
        int n = (int) Math.round((height - 2.0 * border) / sep) + 1;
        int size = m * n;
        DecayingFluorophore dots[] = new DecayingFluorophore[size];
        GaussianBlur gb = new GaussianBlur();
        double initialSig = 0.305 * lambda / numAp;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                dots[i + j * n] = new DecayingFluorophore(border + i * sep, border + j * sep, 255.0, 0.05);
            }
        }
        for (int i = 0; i < length; i++) {
            FloatProcessor image = new FloatProcessor(width, height);
            image.setValue(0.0);
            image.fill();
            for (int j = 0; j < size; j++) {
                dots[j].updateMag();
                int x = (int) Math.round(dots[j].getX());
                int y = (int) Math.round(dots[j].getY());
                double mag = dots[j].getCurrentMag() + image.getPixelValue(x, y);
                image.putPixelValue(x, y, mag);
//                    System.out.println("x: " + x + " y: " + y + " mag: " + mag);
            }
            IJ.saveAs(new ImagePlus("", image), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences\\Simulation\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, initialSig, initialSig, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width / (res * 1000.0)))), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Tracking_Test_Sequences\\Simulation\\BlurredAndScaled_"
                    + indFormat.format(i));
        }
    }

    void runGenerator(int length, int width, int height, Fluorophore dots[],
            double sigma, double finalRes, String outputDir) {
        for (int i = 0; i < length; i++) {
            FloatProcessor image = new FloatProcessor(width, height);
            image.setValue(0.0);
            image.fill();
            int n = dots.length;
            GaussianBlur gb = new GaussianBlur();
            for (int j = 0; j < n; j++) {
                int x = (int) Math.round(dots[j].getX());
                int y = (int) Math.round(dots[j].getY());
                double mag = dots[j].getCurrentMag() + image.getPixelValue(x, y);
                image.putPixelValue(x, y, mag);
                dots[j].updateMag();
            }
            IJ.saveAs(new ImagePlus("", image), "TIF", outputDir + "\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, sigma, sigma, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width / finalRes))),
                    "TIF", outputDir + "\\BlurredAndScaled_" + indFormat.format(i));
        }
    }
}
