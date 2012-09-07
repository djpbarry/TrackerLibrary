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

//    public static void main(String args[]) {
//        TestGenerator tg = new TestGenerator();
//        tg.generateFilledFluorophoreSquare(500, 2000, 2000, 150, 40.0);
//    }

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
            cl.draw2DGaussian(image, g1, 0.0);
            cl.draw2DGaussian(image, g2, 0.0);
            IJ.saveAs(new ImagePlus("", image.duplicate()), "PNG",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\"
                    + indFormat.format(i));
        }
    }

    public void generateMulti(int n, int width, int height, int length) {
        int totalcount = n;
        Co_Localise cl = new Co_Localise();
        double res = Timelapse_Analysis.getSpatialRes();
        MotileGaussian particles[] = new MotileGaussian[n];
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            particles[i] = new MotileGaussian(width * res * r.nextDouble(), height * res * r.nextDouble(),
                    255.0, 2.0, 2.0, 0.1, 0.02, false, true);
        }
        for (int i = 0; i < length; i++) {
            ByteProcessor image = new ByteProcessor(width, height);
            image.setColor(255);
            for (int j = 0; j < n; j++) {
                if (particles[j] != null) {
                    cl.draw2DGaussian(image, particles[j], 0.0);
                    particles[j].updateVelocity();
                    particles[j].updatePosition();
                    double x = particles[j].getX() / res;
                    double y = particles[j].getY() / res;
                    if (x < -2.0 * particles[j].getXSigma()
                            || x > width + 2.0 * particles[j].getXSigma()
                            || y < -2.0 * particles[j].getYSigma()
                            || y > height + 2.0 * particles[j].getYSigma()) {
                        particles[j] = null;
                        /*
                         * particles[j] = new MotileGaussian(width * res *
                         * r.nextDouble(), height * res * r.nextDouble(), 255.0,
                         * 2.0, 2.0, 0.1, 0.02, false, true); totalcount++;
                         */
                    }
                }
            }
            IJ.saveAs(new ImagePlus("", image.duplicate()), "PNG",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\"
                    + indFormat.format(i));
            System.out.println("Frame:\t" + i + "\tTotal Count:\t" + totalcount);
        }
    }

    public void generateFluorophoreCircle(int n, int width, int height, int length, double finalRes) {
        DecayingFluorophore dots[] = new DecayingFluorophore[n];
        Random r = new Random();
        double radius = width / 15.0;
        double res = 125.0 / radius;
        double sigma = (0.305f * 602.0 / 1.4) / res;
        double maxNoise = width / 60.0;
        double theta;
        int cx = width / 2;
        int cy = height / 2;
        GaussianBlur gb = new GaussianBlur();
        for (int i = 0; i < n; i++) {
            theta = i * 2.0 * Math.PI / n;
            double x = cx + radius * Math.cos(theta) + r.nextGaussian() * maxNoise;
            double y = cy + radius * Math.sin(theta) + r.nextGaussian() * maxNoise;
            dots[i] = new DecayingFluorophore(x, y, 255.0, 0.05);
//            System.out.println("x: " + x + " y: " + y + " theta: " + theta);
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
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, sigma, sigma, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width * res / finalRes))), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\BlurredAndScaled_"
                    + indFormat.format(i));
        }
    }

    public void generateFilledFluorophoreCircle(int n, int width, int height, int length, double finalRes) {
        DecayingFluorophore dots[] = new DecayingFluorophore[n];
        Random r = new Random();
        double radius = width / 15.0;
        double res = 125.0 / radius;
        double sigma = (0.305f * 602.0 / 1.4) / res;
        double maxNoise = width / 60.0;
        int cx = width / 2;
        int cy = height / 2;
        GaussianBlur gb = new GaussianBlur();
        for (int i = 0; i < n; i++) {
            double x = cx - radius + 2.0 * radius * r.nextDouble() + r.nextGaussian() * maxNoise;
            double y = cy - radius + 2.0 * radius * r.nextDouble() + r.nextGaussian() * maxNoise;
            dots[i] = new DecayingFluorophore(x, y, 255.0, 0.05);
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
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, sigma, sigma, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width * res / finalRes))), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\BlurredAndScaled_"
                    + indFormat.format(i));
        }
    }

    public void generateFilledFluorophoreSquare(int n, int width, int height, int length, double finalRes) {
        DecayingFluorophore dots[] = new DecayingFluorophore[n];
        Random r = new Random();
        double radius = width / 15.0;
        double res = 125.0 / radius;
        double sigma = (0.305f * 602.0 / 1.4) / res;
        double scope = width - 2 * radius;
        GaussianBlur gb = new GaussianBlur();
        for (int i = 0; i < n; i++) {
            double x = r.nextDouble() * scope + radius;
            double y = r.nextDouble() * scope + radius;
            dots[i] = new DecayingFluorophore(x, y, 255.0, 0.05);
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
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\Original_"
                    + indFormat.format(i));
            gb.blurGaussian(image, sigma, sigma, 0.001);
            image.setInterpolationMethod(ImageProcessor.BICUBIC);
            IJ.saveAs(new ImagePlus("", image.resize((int) Math.round(width * res / finalRes))), "TIF",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\BlurredAndScaled_"
                    + indFormat.format(i));
        }
    }
}
