/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tracking;

import IAClasses.IsoGaussian;
import ij.IJ;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import java.text.DecimalFormat;
import java.util.Random;

/**
 *
 * @author barry05
 */
public class TestGenerator {

//    public static void main(String args[]) {
//        TestGenerator tg = new TestGenerator();
//        tg.generateMulti(1, 250, 250, 100);
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
        DecimalFormat indFormat = new DecimalFormat("000");
        Co_Localise cl = new Co_Localise();
        double res = Timelapse_Analysis.getSpatialRes();
        MotileGaussian particles[] = new MotileGaussian[n];
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            particles[i] = new MotileGaussian(width * res * r.nextDouble(), height * res * r.nextDouble(),
                    255.0, 2.0, 2.0, 0.1, 0.02, false, false);
        }
        for (int i = 0; i < length; i++) {
            ByteProcessor image = new ByteProcessor(width, height);
            image.setColor(255);
            for (int j = 0; j < n; j++) {
                cl.draw2DGaussian(image, particles[j], 0.0);
                particles[j].updateVelocity();
                particles[j].updatePosition();
                double x = particles[j].getX() / res;
                double y = particles[j].getY() / res;
                if (x < -2.0 * particles[j].getXSigma()
                        || x > width + 2.0 * particles[j].getXSigma()
                        || y < -2.0 * particles[j].getYSigma()
                        || y > height + 2.0 * particles[j].getYSigma()) {
                    particles[j] = new MotileGaussian(width * res * r.nextDouble(), height * res * r.nextDouble(),
                            255.0, 2.0, 2.0, 0.1, 0.02, false, false);
                    totalcount++;
                }
            }
            IJ.saveAs(new ImagePlus("", image.duplicate()), "PNG",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\"
                    + indFormat.format(i));
            System.out.println("Frame:\t" + i + "\tTotal Count:\t" + totalcount);
        }
    }
}
