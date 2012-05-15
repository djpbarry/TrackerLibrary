/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package VirusTracker;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import java.text.DecimalFormat;

/**
 *
 * @author barry05
 */
public class TestGenerator {

    /*public static void main(String args[]) {
        TestGenerator tg = new TestGenerator();
        tg.generate();
    }*/

    public TestGenerator() {
    }

    public void generate() {
        DecimalFormat indFormat = new DecimalFormat("000");
        int width = 640, height = 480;
        Co_Localise cl = new Co_Localise();
        double res = Timelapse_Analysis.getSpatialRes();
        for (int i = 0; i < 50; i++) {
            ByteProcessor image = new ByteProcessor(width, height);
            IsoGaussian g1 = new IsoGaussian((i * 0.2 + 320) * res, 240.0 * res, 100.0, 2.0*res, 2.0*res, 0.1);
            IsoGaussian g2 = new IsoGaussian((320 - i * 0.2) * res, 240.0 * res, 100.0, 2.0*res, 2.0*res, 0.1);
            cl.draw2DGaussian(image, g1, 0.0);
            cl.draw2DGaussian(image, g2, 0.0);
            IJ.saveAs(new ImagePlus("", image.duplicate()), "PNG",
                    "C:\\Users\\barry05\\Desktop\\Tracking Test Sequences\\Simulation\\"
                    + indFormat.format(i));
        }
    }
}
