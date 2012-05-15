/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package VirusTracker;

import EMSeg.ProgressDialog;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Toolkit;
import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 *
 * @author barry05
 */
public class SuperResAnalyser extends Co_Localise {

    private int scaleFactor = 10;

    public static void main(String args[]) {
        (new SuperResAnalyser(new ImagePlus("C:\\Users\\barry05\\Desktop\\SuperResTest.tif"))).run(null);
    }

    public SuperResAnalyser(ImagePlus imp) {
        super(imp);
        imp.show();
        this.imp = imp;
        if (imp != null) {
            this.stack = imp.getStack();
        } else {
            stack = null;
        }
    }

    public void run(String arg) {
        if (IJ.getInstance() != null) {
            imp = IJ.getImage();
            stack = imp.getImageStack();
        }
        if (imp == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image stack open.");
            return;
        }
        if (showDialog()) {
            headings = "Image\tChannel 1 (" + channels[channel1]
                    + ") Detections\tColocalised Channel 2 (" + channels[channel2]
                    + ") Detections\t% Colocalisation";
            Timelapse_Analysis.setPreProcess(true);
            Timelapse_Analysis analyser = new Timelapse_Analysis(stack);
            analyser.calcParticleRadius();
            //Timelapse_Analysis.setGaussianRadius(0.139 / Timelapse_Analysis.getSpatialRes());
            //IJ.saveAs(buildOutput(analyser), "TIF", "C:\\Users\\barry05\\Desktop\\SuperResTestOutputII.tif");
            buildOutput(analyser);
        }
    }

    public boolean draw2DGaussian(ImageProcessor image, IsoGaussian g, double tol) {
        if (image == null || g == null) {
            return false;
        }
        double res = Timelapse_Analysis.getSpatialRes();
        int x, y, drawRad;
        double x0 = scaleFactor * g.getX() / res;
        double y0 = scaleFactor * g.getY() / res;
        double xSigma = g.getXSigma();
        double value;
        drawRad = (int) Math.round(xSigma * 3.0);
        double fit = g.getFit();
        if (fit < 0.0) {
            fit = 0.0;
        }
        if (g.getMagnitude() > 0.0 && g.getMagnitude() < 255.0) {
            for (x = (int) Math.floor(x0 - drawRad); x <= x0 + drawRad; x++) {
                for (y = (int) Math.floor(y0 - drawRad); y <= y0 + drawRad; y++) {
                    /* The current pixel value is added so as not to "overwrite" other
                    Gaussians in close proximity: */
                    value = g.getMagnitude() * Math.exp(-(((x - x0) * (x - x0))
                            + ((y - y0) * (y - y0))) / (2 * xSigma * xSigma));
                    value += image.getPixelValue(x, y);
                    image.putPixelValue(x, y, value);
                }
            }
        } else {
            return false;
        }
        return true;
    }

    ImagePlus buildOutput(Timelapse_Analysis analyser) {
        if (stack == null) {
            return null;
        }
        if (analyser == null) {
            analyser = new Timelapse_Analysis(stack);
        }
        DecimalFormat format = new DecimalFormat("000");
        analyser.setColocaliser(this);
        int count;
        int width = imp.getWidth() * scaleFactor, height = imp.getHeight() * scaleFactor;
        //ImageStack outStack = new ImageStack(width, height);
        ProgressDialog dialog = new ProgressDialog(null, "Processing...", false, false);
        dialog.setVisible(true);
        FloatProcessor ch1proc = new FloatProcessor(width, height);
        for (int i = 0; i < stack.getSize(); i++) {
            dialog.updateProgress(i, stack.getSize());
            count = 0;
            ParticleArray curves = analyser.findParticles(coFactor, false, i, i);
            //ImagePlus temp = new ImagePlus("", ch1proc);
            //temp.show();
            //temp.setDisplayRange(0.0, 255.0);
            ArrayList detections = curves.getLevel(0);
            for (int j = 0; j < detections.size(); j++) {
                IsoGaussian c1 = ((IsoGaussian[]) detections.get(j))[0];
                if (draw2DGaussian(ch1proc, c1, curveFitC1)) {
                    count++;
                }
                //temp.updateAndDraw();
            }
            IJ.saveAs(new ImagePlus("", ch1proc.duplicate()), "TIF", "C:\\Users\\barry05\\Desktop\\SuperResTest\\Output_" + format.format(i) + ".tif");
            //outStack.addSlice("" + i, ch1proc.duplicate());
        }
        dialog.dispose();
        if (results != null) {
            results.append("\n" + toString());
            results.setVisible(true);
        }
        //ImagePlus output = new ImagePlus("Detected Particles", outStack);
        //output.setDisplayRange(0.0, 255.0);
        return null;
    }
}
