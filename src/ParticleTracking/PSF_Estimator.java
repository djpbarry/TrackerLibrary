/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import IAClasses.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.*;
import ij.text.TextWindow;
import java.awt.Toolkit;

/**
 *
 * @author barry05
 */
public class PSF_Estimator extends Timelapse_Analysis {

    private double lambda = 488.0d;
    private TextWindow results;
    private String psfTitle = "PSF Estimator v1.0";

//    public static void main(String args[]) {
//        File image = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\Test Data Sets\\PSF Estimator\\Test 1"), null);
//        ImageStack stack = Utils.buildStack(image);
//        ImagePlus imp = new ImagePlus("Stack", stack);
//        PSF_Estimator instance = new PSF_Estimator(imp);
//        if (instance.showDialog()) {
//            instance.analyse();
//        }
//        return;
//    }

    public PSF_Estimator() {
        super();
    }

    public PSF_Estimator(ImagePlus imp) {
        super(imp);
        this.imp = imp;
        this.stack = imp.getImageStack();
    }

    public void analyse() {
        stack = imp.getImageStack();
        if (stack != null) {
            IJ.register(this.getClass());
            results = new TextWindow(psfTitle + " Results", "frame\tx\ty\tA\tsigma_x\tsigma_y\tR^2",
                    new String(), 1000, 500);
            results.append(imp.getTitle() + "\n\n");
            findParticles(1.0, true, 0, stack.getSize() - 1, UserVariables.getCurveFitTol());
            results.setVisible(true);
        }
        return;
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice, double fitTol) {
        if (stack == null) {
            return null;
        }
        xySigEst = (0.21 * lambda / 1.4) / (UserVariables.getSpatialRes() * 1000.0);
        xyPartRad = (int) Math.round(2.0 * xySigEst / 0.95);
        int i, noOfImages = stack.getSize(), width = stack.getWidth(), height = stack.getHeight();
        byte pix[];
        int x, y, pSize = 2 * xyPartRad + 1;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            IJ.showStatus("Analysing Frame " + i);
            IJ.showProgress(i, noOfImages);
            pix = (byte[]) (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte().getPixels());
            FloatProcessor floatProc = preProcess(new ByteProcessor(width, height, pix, null));
            ByteProcessor maxima = Utils.findLocalMaxima(xyPartRad, xyPartRad, FOREGROUND, floatProc, UserVariables.getChan1MaxThresh(), true);
            for (x = 0; x < width; x++) {
                for (y = 0; y < height; y++) {
                    if (maxima.getPixel(x, y) == FOREGROUND) {
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, x, y, floatProc);
                        /*
                         * Remove adjacent Gaussians
                         */
                        NonIsoGaussianFitter fitter = new NonIsoGaussianFitter(xCoords, yCoords, pixValues);
                        fitter.doFit(xySigEst);
                        //if (c1GF.getXsig() < (c1SigmaTol * xySigEst)) {
                        NonIsoGaussian gaussian = new NonIsoGaussian(fitter, UserVariables.getCurveFitTol());
                        results.append(i + "\t" + x + "\t" + y + "\t"
                                + numFormat.format(gaussian.getMagnitude())
                                + "\t" + numFormat.format(gaussian.getxSigma() * UserVariables.getSpatialRes() * 1000.0)
                                + "\t" + numFormat.format(gaussian.getySigma() * UserVariables.getSpatialRes() * 1000.0)
                                + "\t" + numFormat.format(fitter.getRSquared()));
                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        //}
                    }
                }
            }
        }
        return null;
    }

    public boolean showDialog() {
        if (imp == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image stack open.");
            return false;
        }
        GenericDialog gd = new GenericDialog(psfTitle);
        gd.addNumericField("Spatial Resolution", UserVariables.getSpatialRes() * 1000.0, 5, 5, "nm");
        gd.addNumericField("Peak Threshold", UserVariables.getChan1MaxThresh(), 5);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        UserVariables.setSpatialRes(gd.getNextNumber() / 1000.0);
        UserVariables.setChan1MaxThresh(gd.getNextNumber());
        return true;
    }
}
