/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tracking;

import IAClasses.IsoGaussian;
import IAClasses.Gaussian3D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import java.awt.Toolkit;
import java.util.ArrayList;

/**
 *
 * @author barry05
 */
public class Colocalise_3D extends Co_Localise {

    private ImageStack stack1, stack2;
    private String title = "3D Colocaliser v1.0";
    private int particleRadius;
    public final static int FOREGROUND = 255;
    private static double spatialRes = 52.734375 / 1000.0, //Spatial resolution in nm/pixel
            chan1MaxThresh = 150.0, //Threshold value for local maxima in channel 1
            chan2MaxThresh = 150.0, //Threshold value for local maxima in channel 2
            curveFitTol = 0.5d;
    private double xySigEst, zSigEst;

    /*
     * public static void main(String args[]) { (new
     * Colocalise_3D(Utils.buildStack(new
     * File("C:\\Users\\barry05\\Desktop\\Colocalise3DTests\\Test3\\C1")),
     * Utils.buildStack(new
     * File("C:\\Users\\barry05\\Desktop\\Colocalise3DTests\\Test3\\C2")))).run(null);
    }
     */
    public Colocalise_3D() {
        super();
    }

    public Colocalise_3D(ImageStack stack1, ImageStack stack2) {
        super();
        this.stack1 = stack1;
        this.stack2 = stack2;
    }

    public void run(String arg) {
        if (IJ.getInstance() != null) {
            imp = IJ.getImage();
            stack = imp.getImageStack();
        }
        if (stack1 == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image stack open.");
            return;
        }
        if (showDialog()) {
            headings = "Image\tChannel 1 (" + channels[channel1]
                    + ") Detections\tColocalised Channel 2 (" + channels[channel2]
                    + ") Detections\t% Colocalisation";
            (buildOutput()).show();
        }
    }

    public boolean showDialog() {
        if (stack1 == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image open.");
            return false;
        }
        boolean valid = false;
        while (!valid) {
            valid = true;
            GenericDialog dialog = new GenericDialog(title, IJ.getInstance());
            dialog.addMessage("Channel 2 will be co-localised with Channel 1.");
            dialog.addChoice("Channel 1:", channels, channels[RED]);
            dialog.addChoice("Channel 2:", channels, channels[GREEN]);
            dialog.addNumericField("Spatial Resolution:", spatialRes * 1000.0, 1, 5, "nm/pixel");
            dialog.addNumericField("Minimum Peak Size (Ch 1):", chan1MaxThresh, 1, 5, "");
            dialog.addNumericField("Minimum Peak Size (Ch 2):", chan2MaxThresh, 1, 5, "");
            dialog.addNumericField("Curve Fit Tolerance (Ch 1):", curveFitC1, 1, 5, "");
            dialog.addNumericField("Curve Fit Tolerance (Ch 2):", curveFitC2, 1, 5, "");
            dialog.addNumericField("Colocalisation Factor:", coFactor, 1, 5, "");
            dialog.addCheckbox("Include Partial Detections", partialDetect);
            dialog.showDialog();
            if (!dialog.wasCanceled()) {
                channel1 = dialog.getNextChoiceIndex();
                channel2 = dialog.getNextChoiceIndex();
                spatialRes = dialog.getNextNumber() / 1000.0;
                chan1MaxThresh = dialog.getNextNumber();
                chan2MaxThresh = dialog.getNextNumber();
                curveFitC1 = dialog.getNextNumber();
                curveFitC2 = dialog.getNextNumber();
                coFactor = dialog.getNextNumber();
                partialDetect = dialog.getNextBoolean();
                // Check that entries were numeric:
                if (dialog.invalidNumber()) {
                    Toolkit.getDefaultToolkit().beep();
                    IJ.error("Entries must be numeric!");
                    valid = false;
                }
            } else {
                return false;
            }
        }
        return true;
    }

    ImagePlus buildOutput() {
        if (stack1 == null) {
            return null;
        }
        double airyRad = 1.22 * Timelapse_Analysis.LAMBDA / (2.0 * Timelapse_Analysis.NUM_AP); //Airy radius
        particleRadius = (int) Math.round(((2.0 * airyRad / 3.0) + Timelapse_Analysis.virusDiameter)
                / (spatialRes * 2000.0));
        xySigEst = (0.21 * 650.0 / 1.4) / (spatialRes * 1000.0);
        zSigEst = (0.66 * 1.5 * 650.0 / (1.4 * 1.4)) / 50.0;
        double displaymax = 0.0;
        int colocalisation, count;
        int width = stack1.getWidth(), height = stack2.getHeight();
        ImageStack outStack = new ImageStack(width, height);
        ParticleArray curves = findParticles(coFactor);
        ArrayList detections = curves.getLevel(0);
        for (int i = 0; i < stack1.getSize(); i++) {
            colocalisation = 0;
            count = 0;
            FloatProcessor ch1proc = new FloatProcessor(width, height);
            FloatProcessor ch2proc = new FloatProcessor(width, height);
            for (int j = 0; j < detections.size(); j++) {
                Gaussian3D c1 = (Gaussian3D) ((IsoGaussian[]) detections.get(j))[0];
                if (draw2DGaussian(ch1proc, c1, curveFitC1, i)) {
                    if (c1.getMagnitude() > displaymax) {
                        displaymax = c1.getMagnitude();
                    }
                    count++;
                    Gaussian3D c2 = (Gaussian3D) ((IsoGaussian[]) detections.get(j))[1];
                    if (draw2DGaussian(ch2proc, c2, curveFitC2, i)) {
                        if (c2.getMagnitude() > displaymax) {
                            displaymax = c2.getMagnitude();
                        }
                        colocalisation++;
                    }
                }
            }
            if (results == null) {
                results = new TextWindow(title + " Results", headings, null, 1000, 500);
            }
            String imagename = null;
            if (imp != null) {
                imagename = imp.getTitle();
            }
            results.append(imagename + " (Slice " + i + ")\t" + count + "\t" + colocalisation
                    + "\t" + numFormat.format(100.0 * colocalisation / count));
            ColorProcessor output = new ColorProcessor(width, height);
            output.setRGB(outPix(ch1proc, ch2proc, RED), outPix(ch1proc, ch2proc, GREEN),
                    outPix(ch1proc, ch2proc, BLUE));
            outStack.addSlice("" + i, output);
        }
        if (results != null) {
            results.append("\n" + toString());
            results.setVisible(true);
        }
        ImagePlus output = new ImagePlus("Detected Particles", outStack);
        if (displaymax > 255) {
            displaymax = 255;
        }
        //output.setDisplayRange(0.0, displaymax);
        return output;
    }

    public ParticleArray findParticles(double searchScale) {
        if (stack1 == null) {
            return null;
        }
        int width = stack1.getWidth(), height = stack1.getHeight(),
                size = width * height;
        int c1X, c1Y, c1Z, pSize = 2 * particleRadius + 1;
        int c2Points[][];
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[] zCoords = new double[2 * pSize + 1];
        double[][][] pixValues = new double[pSize][pSize][2 * pSize + 1];
        ParticleArray particles = new ParticleArray(1);
        ImageStack thisC1Max = findLocalMaxima3D(particleRadius, particleRadius,
                particleRadius, FOREGROUND, stack1, chan1MaxThresh, true),
                thisC2Max = null;
        if (stack2 != null) {
            thisC2Max = findLocalMaxima3D(particleRadius, particleRadius,
                    particleRadius, FOREGROUND, stack2, chan2MaxThresh, true);
        }
        for (c1Z = 0; c1Z < thisC1Max.getSize(); c1Z++) {
            ImageProcessor currentMaxima = thisC1Max.getProcessor(c1Z + 1);
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (currentMaxima.getPixel(c1X, c1Y) == FOREGROUND) {
                        Gaussian3D c1Gaussian = null, c2Gaussian = null;
                        extractValues3D(xCoords, yCoords, zCoords, pixValues, c1X, c1Y, c1Z, stack1);
                        GaussianFitter3D c1GF = new GaussianFitter3D(xCoords, yCoords, zCoords, pixValues);
                        c1GF.doFit(xySigEst);
                        double c1params[] = c1GF.getParams();
                        if (c1GF.getRSquared() > curveFitTol) {
                            c1Gaussian = new Gaussian3D((c1params[5] + c1X - particleRadius) * spatialRes,
                                    (c1params[6] + c1Y - particleRadius) * spatialRes,
                                    (c1params[7] + c1Z - particleRadius) * spatialRes,
                                    c1params[4], c1params[1], c1params[2], c1params[3],
                                    c1GF.getRSquared() - curveFitTol);
                        } else {
                            c1Gaussian = new Gaussian3D(c1X * spatialRes, c1Y * spatialRes,
                                    c1Z * spatialRes, stack1.getProcessor(c1Z + particleRadius + 1).getPixelValue(c1X, c1Y),
                                    xySigEst, xySigEst, xySigEst,
                                    c1GF.getRSquared() - curveFitTol);
                        }
                        c2Points = searchNeighbourhood(c1X, c1Y, c1Z, (int) Math.round(particleRadius * searchScale),
                                FOREGROUND, thisC2Max);
                        if (c2Points != null) {
                            extractValues3D(xCoords, yCoords, zCoords, pixValues, c2Points[0][0], c2Points[0][1], c2Points[0][2], stack2);
                            GaussianFitter3D c2GF = new GaussianFitter3D(xCoords, yCoords, zCoords, pixValues);
                            c2GF.doFit(xySigEst);
                            double c2params[] = c2GF.getParams();
                            // p[0]=theta, p[1]=xSigma, p[2]=ySigma, p[3]=A, p[4]=x0, p[5]=y0, p[6]=offset
                            c2Gaussian = new Gaussian3D((c2params[5] + c2Points[0][0] - particleRadius) * spatialRes,
                                    (c2params[6] + c2Points[0][1] - particleRadius) * spatialRes,
                                    (c2params[7] + c2Points[0][2] - particleRadius) * spatialRes,
                                    c2params[4], c2params[1], c2params[2], c2params[3],
                                    c2GF.getRSquared() - curveFitTol);
                        }
                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        if (c1Gaussian != null) {
                            particles.addDetection(0, c1Gaussian, c2Gaussian);
                        }

                    }
                }
            }
        }
        return particles;
    }

    public ImageStack findLocalMaxima3D(int kWidth, int kHeight, int kDepth,
            int drawValue, ImageStack stack, double maxThresh, boolean varyBG) {
        if (stack == null) {
            return null;
        }
        int i, j, k, x, y, z, width = stack.getWidth(), height = stack.getHeight(),
                depth = stack.getSize();
        double max, thispix, min;
        ImageStack maxima = new ImageStack(width, height);
        for (z = 0; z < depth; z++) {
            ByteProcessor bproc = new ByteProcessor(width, height);
            bproc.setValue(drawValue);
            for (x = kWidth; x < width - kWidth; x++) {
                for (y = kHeight; y < height - kHeight; y++) {
                    for (min = Double.MAX_VALUE, max = -Double.MAX_VALUE, k = z - kDepth; k <= z + kDepth; k++) {
                        if (k > 0 && k < depth) {
                            ImageProcessor current = stack.getProcessor(k + 1);
                            for (i = x - kWidth; i <= x + kWidth; i++) {
                                for (j = y - kHeight; j <= y + kHeight; j++) {
                                    thispix = current.getPixelValue(i, j);
                                    if ((thispix > max) && !((x == i) && (y == j) && (k == z))) {
                                        max = thispix;
                                    }
                                    if ((thispix < min) && !((x == i) && (y == j) && (k == z))) {
                                        min = thispix;
                                    }
                                }
                            }
                        }
                    }
                    double pix = stack.getProcessor(z + 1).getPixelValue(x, y);
                    double diff;
                    if (varyBG) {
                        diff = pix - min;
                    } else {
                        diff = pix;
                    }
                    if ((pix >= max) && (diff > maxThresh)) {
                        bproc.drawPixel(x, y);
                    }
                }
            }
            maxima.addSlice("", bproc);
        }
        return maxima;
    }

    public static void extractValues3D(double[] xCoords, double[] yCoords, double[] zCoords,
            double[][][] values, int xc, int yc, int zc, ImageStack stack) {
        int i, j, k, x, y, z, w, h, s;
        if (stack == null || xCoords == null || yCoords == null || values == null) {
            return;
        }
        w = stack.getWidth();
        h = stack.getHeight();
        s = stack.getSize();
        if ((xc < 0) || (xc >= w) || (yc < 0) || (yc >= h) || (zc < 0) || (zc >= s)) {
            return;
        }
        int xyradius = (xCoords.length - 1) / 2;
        int zradius = (zCoords.length - 1) / 2;
        for (z = zc - zradius, k = 0; (z <= zc + zradius) && (z < s) && (z >= 0); z++) {
            ImageProcessor current = stack.getProcessor(z + 1);
            zCoords[k] = z;
            for (x = xc - xyradius, i = 0; (x <= xc + xyradius) && (x < w) && (x >= 0); x++) {
                xCoords[i] = x;
                for (y = yc - xyradius, j = 0; (y <= yc + xyradius) && (y < h) && (y >= 0); y++) {
                    yCoords[j] = y;
                    values[i][j][k] = current.getPixelValue(x, y);
                    j++;
                }
                i++;
            }
            k++;
        }
        return;
    }

    public int[][] searchNeighbourhood(int x, int y, int z, int radius, int value,
            ImageStack stack) {
        if (stack == null || x < 0 || x >= stack.getWidth() || y < 0 || y >= stack.getHeight()) {
            return null;
        }
        ArrayList pixels = new ArrayList();
        double currentDist, minDist = Double.MAX_VALUE;
        for (int k = z - radius; k <= z + radius; k++) {
            if (k > 0 && k < stack.getSize()) {
                ImageProcessor currentSlice = stack.getProcessor(k + 1);
                for (int i = x - radius; i <= x + radius; i++) {
                    for (int j = y - radius; j <= y + radius; j++) {
                        if (currentSlice.getPixel(i, j) == value) {
                            currentDist = calcDistance(i, j, k, x, y, z);
                            if (currentDist < minDist) {
                                double p[] = {i, j, k, currentDist};
                                pixels.add(p);
                                minDist = currentDist;
                            }
                        }
                    }
                }
            }
        }
        if (pixels.size() > 0) {
            int points[][] = new int[pixels.size()][3];
            int l = 0;
            while (!pixels.isEmpty()) {
                double currentMin = Double.MAX_VALUE;
                int minIndex = -1;
                for (int m = 0; m < pixels.size(); m++) {
                    double current[] = (double[]) pixels.get(m);
                    if (current[3] < currentMin) {
                        currentMin = current[3];
                        minIndex = m;
                    }
                }
                points[l][0] = (int) (Math.round(((double[]) pixels.get(minIndex))[0]));
                points[l][1] = (int) (Math.round(((double[]) pixels.get(minIndex))[1]));
                points[l][2] = (int) (Math.round(((double[]) pixels.get(minIndex))[2]));
                pixels.remove(minIndex);
                l++;
            }
            return points;
        }
        return null;
    }

    public double calcDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
        return Math.sqrt(Math.pow((x2 - x1), 2.0) + Math.pow((y2 - y1), 2.0) + Math.pow((z2 - z1), 2.0));
    }

    public boolean draw2DGaussian(ImageProcessor image, Gaussian3D g, double tol, int z) {
        if (image == null || g == null || (!partialDetect && (g.getFit() < -tol))) {
            return false;
        }
        int x, y, drawRad;
        double x0 = g.getX() / spatialRes;
        double y0 = g.getY() / spatialRes;
        double z0 = g.getZ0() / spatialRes;
        double xSigma = g.getXSigma(), ySigma = g.getYSigma(), zSigma = g.getzSigma();
        double value;
        drawRad = (int) Math.round(xSigma * 3.0);
        if (g.getFit() < -tol) {
            image.setColor(100);
            image.drawOval((int) Math.round(x0 - drawRad), (int) Math.round(y0 - drawRad),
                    2 * drawRad + 1, 2 * drawRad + 1);
        } else {
            for (x = (int) Math.floor(x0 - drawRad); x <= x0 + drawRad; x++) {
                for (y = (int) Math.floor(y0 - drawRad); y <= y0 + drawRad; y++) {
                    /*
                     * The current pixel value is added so as not to "overwrite"
                     * other Gaussians in close proximity:
                     */
                    value = g.getMagnitude() * Math.exp(-(Math.pow(x - x0, 2.0) / (2.0 * xSigma * xSigma)
                            + Math.pow(y - y0, 2.0) / (2.0 * ySigma * ySigma)
                            + Math.pow(z - z0, 2.0) / (2.0 * zSigma * zSigma)));
                    value += image.getPixelValue(x, y);
                    image.putPixelValue(x, y, value);
                }
            }
        }
        return true;
    }
}
