/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import AnaMorf.Utilities;
import IAClasses.IsoGaussian;
import IAClasses.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import ij.text.TextWindow;
import java.awt.Rectangle;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author barry05
 */
public class Volume_Analysis extends Timelapse_Analysis {

    private int outputsize = 51, midpoint = (outputsize - 1) / 2;

//    public static void main(String args[]) {
//        File image = Utilities.getFolder(new File("C:\\Users\\barry05\\Desktop\\Tracking Test Sequences"), null);
//        ImageStack stack = Utils.buildStack(image);
//        ImagePlus imp = new ImagePlus("Stack", stack);
//        Volume_Analysis instance = new Volume_Analysis(imp);
//        if (instance.showDialog()) {
//            instance.analyse();
//        }
//        return;
//    }
    public Volume_Analysis() {
        super();
    }

    public Volume_Analysis(ImagePlus imp) {
        super(imp);
        this.imp = imp;
        this.stack = imp.getImageStack();
    }

    public void analyse() {
        if (stack != null) {
            IJ.register(this.getClass());
            int i, count;
            int width = stack.getWidth(), height = stack.getHeight();
            if (width > VIS_SIZE || height > VIS_SIZE) {
                scale = 1.0;
            } else {
                scale = VIS_SIZE / Math.max(width, height);
            }

            findParticles(1.0, true, 0, stack.getSize() - 1);

            TextWindow results = new TextWindow(title + " Results", "X\tY\tFrame\tChannel 1 ("
                    + Co_Localise.channels[colocaliser.getChannel1()]
                    + ")\tChannel 2 (" + Co_Localise.channels[colocaliser.getChannel2()]
                    + ")\tChannel 2 " + '\u03C3' + "x\tChannel 2 " + '\u03C3' + "y\t" + '\u03B8',
                    null, 1000, 500);
            results.append(imp.getTitle() + "\n\n");
            TextWindow resultSummary = new TextWindow(title + " Results Summary",
                    "Particle\tType\t% Colocalisation\tDuration (s)\tDisplacement (" + IJ.micronSymbol
                    + "m)\tVelocity (" + IJ.micronSymbol + "m/s)\tDirectionality\tDiffusion Coefficient ("
                    + IJ.micronSymbol + "m^2/s)" + "\tFractal Dimension"
                    + "\tFluorescence Ratio ("
                    + Co_Localise.channels[colocaliser.getChannel2()] + "/"
                    + Co_Localise.channels[colocaliser.getChannel1()]
                    + ")\tAngle Spread\tStep Spread\tDC\tCurvature\tC2 Fluor Area\tC2 Fluor Skew",
                    null, 1200, 500);
            resultSummary.append(imp.getTitle() + "\n\n");

            int n = trajectories.size();
            for (i = 0; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                if (!(traj.getSize() > minTrajLength && ((traj.getType() == ParticleTrajectory.COLOCAL)
                        || ((traj.getType() == ParticleTrajectory.NON_COLOCAL) && !colocal)))) {
                    trajectories.remove(i);
                    i--;
                    n--;
                }
            }
            n = trajectories.size();
            mapTrajectories(stack, monoChrome, trajectories, scale, spatialRes, xyPartRad, minTrajLength, timeRes, true);
            ArrayList distributions = new ArrayList();
            int cropRad = 4 * xyPartRad + 1;
            for (i = 0, count = 1; i < n; i++) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
                int s = traj.getSize();
                int t = traj.getType();
                if (s > minTrajLength && ((t == ParticleTrajectory.COLOCAL)
                        || ((t == ParticleTrajectory.NON_COLOCAL) && !colocal))) {
                    Particle current = traj.getEnd();
                    double xsum = 0.0, ysum = 0.0;
                    int peakTime = (int) Math.round(traj.getPeakTime() / timeRes);
                    while (current != null) {
                        xsum += current.getX() / spatialRes;
                        ysum += current.getY() / spatialRes;
                        current = current.getLink();
                    }
                    int x = (int) Math.round(xsum / s);
                    int y = (int) Math.round(ysum / s);
                    if (!((x < cropRad) || (y < cropRad) || (stack.getWidth() - x < cropRad)
                            || (stack.getHeight() - y < cropRad))) {
                        ImageStack dist = new ImageStack(cropRad, cropRad);
                        Rectangle roi = new Rectangle(x - 2 * xyPartRad, y - 2 * xyPartRad, cropRad, cropRad);
                        for (int k = peakTime + 1; k <= stack.getSize() && k - peakTime < midpoint + 1; k++) {
                            ImageProcessor currentIP = stack.getProcessor(k);
                            currentIP.setRoi(roi);
                            dist.addSlice("" + count, currentIP.crop());
                        }
                        while (dist.getSize() < midpoint + 1) {
                            dist.addSlice("" + count, new ColorProcessor(cropRad, cropRad));
                        }
                        for (int k = peakTime; k > 0 && peakTime - k < midpoint - 1; k--) {
                            ImageProcessor currentIP = stack.getProcessor(k);
                            currentIP.setRoi(roi);
                            dist.addSlice("" + count, currentIP.crop(), 0);
                        }
                        while (dist.getSize() < outputsize) {
                            dist.addSlice("" + count, new ColorProcessor(cropRad, cropRad), 0);
                        }
                        distributions.add(dist);
                        count++;
                        if (intensPlot) {
                            plotIntensity(i, count);
                        }
                        if (trajPlot) {
                            plotTrajectory(width, height, i, count);
                        }
                        printData(i, resultSummary, count);
                        traj.printTrajectory(count, results, numFormat);
                    }
                }
            }
            double distVals[][] = new double[outputsize][cropRad * cropRad];
            for (int d = 0; d < outputsize; d++) {
                Arrays.fill(distVals[d], 0.0);
            }
            double maxval = -Double.MAX_VALUE;
            for (int m = 0; m < distributions.size(); m++) {
                ImageStack thisStack = (ImageStack) distributions.get(m);
                for (int c = 0; c < thisStack.getSize(); c++) {
                    ImageProcessor slice = (ColorProcessor) thisStack.getProcessor(c + 1);
                    for (int x = 0; x < cropRad; x++) {
                        for (int y = 0; y < cropRad; y++) {
                            distVals[c][x + y * cropRad] += (slice.getPixel(x, y, null))[1];
                            if (distVals[c][x + y * cropRad] > maxval) {
                                maxval = distVals[c][x + y * cropRad];
                            }
                        }
                    }
                }
            }
            ImageStack output = new ImageStack(cropRad, cropRad);
            for (int b = 0; b < distVals.length; b++) {
                output.addSlice("" + b, new FloatProcessor(cropRad, cropRad, distVals[b]));
            }
            ImagePlus outimp = new ImagePlus("Sum of Distributions", output);
            outimp.setDisplayRange(0.0, maxval);
            outimp.show();
            results.append(toString());
            results.setVisible(true);
            resultSummary.setVisible(true);
        }
        return;
    }

    public ParticleArray findParticles(double searchScale, boolean update, int startSlice, int endSlice) {
        if (stack == null) {
            return null;
        }
        xySigEst = (0.21 * 650.0 / 1.4) / (spatialRes * 1000.0);
        xyPartRad = (int) Math.round(2.0 * xySigEst / 0.95);
        int i, noOfImages = stack.getSize(), width = stack.getWidth(), height = stack.getHeight(),
                size = width * height, arraySize = endSlice - startSlice + 1;
        byte c1Pix[] = new byte[size], c2Pix[] = new byte[size],
                c3Pix[] = new byte[size];
        int c1X, c1Y, pSize = 2 * xyPartRad + 1;
        double[] xCoords = new double[pSize];
        double[] yCoords = new double[pSize];
        double[][] pixValues = new double[pSize][pSize];
        ParticleArray particles = new ParticleArray(arraySize);
        for (i = startSlice; i < noOfImages && i <= endSlice; i++) {
            IJ.freeMemory();
            IJ.showStatus("Analysing Frame " + i);
            IJ.showProgress(i, noOfImages);
            if (!monoChrome) {
                ColorProcessor colProc = new ColorProcessor(width, height);
                colProc.setRGB(colocaliser.getPixels(colocaliser.getChannel1(), i),
                        colocaliser.getPixels(colocaliser.getChannel2(), i),
                        colocaliser.getPixels(-1, i));
                ((ColorProcessor) colProc).getRGB(c1Pix, c2Pix, c3Pix);
            } else {
                c1Pix = (byte[]) (new TypeConverter(stack.getProcessor(i + 1).duplicate(), true).convertToByte().getPixels());
                c2Pix = null;
            }
            FloatProcessor chan1Proc = preProcess(new ByteProcessor(width, height, c1Pix, null));
            ByteProcessor thisC1Max = Utils.findLocalMaxima(xyPartRad, xyPartRad, FOREGROUND, chan1Proc, chan1MaxThresh, true);
            for (c1X = 0; c1X < width; c1X++) {
                for (c1Y = 0; c1Y < height; c1Y++) {
                    if (thisC1Max.getPixel(c1X, c1Y) == FOREGROUND) {
                        IsoGaussian c1Gaussian = null;
                        /*
                         * Search for local maxima in green image within
                         * <code>xyPartRad</code> pixels of maxima in red image:
                         */
                        Utils.extractValues(xCoords, yCoords, pixValues, c1X, c1Y, chan1Proc);
                        /*
                         * Remove adjacent Gaussians
                         */
                        IsoGaussianFitter c1GF = new IsoGaussianFitter(xCoords, yCoords, pixValues);
                        c1GF.doFit(xySigEst);
                        //if (c1GF.getXsig() < (c1SigmaTol * xySigEst)) {
                        if (c1GF.getRSquared() > curveFitTol) {
                            c1Gaussian = new IsoGaussian((c1GF.getX0() + c1X - xyPartRad) * spatialRes,
                                    (c1GF.getY0() + c1Y - xyPartRad) * spatialRes, c1GF.getMag(),
                                    c1GF.getXsig(), c1GF.getYsig(), c1GF.getRSquared() - curveFitTol);
                        } else {
                            c1Gaussian = new IsoGaussian(c1X * spatialRes, c1Y * spatialRes, chan1Proc.getPixelValue(c1X, c1Y),
                                    xySigEst, xySigEst, c1GF.getRSquared() - curveFitTol);
                        }
                        /*
                         * A particle has been isolated - trajectories need to
                         * be updated:
                         */
                        if (c1Gaussian != null) {
                            particles.addDetection(i - startSlice, new Particle(i - startSlice, c1Gaussian, null, null, -1));
                        }
                        //}
                    }
                }
            }
        }
        if (update) {
            updateTrajectories(particles, timeRes, trajMaxStep, chan1MaxThresh, hystDiff, spatialRes, true);
        }
        return particles;
    }

    public boolean printData(int particleNumber, TextWindow output, int label) {
        if (trajectories.size() < 1) {
            return false;
        }
        DecimalFormat decFormat = new DecimalFormat("0.000");
        DecimalFormat msdFormat = new DecimalFormat("0.000000");
        ParticleTrajectory traj = (ParticleTrajectory) (trajectories.get(particleNumber));
        if (traj == null) {
            return false;
        }
        traj.smooth();
        traj.calcMSD(label, -1, msdPlot);
        traj.calcAngleSpread();
        traj.calcStepSpread();
        traj.calcDirectionality();
        traj.calcFluorSpread();
        double fluorMaj = Math.max(traj.getxFluorSpread(), traj.getyFluorSpread());
        double fluorMin = Math.min(traj.getxFluorSpread(), traj.getyFluorSpread());
        double fluorArea = Math.PI * 4.0 * fluorMin * fluorMaj * spatialRes * spatialRes;
        double displacement = traj.getDisplacement();
        double duration = traj.getSize() * timeRes;
        int type = traj.getType();
        String trajType = null;
        switch (type) {
            case ParticleTrajectory.COLOCAL:
                trajType = "Colocalised";
                break;
            case ParticleTrajectory.NON_COLOCAL:
                trajType = "Non-Colocalised";
                break;
            case ParticleTrajectory.UNKNOWN:
                trajType = "Unknown";
        }
        output.append(label + "\t" + trajType + "\t"
                + decFormat.format(traj.getDualScore() * 100.0 / traj.getSize()) + "\t"
                + decFormat.format(duration) + "\t"
                + decFormat.format(displacement)
                + "\t" + decFormat.format(displacement / duration) + "\t"
                + decFormat.format(traj.getDirectionality()) + "\t"
                + msdFormat.format(traj.getDiffCoeff()) + "\t"
                + decFormat.format(traj.getBoxCountFD()) + "\t"
                + decFormat.format(traj.getFluorRatio()) + "\t"
                + decFormat.format(traj.getAngleSpread()) + "\t"
                + decFormat.format(traj.getStepSpread()) + "\t"
                + decFormat.format(traj.getLogDC()) + "\t"
                + decFormat.format(traj.getMeanKappa()) + "\t"
                + decFormat.format(fluorArea) + "\t"
                + decFormat.format(fluorMin / fluorMaj));
        return true;
    }
}