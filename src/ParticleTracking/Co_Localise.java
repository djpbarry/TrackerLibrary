package ParticleTracking;

import IAClasses.IsoGaussian;
import IAClasses.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.*;
import ij.text.TextWindow;
import java.awt.Toolkit;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

public class Co_Localise implements PlugIn {

    protected ImagePlus imp;
    protected ImageStack stack;
    protected final String title = "Colocaliser v1.06";
//    public static final String[] channels = {"Red", "Green", "Blue"};
    protected String resultsHeadings, coordHeadings;
    protected static double coFactor = 2.5;
//    protected static int channel1 = 0, channel2 = 1;
    public static final int RED = 0, GREEN = 1, BLUE = 2;
    protected DecimalFormat numFormat = new DecimalFormat("0.0");
    protected static boolean partialDetect = false;
    protected TextWindow results = null, particleCoords = null;
    protected boolean findTails = false;

//    public static void main(String args[]) {
//        (new Co_Localise(new ImagePlus("C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Co_Localiser_Test\\Co_Localiser_Test.png"))).run(null);
//    }
    public Co_Localise() {
    }

    public Co_Localise(ImagePlus imp) {
        this.imp = imp;
        if (imp != null) {
            this.stack = imp.getStack();
        } else {
            stack = null;
        }
    }

    @Override
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
        if (!(imp.getImageStack().getProcessor(1) instanceof ColorProcessor)) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("Colour (RGB) image required.");
            return;
        }
        if (showDialog()) {
            resultsHeadings = "Image\tChannel 1 (" + UserVariables.channels[UserVariables.getC1Index()]
                    + ") Detections\tColocalised Channel 2 (" + UserVariables.channels[UserVariables.getC2Index()]
                    + ") Detections\t% Colocalisation\t"
                    + "\u0394 (nm)";
            coordHeadings = "C0_X\tC0_Y\tC1_X\tC1_Y";
            UserVariables.setPreProcess(true);
            Timelapse_Analysis analyser = new Timelapse_Analysis(stack);
            analyser.calcParticleRadius(UserVariables.getSpatialRes());
            UserVariables.setnMax(1);
            //Timelapse_Analysis.setGaussianRadius(0.139 / Timelapse_Analysis.getSpatialRes());
            (buildOutput(analyser)).show();
        }
    }

    public boolean showDialog() {
        if (imp == null) {
            Toolkit.getDefaultToolkit().beep();
            IJ.error("No image open.");
            return false;
        }
        boolean valid = false;
        while (!valid) {
            valid = true;
            GenericDialog dialog = new GenericDialog(title, IJ.getInstance());
            dialog.addMessage("Channel 2 will be co-localised with Channel 1.");
            dialog.addChoice("Channel 1:", UserVariables.channels, UserVariables.channels[RED]);
            dialog.addChoice("Channel 2:", UserVariables.channels, UserVariables.channels[GREEN]);
            dialog.addNumericField("Spatial Resolution:", UserVariables.getSpatialRes() * 1000.0, 3, 5, "nm/pixel");
            dialog.addNumericField("Minimum Peak Size (Ch 1):", UserVariables.getChan1MaxThresh(), 3, 5, "");
            dialog.addNumericField("Minimum Peak Size (Ch 2):", UserVariables.getChan2MaxThresh(), 3, 5, "");
            dialog.addNumericField("Curve Fit Tolerance (Ch 1):", UserVariables.getCurveFitTol(), 3, 5, "");
            dialog.addNumericField("Curve Fit Tolerance (Ch 2):", UserVariables.getCurveFitTol(), 3, 5, "");
            dialog.addNumericField("Colocalisation Factor:", coFactor, 3, 5, "");
            dialog.addCheckbox("Include Partial Detections", partialDetect);
            dialog.showDialog();
            if (!dialog.wasCanceled()) {
                UserVariables.setC1Index(dialog.getNextChoiceIndex());
                UserVariables.setC2Index(dialog.getNextChoiceIndex());
                UserVariables.setSpatialRes(dialog.getNextNumber() / 1000.0);
                UserVariables.setChan1MaxThresh(dialog.getNextNumber());
                UserVariables.setChan2MaxThresh(dialog.getNextNumber());
                UserVariables.setCurveFitTol(dialog.getNextNumber());
                /*
                 * Timelapse_Analysis.setC1SigmaTol(sigmaTolC1);
                 * Timelapse_Analysis.setC2SigmaTol(sigmaTolC2);
                 */
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

    public byte[] getPixels(int channel, int frame) {
        if (imp == null) {
            return null;
        }
        int width = imp.getWidth(), height = imp.getHeight();
        int size = width * height;
        ColorProcessor processor = (ColorProcessor) stack.getProcessor(frame + 1);
        byte redPix[] = new byte[size], greenPix[] = new byte[size],
                bluePix[] = new byte[size], emptyPix[] = new byte[size];
        Arrays.fill(emptyPix, (byte) 0);
        processor.getRGB(redPix, greenPix, bluePix);
        switch (channel) {
            case RED:
                return redPix;
            case GREEN:
                return greenPix;
            case BLUE:
                return bluePix;
            default:
                return emptyPix;
        }
    }

    byte[] outPix(ImageProcessor ch1, ImageProcessor ch2, int colour) {
        ImageProcessor fch1 = (new TypeConverter(ch1, true)).convertToByte();
        ImageProcessor fch2 = (new TypeConverter(ch2, true)).convertToByte();
        if (fch1 == null || fch2 == null) {
            return null;
        } else if (UserVariables.getC1Index() == colour) {
            return (byte[]) fch1.getPixels();
        } else if (UserVariables.getC2Index() == colour) {
            return (byte[]) fch2.getPixels();
        } else {
            return (byte[]) (new ByteProcessor(fch1.getWidth(), fch1.getHeight())).getPixels();
        }
    }

    ImagePlus buildOutput(Timelapse_Analysis analyser) {
        if (stack == null) {
            return null;
        }
        if (analyser == null) {
            analyser = new Timelapse_Analysis(stack);
        }
        double displaymax = 0.0;
        int colocalisation, count;
        int width = imp.getWidth(), height = imp.getHeight();
        ImageStack outStack = new ImageStack(width, height);
        double res = UserVariables.getSpatialRes();
        double sepsum;
        for (int i = 0; i < stack.getSize(); i++) {
//            ByteProcessor tailImage = new ByteProcessor(stack.getWidth(), stack.getHeight());
//            tailImage.setPixels(getPixels(channel2, i));
            colocalisation = 0;
            count = 0;
//            tails = 0;
            sepsum = 0.0;
            ParticleArray curves = analyser.findParticles(coFactor, false, i, i);
            FloatProcessor ch1proc = new FloatProcessor(width, height);
            FloatProcessor ch2proc = new FloatProcessor(width, height);
            ArrayList detections = curves.getLevel(0);
            for (int j = 0; j < detections.size(); j++) {
                IsoGaussian c1 = ((Particle) detections.get(j)).getC1Gaussian();
                String coordString = " ";
                if (particleCoords == null) {
                    particleCoords = new TextWindow(title + " Particle Coordinates", coordHeadings, new String(), 1000, 500);
                }
                if (Utils.draw2DGaussian(ch1proc, c1, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(), partialDetect, false)) {
                    if (c1.getMagnitude() > displaymax) {
                        displaymax = c1.getMagnitude();
                    }
                    count++;
//                    if (findTails) {
//                        TailTracer tracer = new TailTracer(tailImage.duplicate(), 0.25, c1.getX() / res,
//                                c1.getY() / res, 2, 0.5, 20, 5);
//                        if (tracer.trace()) {
//                            double[] xTail = tracer.getX();
//                            double[] yTail = tracer.getY();
//                            double[] intens = tracer.getTailIntens();
//                            for (int k = 0; k < xTail.length; k++) {
//                                int x = (int) Math.round(xTail[k]);
//                                int y = (int) Math.round(yTail[k]);
//                                double val = 255.0 * intens[k];
//                                if (val > displaymax) {
//                                    displaymax = val;
//                                }
//                                ch2proc.setValue(val);
//                                ch2proc.drawPixel(x, y);
//                            }
//                            tails++;
//                        }
//                    } else {
                    IsoGaussian c2 = ((Particle) detections.get(j)).getC2Gaussian();
                    if (Utils.draw2DGaussian(ch2proc, c2, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(),
                            partialDetect, false)) {
                        if (c2.getMagnitude() > displaymax) {
                            displaymax = c2.getMagnitude();
                        }
                        colocalisation++;
                        sepsum += Utils.calcDistance(c1.getX(), c1.getY(), c2.getX(), c2.getY());
                        coordString = String.valueOf(c1.getX()) + "\t" + String.valueOf(c1.getY())
                                + "\t" + String.valueOf(c2.getX()) + "\t" + String.valueOf(c2.getY());
//                        }
                    } else {
                        coordString = String.valueOf(c1.getX()) + "\t" + String.valueOf(c1.getY()) + "\t \t ";
                    }
                }
                particleCoords.append(coordString);
            }
            if (results == null) {
                results = new TextWindow(title + " Results", resultsHeadings, new String(), 1000, 500);
            }
            results.append(imp.getTitle() + " (Slice " + i + ")\t" + count + "\t" + colocalisation
                    + "\t" + numFormat.format(100.0 * colocalisation / count)
                    + "\t" + numFormat.format(1000.0 * res * sepsum / count));

            ColorProcessor output = new ColorProcessor(width, height);
            output.setRGB(outPix(ch1proc, ch2proc, RED), outPix(ch1proc, ch2proc, GREEN),
                    outPix(ch1proc, ch2proc, BLUE));
            outStack.addSlice("" + i, output);
        }
        if (results != null) {
            results.append("\n" + toString());
            results.setVisible(true);
        }
        if (particleCoords != null) {
            particleCoords.setVisible(true);
        }
        ImagePlus output = new ImagePlus("Detected Particles", outStack);
        if (displaymax > 255) {
            displaymax = 255;
        }
        output.setDisplayRange(0.0, displaymax);
        return output;
    }
}
