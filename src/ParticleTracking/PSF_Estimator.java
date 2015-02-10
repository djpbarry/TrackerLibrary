/*
 * Copyright (C) 2015 David Barry <david.barry at cancer.org.uk>
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

import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import static ParticleTracking.Colocalisation_Analysis.coFactor;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.text.TextWindow;
import java.util.ArrayList;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class PSF_Estimator extends Colocalisation_Analysis {

    public static void main(String args[]) {
        ImagePlus inputs[] = new ImagePlus[2];
        inputs[0] = new ImagePlus("C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Co_Localiser_Test\\Co_Localiser_Test_Red.png");
        inputs[1] = new ImagePlus("C:\\Users\\barry05\\Desktop\\Test_Data_Sets\\Co_Localiser_Test\\Co_Localiser_Test_Green.png");
        (new PSF_Estimator(inputs)).run(null);
    }

    public PSF_Estimator(ImagePlus[] inputs) {
        super(inputs);
    }

    ImagePlus buildOutput(Analyse_Movie analyser) {
        if (stacks == null) {
            return null;
        }
        if (analyser == null) {
            analyser = new Analyse_Movie(stacks);
        }
        int colocalisation, count;
        int width = stacks[0].getWidth(), height = stacks[0].getHeight();
        ImageStack outStack = new ImageStack(width, height);
        double sepsum;
        ProgressDialog progress = new ProgressDialog(null, "Analysing Stacks...", false, title, false);
        progress.setVisible(true);
        for (int i = 0; i < stacks[0].getSize(); i++) {
            progress.updateProgress(i, stacks[0].getSize());
            colocalisation = 0;
            count = 0;
            sepsum = 0.0;
            ParticleArray curves = analyser.findParticles(coFactor, false, i, i, UserVariables.getCurveFitTol(), stacks[0], stacks[1], false, true);
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
                    count++;
                    IsoGaussian c2 = ((Particle) detections.get(j)).getC2Gaussian();
                    if (Utils.draw2DGaussian(ch2proc, c2, UserVariables.getCurveFitTol(), UserVariables.getSpatialRes(),
                            partialDetect, false)) {
                        colocalisation++;
                        sepsum += Utils.calcDistance(c1.getX(), c1.getY(), c2.getX(), c2.getY());
                        coordString = String.valueOf(c1.getX()) + "\t" + String.valueOf(c1.getY())
                                + "\t" + String.valueOf(c2.getX()) + "\t" + String.valueOf(c2.getY());
                    } else {
                        coordString = String.valueOf(c1.getX()) + "\t" + String.valueOf(c1.getY()) + "\t \t ";
                    }
                }
                particleCoords.append(coordString);
            }
            if (results == null) {
                results = new TextWindow(title + " Results", resultsHeadings, new String(), 1000, 500);
            }
            results.append("Slice " + i + "\t" + count + "\t" + colocalisation
                    + "\t" + numFormat.format(100.0 * colocalisation / count)
                    + "\t" + numFormat.format(1000.0 * UserVariables.getSpatialRes() * sepsum / count));

            ColorProcessor output = new ColorProcessor(width, height);
            output.setRGB(outPix(ch1proc, ch2proc, RED), outPix(ch1proc, ch2proc, GREEN),
                    outPix(ch1proc, ch2proc, BLUE));
            outStack.addSlice("" + i, output);
        }
        progress.dispose();
        if (results != null) {
            results.append("\n" + toString());
            results.setVisible(true);
        }
        if (particleCoords != null) {
            particleCoords.setVisible(true);
        }
        ImagePlus output = new ImagePlus("Detected Particles", outStack);
        return output;
    }

}
