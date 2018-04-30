/*
 * Copyright (C) 2018 David Barry <david.barry at crick dot ac dot uk>
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
package Trajectory.DiffusionAnalysis;

import ij.IJ;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author David Barry <david.barry at crick dot ac dot uk>
 */
public class DiffusionAnalyser {

    private static Plot msdPlot;
    private static ArrayList<DescriptiveStatistics> globalMSD = new ArrayList<>();
    private static String plotLegend = "";
    private double diffCoeff;
    private final float D_SCALING = 4.0f;

    public DiffusionAnalyser() {

    }

    public double calcMSD(int seg, int label, double[][] points, int MIN_POINTS_TO_AVERAGE, double timeRes) {
        int maxLength;
        double xval, yval;
        double[] xPoints = points[0], yPoints = points[1], tPoints = points[2];
        if (xPoints == null) {
            return Double.NaN;
        }
        int length = xPoints.length;
        if (seg > 0) {
            maxLength = seg;
        } else {
            maxLength = length;
        }
        ArrayList<Double> timesteps = new ArrayList<>();
        ArrayList<Double> msd = new ArrayList<>();
        for (int i = 0; i < maxLength; i++) {
            DescriptiveStatistics thisMSD = new DescriptiveStatistics();
            for (int j = 0; i + j < maxLength; j++) {
                int dt = (int) Math.round(tPoints[i + j] - tPoints[j]);
                if (dt == i) {
                    xval = Math.pow(xPoints[i + j] - xPoints[j], 2.0);
                    yval = Math.pow(yPoints[i + j] - yPoints[j], 2.0);
                    thisMSD.addValue(xval + yval);
                } else {
                    IJ.wait(0);
                }
            }
            long N = thisMSD.getN();
            if (N >= MIN_POINTS_TO_AVERAGE) {
                timesteps.add(i / timeRes);
                msd.add(thisMSD.getMean());
            }
        }
        if (!(msd.size() > 0)) {
            return Double.NaN;
        }
        if (msdPlot == null) {
            msdPlot = new Plot("Mean Square Displacement",
                    "Time (s)", "Mean Square Displacement (" + IJ.micronSymbol + "m^2)");
            msdPlot.setLineWidth(3);
        }
        double[] tsA = new double[timesteps.size()];
        double[] msdA = new double[msd.size()];
        for (int i = 0; i < timesteps.size(); i++) {
            tsA[i] = timesteps.get(i);
        }
        for (int i = 0; i < msd.size(); i++) {
            msdA[i] = msd.get(i);
            if (globalMSD.size() <= i) {
                globalMSD.add(new DescriptiveStatistics());
            }
            globalMSD.get(i).addValue(msd.get(i));
        }
        Random r = new Random();
        msdPlot.setColor(new Color(r.nextFloat(), r.nextFloat(), r.nextFloat()));
        msdPlot.addPoints(timesteps, msd, Plot.CONNECTED_CIRCLES);
        msdPlot.setLimitsToFit(false);
        plotLegend = ((plotLegend.concat("Particle ")).concat(String.valueOf(label))).concat("\n");
        msdPlot.addLegend(plotLegend);
        msdPlot.draw();
        msdPlot.show();
        CurveFitter fitter = new CurveFitter(tsA, msdA);
        fitter.doFit(CurveFitter.STRAIGHT_LINE);
        return (fitter.getParams())[1] / D_SCALING;
    }

    public double getDiffCoeff() {
        return diffCoeff;
    }

    public static Plot getMsdPlot() {
        return msdPlot;
    }

}
