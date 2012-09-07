/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package VirusTracking;

import IAClasses.IsoGaussian;
import ij.process.ImageProcessor;

/**
 *
 * @author barry05
 */
public class RandomDistribution extends IsoGaussian {

    private ImageProcessor pixels;

    public RandomDistribution(ImageProcessor pixels) {
        super();
        this.pixels = pixels;
    }

    public ImageProcessor getPixels() {
        return pixels;
    }
}
