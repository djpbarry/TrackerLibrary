/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.trackerlibrary.ParticleTracking;

import ij.process.ImageProcessor;
import net.calm.iaclasslibrary.Particle.IsoGaussian;

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
