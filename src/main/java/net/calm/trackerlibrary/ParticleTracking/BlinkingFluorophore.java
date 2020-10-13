/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.trackerlibrary.ParticleTracking;

import java.util.Random;

/**
 *
 * @author barry05
 */
public class BlinkingFluorophore extends Fluorophore {

    private double onProb;
    private Random rand = new Random();

    public BlinkingFluorophore(double x, double y, double mag, double onProb) {
        super(x, y, mag, 0.05);
        this.onProb = onProb;
    }

    public void updateMag() {
        if (rand.nextDouble() < onProb) {
            currentMag = 255.0;
        } else {
            currentMag = 0.0;
        }
    }

    public double getOnProb() {
        return onProb;
    }
}
