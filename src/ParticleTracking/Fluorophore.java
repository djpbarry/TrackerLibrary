/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import java.util.Random;

/**
 *
 * @author barry05
 */
public class Fluorophore {

    protected double x;
    protected double y;
    protected double initialMag;
    protected double currentMag;
    protected double thresh;
    private Random rand;

    public Fluorophore(double x, double y, double mag, double thresh) {
        this.x = x;
        this.y = y;
        this.initialMag = mag;
        this.currentMag = mag;
        this.thresh = thresh;
        rand = new Random();
    }

    public void updateMag() {
        if (rand.nextDouble() < thresh) {
            currentMag = 0.0;
        }
    }

    public void updateMag(double newMag) {
        currentMag = newMag;
    }

    public double getCurrentMag() {
        return currentMag;
    }

    public double getInitialMag() {
        return initialMag;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }
}
