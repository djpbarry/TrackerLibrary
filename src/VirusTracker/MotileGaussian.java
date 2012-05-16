/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package VirusTracker;

import java.util.Random;

/**
 *
 * @author barry05
 */
public class MotileGaussian extends IsoGaussian {

    boolean persistent, changeState;
    private double sens, xvel, yvel;
    Random r = new Random();

    public MotileGaussian(double x0, double y0, double a, double xsig, double ysig,
            double fit, double sens, boolean persistent, boolean changeState) {
        super(x0, y0, a, xsig, ysig, fit);
        this.sens = sens;
        this.persistent = persistent;
        this.changeState = changeState;
        if (persistent) {
            xvel = yvel = 0.1;
            if (r.nextBoolean()) {
                xvel *= -1.0;
            }
            if (r.nextBoolean()) {
                yvel *= -1.0;
            }
        }
    }

    public void updatePosition() {
        this.x0 += xvel;
        this.y0 += yvel;
    }

    public void updateVelocity() {
        double inc = r.nextGaussian() * sens;
        if (r.nextBoolean()) {
            inc *= -1.0;
        }
        if (persistent) {
            xvel += inc;
        } else {
            xvel = inc;
        }
        inc = r.nextGaussian() * sens;
        if (r.nextBoolean()) {
            inc *= -1.0;
        }
        if (persistent) {
            yvel += inc;
        } else {
            yvel = inc;
        }
        if (changeState && r.nextDouble() < 0.05) {
            persistent = !persistent;
        }
    }
}
