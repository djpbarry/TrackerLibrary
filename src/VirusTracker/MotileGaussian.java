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

    private double sens, xvel, yvel;
    Random r = new Random();

    public MotileGaussian(double x0, double y0, double a, double xsig, double ysig,
            double fit, double sens) {
        super(x0, y0, a, xsig, ysig, fit);
        this.sens = sens;
        if (r.nextBoolean()) {
            xvel = 0.1;
        } else {
            xvel = -0.1;
        }
        if (r.nextBoolean()) {
            yvel = 0.1;
        } else {
            yvel = -0.1;
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
        xvel += inc;
        inc = r.nextGaussian() * sens;
        if (r.nextBoolean()) {
            inc *= -1.0;
        }
        yvel += inc;
    }
}
