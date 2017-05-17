/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import Particle.IsoGaussian;
import java.util.Random;

/**
 *
 * @author barry05
 */
public class MotileGaussian extends IsoGaussian {

    boolean persistent, changeState;
    protected double rad, theta;
    private double sens, scale = 1.0, initvel = 0.5;
    Random r = new Random();

    public MotileGaussian(double x0, double y0, double a, double xsig, double ysig,
            double fit, double sens, boolean persistent, boolean changeState) {
        super(x0, y0, a, xsig, ysig, fit);
        this.sens = sens;
        this.persistent = persistent;
        this.changeState = changeState;
//        rad = initvel + r.nextGaussian() * sens;
        rad = initvel;
        theta = r.nextDouble() * 2.0 * Math.PI;
        if (!persistent) {
            rad *= scale;
        }
    }

    public void updatePosition() {
        this.x += rad * Math.cos(theta);
        this.y += rad * Math.sin(theta);
    }

    public void updateVelocity() {
        double inc = r.nextGaussian() * sens;
//        if (r.nextBoolean()) {
//            inc *= -1.0;
//        }
//        if (persistent) {
//            rad += inc;
//        } else {
//            rad += inc * scale;
//        }
        inc = r.nextGaussian() * Math.PI * 2.0;
        if (r.nextBoolean()) {
            inc *= -1.0;
        }
        if (persistent) {
            theta += inc * sens;
        } else {
            theta += inc;
        }
        if (changeState && r.nextDouble() < 0.05) {
            persistent = !persistent;
        }
    }

    public double[] projectPosition(boolean positive, double dist) {
        double projectedPos[] = new double[2];
        int pos = positive ? 1 : -1;
        projectedPos[0] = this.x + dist * Math.cos(theta) * pos;
        projectedPos[1] = this.y + dist * Math.sin(theta) * pos;
        return projectedPos;
    }

    public Object clone() {
        return new MotileGaussian(this.x, this.y, this.magnitude, this.xSigma, this.ySigma,
                this.fit, this.sens, this.persistent, this.changeState);
    }
}
