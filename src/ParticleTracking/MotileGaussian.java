/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ParticleTracking;

import Particle.IsoGaussian;
import ij.IJ;
import java.util.Random;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 *
 * @author barry05
 */
public class MotileGaussian extends IsoGaussian {

    boolean persistent, changeState;
    protected double rad, theta;
    private double sens, scale = 1.0, initvel = 0.25;
    Random r = new Random();
    private double D;

    public MotileGaussian(double x0, double y0, double a, double xsig, double ysig, double fit, double sens, boolean persistent, boolean changeState, double D, double initVel) {
        super(x0, y0, a, xsig, ysig, fit);
        this.sens = sens;
        this.persistent = persistent;
        this.changeState = changeState;
//        rad = initvel + r.nextGaussian() * sens;
        this.initvel = initVel;
        rad = initvel;
        theta = r.nextDouble() * 2.0 * Math.PI;
        if (!persistent) {
            rad *= scale;
        }
        this.D = D;
    }

    public void updatePosition() {
        if (persistent) {
            this.x += rad * Math.cos(theta);
            this.y += rad * Math.sin(theta);
        } else {
            updateBrownianCoords();
        }
    }

    void updateBrownianCoords() {
        NormalDistribution nd = new NormalDistribution(0.0, Math.sqrt(4.0 * D));
        this.x += nd.sample();
        this.y += nd.sample();
        D = D * 1.001;
        IJ.log(String.valueOf(D));
    }

    public void updateVelocity() {
        if (persistent) {
//            double inc = r.nextGaussian() * sens;
            double inc = r.nextGaussian() * Math.PI * 2.0;
            if (r.nextBoolean()) {
                inc *= -1.0;
            }
//        if (persistent) {
            theta += inc * sens;
        }
//        } else {
//            theta += inc;
//        }
        if (changeState && r.nextDouble() < 0.2) {
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
                this.fit, this.sens, this.persistent, this.changeState, 0.001, initvel);
    }
}
