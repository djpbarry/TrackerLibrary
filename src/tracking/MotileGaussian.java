/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tracking;

import IAClasses.IsoGaussian;
import java.util.Random;

/**
 *
 * @author barry05
 */
public class MotileGaussian extends IsoGaussian {

    boolean persistent, changeState;
    private double sens, rad, theta, scale = 1.0, initvel = 0.75;
    Random r = new Random();

    public MotileGaussian(double x0, double y0, double a, double xsig, double ysig,
            double fit, double sens, boolean persistent, boolean changeState) {
        super(x0, y0, a, xsig, ysig, fit);
        this.sens = sens;
        this.persistent = persistent;
        this.changeState = changeState;
        rad = initvel + r.nextGaussian() * sens;
        theta = r.nextGaussian() * 2.0 * Math.PI;
        if (!persistent) {
            rad *= scale;
        }
    }

    public void updatePosition() {
        this.x0 += rad * Math.cos(theta);
        this.y0 += rad * Math.sin(theta);
    }

    public void updateVelocity() {
        double inc = r.nextGaussian() * sens;
        if (r.nextBoolean()) {
            inc *= -1.0;
        }
        if (persistent) {
            rad += inc;
        } else {
            rad += inc * scale;
        }
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
}
