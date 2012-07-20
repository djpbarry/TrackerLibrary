/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tracking;

import IAClasses.IsoGaussian;

/**
 *
 * @author barry05
 */
public class NonIsoGaussian extends IsoGaussian {

    private double theta, a, b, c;

    public NonIsoGaussian(double x0, double y0, double a, double xsig, double ysig, double theta, double fit) {
        this.x0 = x0;
        this.y0 = y0;
        this.magnitude = a;
        this.xSigma = xsig;
        this.ySigma = ysig;
        this.fit = fit;
        this.theta = theta;
        this.a = Math.pow(Math.cos(theta), 2.0) / (2.0 * Math.pow(xSigma, 2.0))
                + Math.pow(Math.sin(theta), 2.0) / (2.0 * Math.pow(ySigma, 2.0));
        this.b = -Math.sin(2.0 * theta) / (4.0 * Math.pow(xSigma, 2.0))
                + Math.sin(2.0 * theta) / (4.0 * Math.pow(ySigma, 2.0));
        this.c = Math.pow(Math.cos(theta), 2.0) / (2.0 * Math.pow(ySigma, 2.0))
                + Math.pow(Math.sin(theta), 2.0) / (2.0 * Math.pow(xSigma, 2.0));
    }

    public NonIsoGaussian(NonIsoGaussianFitter fitter, double fitTol) {
        super();
        double p[] = fitter.getParams();
        this.xSigma = p[1];
        this.ySigma = p[2];
        this.magnitude = p[3];
        this.x0 = p[4];
        this.y0 = p[5];
        this.a = Math.pow(Math.cos(p[0]), 2.0) / (2.0 * Math.pow(p[1], 2.0))
                + Math.pow(Math.sin(p[0]), 2.0) / (2.0 * Math.pow(p[2], 2.0));
        this.b = -Math.sin(2.0 * p[0]) / (4.0 * Math.pow(p[1], 2.0))
                + Math.sin(2.0 * p[0]) / (4.0 * Math.pow(p[2], 2.0));
        this.c = Math.pow(Math.cos(p[0]), 2.0) / (2.0 * Math.pow(p[2], 2.0))
                + Math.pow(Math.sin(p[0]), 2.0) / (2.0 * Math.pow(p[1], 2.0));
    }

    public double evaluate(double x, double y) {
        return magnitude * Math.exp(-(a * Math.pow(x - x0, 2.0) + 2 * b * (x - x0)
                * (y - y0) + c * Math.pow(y - y0, 2.0)));
    }

    public double getTheta() {
        return theta;
    }
}
