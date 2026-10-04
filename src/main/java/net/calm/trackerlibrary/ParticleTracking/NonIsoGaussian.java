/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.calm.trackerlibrary.ParticleTracking;

import ij.process.FloatProcessor;
import net.calm.iaclasslibrary.Math.Optimisation.NonIsoGaussianFitter;
import net.calm.iaclasslibrary.Particle.IsoGaussian;

/**
 *
 * @author barry05
 */
public class NonIsoGaussian extends IsoGaussian {

    private double theta, a, b, c;

    public NonIsoGaussian(double x0, double y0, double a, double xsig, double ysig, double theta, double fit) {
        this.x = x0;
        this.y = y0;
        this.magnitude = a;
        this.xSigma = xsig;
        this.ySigma = ysig;
        this.fit = fit;
        this.theta = theta;
        double cosTheta = Math.cos(theta);
        double sinTheta = Math.sin(theta);
        double sin2Theta = Math.sin(2.0 * theta);
        double xSigmaSq = xSigma * xSigma;
        double ySigmaSq = ySigma * ySigma;
        this.a = cosTheta * cosTheta / (2.0 * xSigmaSq) + sinTheta * sinTheta / (2.0 * ySigmaSq);
        this.b = -sin2Theta / (4.0 * xSigmaSq) + sin2Theta / (4.0 * ySigmaSq);
        this.c = cosTheta * cosTheta / (2.0 * ySigmaSq) + sinTheta * sinTheta / (2.0 * xSigmaSq);
    }

    public NonIsoGaussian(NonIsoGaussianFitter fitter) {
        super();
        double p[] = fitter.getParams();
        this.xSigma = p[1];
        this.ySigma = p[2];
        this.magnitude = p[3];
        this.x = p[4];
        this.y = p[5];
        double cosP = Math.cos(p[0]);
        double sinP = Math.sin(p[0]);
        double sin2P = Math.sin(2.0 * p[0]);
        double p1Sq = p[1] * p[1];
        double p2Sq = p[2] * p[2];
        this.a = cosP * cosP / (2.0 * p1Sq) + sinP * sinP / (2.0 * p2Sq);
        this.b = -sin2P / (4.0 * p1Sq) + sin2P / (4.0 * p2Sq);
        this.c = cosP * cosP / (2.0 * p2Sq) + sinP * sinP / (2.0 * p1Sq);
    }

    public double evaluate(double x, double y) {
        double dx = x - this.x;
        double dy = y - this.y;
        return magnitude * Math.exp(-(a * dx * dx + 2 * b * dx * dy + c * dy * dy));
    }

    public double getTheta() {
        return theta;
    }

    public double getX0() {
        return x;
    }

    public double getxSigma() {
        return xSigma;
    }

    public double getY0() {
        return y;
    }

    public double getySigma() {
        return ySigma;
    }

    public void draw(FloatProcessor image, double res) {
        int width = image.getWidth();
        int height = image.getHeight();
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                image.putPixelValue(i, j, evaluate(i * res, j * res));
            }
        }
    }
}
