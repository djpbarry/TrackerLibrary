/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package VirusTracker;

import ij.IJ;
import ij.measure.CurveFitter;
import ij.plugin.filter.GaussianBlur;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.TypeConverter;
import java.awt.Rectangle;
import java.util.ArrayList;

/**
 *
 * @author barry05
 */
public class TailTracer {

    ArrayList<Double> x = new ArrayList<Double>();
    ArrayList<Double> y = new ArrayList<Double>();
    ArrayList<Double> tailIntens = new ArrayList<Double>();
    ImageProcessor im;
    double tf;
    double x1;
    double y1;
    int pn;
    double dxy;
    int nlines;
    int npoints;

    public TailTracer() {
    }

    public TailTracer(ImageProcessor im, double tf, double x1, double y1,
            int pn, double dxy, int nlines, int npoints) {
        this.im = im;
        this.tf = tf;
        this.x1 = x1;
        this.y1 = y1;
        this.pn = pn;
        this.dxy = dxy;
        this.nlines = nlines;
        this.npoints = npoints;
    }

    double[] normalizedVector(double x, double y) {
// Normalized vector
        double sn[] = new double[2];
        double n = Math.pow(x * x + y * y, 0.5);
        sn[0] = x / n;
        sn[1] = y / n;
        return sn;
    }

    double[] normalVector(double x1, double y1, double x2, double y2) {
        // Normal vector
        double n[] = new double[2];
        double qx = x2 - x1;
        double qy = y2 - y1;
        double pmag = Math.pow(qx * qx + qy * qy, 0.5);
        n[0] = qy / pmag;
        n[1] = -qx / pmag;
        return n;
    }

    double[] candidatesIntensityCentre(double[] xc, double[] yc) {
        double xcc = 0.0, ycc = 0.0;
        int count = 0;
        for (int i = 0; i < xc.length; i++) {
            double it = im.getInterpolatedValue(xc[i], yc[i]);
            if (it > tf) {
                xcc += xc[i];
                ycc += yc[i];
                count++;
            }
        }
        if (count < 1) {
            return null;
        }
        double cn[] = new double[2];
        cn[0] = xcc / count;
        cn[1] = ycc / count;
        return cn;
    }

    double[] intersections(double m1, double c1, double m2, double c2) {
        double point[] = new double[2];
        point[0] = (c2 - c1) / (m1 - m2);
        point[1] = m1 * point[0] + c1;
        return point;
    }

    double[] intersections2(double a, double b, double c2,
            double x1, double x2, double y1, double y2) {
        double point1[] = new double[2];
        double point2[] = new double[2];
        double m = (y2 - y1) / (x2 - x1);
        double c1 = y1 - m * x1;
        point1[0] = (m - b + Math.pow((b - m) * (b - m) - 4 * a * (c2 - c1), 0.5)) / (2 * a);
        point2[0] = (m - b - Math.pow((b - m) * (b - m) - 4 * a * (c2 - c1), 0.5)) / (2 * a);
        point1[1] = m * point1[0] + c1;
        point2[1] = m * point2[0] + c1;
        Rectangle r;
        int ip1x = (int) Math.round(x1);
        int ip2x = (int) Math.round(x2);
        int ip1y = (int) Math.round(y1);
        int ip2y = (int) Math.round(y2);
        if (x1 < x2) {
            if (y1 < y2) {
                r = new Rectangle(ip1x, ip1y, ip2x - ip1x, ip2y - ip1y);
            } else {
                r = new Rectangle(ip1x, ip2y, ip2x - ip1x, ip1y - ip2y);
            }
        } else {
            if (y1 < y2) {
                r = new Rectangle(ip2x, ip1y, ip1x - ip2x, ip2y - ip1y);
            } else {
                r = new Rectangle(ip2x, ip2y, ip1x - ip2x, ip1y - ip2y);
            }
        }
        if (r.contains(point1[0], point1[1])) {
            return point1;
        } else {
            return point2;
        }
    }

    public boolean trace() {
        im = (new TypeConverter(im, true)).convertToFloat(null);
        (new GaussianBlur()).blur(im, 2.0);
        int ne = 10 * nlines;
        ArrayList<Double[][]> hc = new ArrayList<Double[][]>();
        for (int i = 0; i < nlines * 2 - 1; i++) {
            hc.add(null);
        }
        int lpix = nlines; // index of current candidate point
        // Normalize
        ImageStatistics stats = im.getStatistics();
        im.subtract(stats.min);
        im.multiply(1.0 / (stats.max - stats.min));
        stats = im.getStatistics();
        // Starting vector
        double xy2[] = findStartingVector();
        if (xy2 == null) {
            return false;
        }
        double xs = xy2[0] - x1;
        double ys = xy2[1] - y1;
        double sn[] = normalizedVector(xs, ys);
        for (int i = 0; i < pn; i++) {
            double xcoord = sn[0] * dxy * (i + 1) + x1;
            double ycoord = sn[0] * dxy * (i + 1) + y1;
            x.add(new Double(xcoord));
            y.add(new Double(ycoord));
            tailIntens.add(new Double(im.getInterpolatedValue(xcoord, ycoord)));
        }
        // Tracing
        int stop = 0;
        while (stop == 0) {
            double xl[][] = new double[2 * nlines + 1][2 * npoints + 1];
            double yl[][] = new double[2 * nlines + 1][2 * npoints + 1];
            int N = 0;
            int s = x.size();
            if (s == 21) {
                IJ.wait(2);
            }
            // Create lines  ||||o||||
            for (int k = s - nlines; k <= s + nlines; k++) {
                double nr[], n[];
                int ix;
                if (k > s) {
                    ix = s - 1;
                } else if (k <= 2) {
                    ix = 1;
                } else {
                    ix = k - 1;
                }
                double xn1 = ((Double) x.get(ix)).doubleValue();
                double xnm1 = ((Double) x.get(ix - 1)).doubleValue();
                double yn1 = ((Double) y.get(ix)).doubleValue();
                double ynm1 = ((Double) y.get(ix - 1)).doubleValue();
                nr = normalizedVector(xn1 - xnm1, yn1 - ynm1);
                n = normalVector(xn1, yn1, xnm1, ynm1);
                // -o-o-o-o-o-o-o-
                int m = 0;
                for (int j = -npoints; j <= npoints; j++) {
                    if (k > s) {
                        xl[N][m] = xn1 + (k - s) * dxy * nr[0] + j * dxy * n[0];
                        yl[N][m] = yn1 + (k - s) * dxy * nr[1] + j * dxy * n[1];
                    } else if (k <= 2) {
                        xl[N][m] = xn1 + (k - 2) * dxy * nr[0] + j * dxy * n[0];
                        yl[N][m] = yn1 + (k - 2) * dxy * nr[1] + j * dxy * n[1];
                    } else {
                        xl[N][m] = xn1 + j * dxy * n[0];
                        yl[N][m] = yn1 + j * dxy * n[1];
                    }
                    m++;
                }
                N++;
            }
            double temp[][] = new double[xl.length][2];
            //Intensity for each line
            int iccount = 0;
            for (int i = 0; i < 2 * nlines + 1; i++) { //size(xl,2)
                temp[i] = candidatesIntensityCentre(xl[i], yl[i]);
                if (temp[i] != null) {
                    iccount++;
                }
            }
            double ic[][];
            if (iccount > 0) {
                ic = new double[iccount][2];
            } else {
                ic = null;
            }
            iccount = 0;
            for (int i = 0; i < xl.length; i++) {
                if (temp[i] != null) {
                    ic[iccount][0] = temp[i][0];
                    ic[iccount][1] = temp[i][1];
                    iccount++;
                }
            }
            // Fit polyline
            // if just started - > shipt points
            if (s == pn) {
                // Indexes of x and y to be used for fitting
                int xyix[];
                if (s < nlines + 1) {
                    xyix = new int[s];
                    for (int i = 0; i < s; i++) {
                        xyix[i] = i;
                    }
                } else {
                    for (int i = s - nlines; i < s; i++) {
                        xyix = new int[2 * s - nlines];
                        xyix[i] = i;
                    }
                }
                // intensity centers + xy
                int icl;
                if (ic != null) {
                    icl = ic.length;
                } else {
                    icl = 0;
                }
                double xw[] = new double[icl + x.size()];
                double yw[] = new double[icl + y.size()];
                double xmax = -Double.MAX_VALUE, xmin = -xmax, ymax = xmax, ymin = xmin;
                for (int i = 0; i < xw.length; i++) {
                    if (i < icl) {
                        xw[i] = ic[i][0];
                        yw[i] = ic[i][1];
                    } else {
                        xw[i] = ((Double) x.get(i - icl)).doubleValue();
                        yw[i] = ((Double) y.get(i - icl)).doubleValue();
                    }
                    if (xw[i] > xmax) {
                        xmax = xw[i];
                    } else if (xw[i] < xmin) {
                        xmin = xw[i];
                    }
                    if (yw[i] > ymax) {
                        ymax = yw[i];
                    } else if (yw[i] < ymin) {
                        ymin = yw[i];
                    }
                }
                // Fit polyline
                double xfit[], yfit[];
                double p[];
                CurveFitter cf = new CurveFitter(xw, yw);
                cf.doFit(CurveFitter.POLY2);
                p = cf.getParams();
                xfit = new double[ne];
                yfit = new double[ne];
                for (int i = 0; i < ne; i++) {
                    xfit[i] = xmin + i * (xmax - xmin) / ne;
                    yfit[i] = cf.f(p, xfit[i]);
                }
                int k = 0;
                for (int i = lpix - pn + 1; i <= lpix + 1; i++) { // was no +1
                    double point[] = intersections2(p[2], p[1], p[0], xl[i][0],
                            xl[i][2 * npoints], yl[i][0], yl[i][2 * npoints]);
                    if (k < x.size()) {
                        x.set(k, new Double(point[0]));
                        y.set(k, new Double(point[1]));
                        tailIntens.set(k, new Double(im.getInterpolatedValue(point[0], point[1])));
                    } else {
                        x.add(new Double(point[0]));
                        y.add(new Double(point[1]));
                        tailIntens.add(new Double(im.getInterpolatedValue(point[0], point[1])));
                    }
                    k++;
                }
                // Intersection of polyline and lines
                Double array[][];
                for (int i = 1; i < 2 * nlines; i++) {
                    double point[] = intersections2(p[2], p[1], p[0], xl[i][2 * npoints],
                            xl[i][0], yl[i][2 * npoints], yl[i][0]);
                    Double hcp[][] = (Double[][]) hc.get(i - 1);
                    if (hcp != null) {
                        array = new Double[hcp.length + 1][2];
                        for (int j = 0; j < hcp.length; j++) {
                            array[j][0] = hcp[j][0];
                            array[j][1] = hcp[j][1];
                        }
                        array[hcp.length][0] = new Double(point[0]);
                        array[hcp.length][1] = new Double(point[1]);
                    } else {
                        array = new Double[1][2];
                        array[0][0] = new Double(point[0]);
                        array[0][1] = new Double(point[1]);
                    }
                    hc.set(i - 1, array);
                }
                // for the rest of points        
            } else {
                // Shift the histogram table hc(1:n-1) = hc(2:n) and hc(n) = [];
                int h = hc.size();
                for (int i = 0; i < h - 1; i++) {
                    hc.set(i, hc.get(i + 1));
                }
                hc.set(h - 1, null);
                // Indexes of x and y to be used for fitting
                int xyix[];
                if (s < nlines + 1) {
                    xyix = new int[s];
                    for (int i = 0; i < s; i++) {
                        xyix[i] = i;
                    }
                } else {
                    xyix = new int[nlines];
                    for (int i = s - nlines; i < s; i++) {
                        xyix[i - s + nlines] = i;
                    }
                }
                // Intensity centers + xy
                int icl;
                if (ic != null) {
                    icl = ic.length;
                } else {
                    icl = 0;
                }
                double xw[] = new double[icl + xyix.length];
                double yw[] = new double[icl + xyix.length];
                double xmax = -Double.MAX_VALUE, xmin = -xmax, ymax = xmax, ymin = xmin;
                for (int i = 0; i < xw.length; i++) {
                    if (i < icl) {
                        xw[i] = ic[i][0];
                        yw[i] = ic[i][1];
                    } else {
                        xw[i] = ((Double) x.get(xyix[i - icl])).doubleValue();
                        yw[i] = ((Double) y.get(xyix[i - icl])).doubleValue();
                    }
                    if (xw[i] > xmax) {
                        xmax = xw[i];
                    } else if (xw[i] < xmin) {
                        xmin = xw[i];
                    }
                    if (yw[i] > ymax) {
                        ymax = yw[i];
                    } else if (yw[i] < ymin) {
                        ymin = yw[i];
                    }
                }
                // Fit polyline
                double xfit[], yfit[];
                double p[];
                CurveFitter cf = new CurveFitter(xw, yw);
                cf.doFit(CurveFitter.POLY2);
                p = cf.getParams();
                xfit = new double[ne];
                yfit = new double[ne];
                for (int i = 0; i < ne; i++) {
                    xfit[i] = xmin + i * (xmax - xmin) / ne;
                    yfit[i] = cf.f(p, xfit[i]);
                }
                // Intersection of polyline and lines
                Double array[][];
                for (int i = 1; i < 2 * nlines; i++) {
                    double point[] = intersections2(p[2], p[1], p[0], xl[i][2 * npoints],
                            xl[i][0], yl[i][2 * npoints], yl[i][0]);
                    if (point[0] >= xmin && point[0] <= xmax && point[1] >= ymin && point[1] <= ymax) {
                        Double hcp[][] = (Double[][]) hc.get(i - 1);
                        if (hcp != null) {
                            array = new Double[hcp.length + 1][2];
                            for (int j = 0; j < hcp.length; j++) {
                                array[j][0] = hcp[j][0];
                                array[j][1] = hcp[j][1];
                            }
                            array[hcp.length][0] = new Double(point[0]);
                            array[hcp.length][1] = new Double(point[1]);
                        } else {
                            array = new Double[1][2];
                            array[0][0] = new Double(point[0]);
                            array[0][1] = new Double(point[1]);
                        }
                        hc.set(i - 1, array);
                    }
                }
                // All possible positions 
                Double xyhc[][] = (Double[][]) hc.get(lpix);
                // STOP
                if (xyhc != null) {
                    // Add new x,y
                    double xsum = 0.0, ysum = 0.0;
                    for (int i = 0; i < xyhc.length; i++) {
                        xsum += xyhc[i][0].doubleValue();
                        ysum += xyhc[i][1].doubleValue();
                    }
                    x.add(new Double(xsum / xyhc.length));
                    y.add(new Double(ysum / xyhc.length));
                    tailIntens.add(new Double(im.getInterpolatedValue(xsum / xyhc.length, ysum / xyhc.length)));
                } else {
                    stop = 1;
                }
            }
        }
        return true;
    }

    double[] findStartingVector() {
        double maxIntens = -Double.MAX_VALUE;
        double point[] = {-1.0, -1.0};
        for (int i = -3; i <= 3; i++) {
            for (int j = -3; j <= 3; j++) {
                double current = im.getInterpolatedValue(x1 + i, y1 + j);
                if (current > maxIntens && current > tf) {
                    maxIntens = current;
                    point[0] = x1 + i;
                    point[1] = y1 + j;
                }
            }
        }
        if (point[0] > 0.0 && point[1] > 0.0) {
            return point;
        } else {
            return null;
        }
    }

    public double[] getX() {
        int size = x.size();
        double array[] = new double[size];
        for (int i = 0; i < size; i++) {
            array[i] = ((Double) x.get(i)).doubleValue();
        }
        return array;
    }

    public double[] getY() {
        int size = y.size();
        double array[] = new double[size];
        for (int i = 0; i < size; i++) {
            array[i] = ((Double) y.get(i)).doubleValue();
        }
        return array;
    }

    public double[] getTailIntens() {
        int size = tailIntens.size();
        double array[] = new double[size];
        for (int i = 0; i < size; i++) {
            array[i] = ((Double) tailIntens.get(i)).doubleValue();
        }
        return array;
    }
}
