package ParticleTracking;

import IAClasses.IsoGaussian;

/**
 * Represents a detected particle in an individual image or frame.
 *
 * @author David J Barry
 * @version 2.0, FEB 2011
 */
public class Particle {

    private int iD;
    private double t;
    private double x, y;
    private Particle link;
    private IsoGaussian c1Gaussian, c2Gaussian;

    public Particle() {
        t = 0.0;
        x = 0.0;
        y = 0.0;
        iD = -1;
        link = null;
        c1Gaussian = null;
        c2Gaussian = null;
    }

    /**
     * Creates a new Particle
     *
     * @param time position within an image stack
     * @param c1Gaussian IsoGaussian representation of particle in red channel
     * @param c2Gaussian IsoGaussian representation of particle in green channel
     * @param newLink the last
     * <code>Particle</code> in the current
     * <code>ParticleTrajectory</code> to which this particle should be linked.
     * Set to
     * <code>null</code> if this is the first particle in a new trajectory.
     */
    public Particle(double time, IsoGaussian c1Gaussian, IsoGaussian c2Gaussian,
            Particle newLink, int iD) {
        this.iD = iD;
        this.t = time;
        this.c1Gaussian = c1Gaussian;
        this.c2Gaussian = c2Gaussian;
        if (c1Gaussian != null) {
            this.x = c1Gaussian.getX();
            this.y = c1Gaussian.getY();
        } else {
            x = y = Double.NaN;
        }
        this.link = newLink;
    }

    /**
     * Returns the particle's mean x-position, calculated as the average of the
     * x-positions in the red and green channels.
     */
    public double getX() {
        return x;
    }

    /**
     * Returns the particle's mean y-position, calculated as the average of the
     * y-positions in the red and green channels.
     */
    public double getY() {
        return y;
    }

    /**
     * This particle's z-position within an image stack.
     */
    public double getTimePoint() {
        return t;
    }

    /**
     * Returns the peak intensity of this particle in the red channel.
     */
    public double getC1Intensity() {
        if (c1Gaussian != null) {
            return c1Gaussian.getMagnitude();
        } else {
            return 0.0;
        }
    }

    /**
     * Returns the peak intensity of this particle in the green channel.
     */
    public double getC2Intensity() {
        if (c2Gaussian != null) {
            return c2Gaussian.getMagnitude();
        } else {
            return 0.0;
        }
    }

    /**
     * Returns the
     * <code>Particle</code> to which this particle is linked.
     */
    public Particle getLink() {
        return link;
    }

    public IsoGaussian getC2Gaussian() {
        return c2Gaussian;
    }

    public IsoGaussian getC1Gaussian() {
        return c1Gaussian;
    }

    protected Object clone() {
        return new Particle(t, (IsoGaussian) c1Gaussian.clone(), (IsoGaussian) c2Gaussian.clone(), (Particle) link.clone(), -1);
    }

    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Particle other = (Particle) obj;
        if (Double.doubleToLongBits(this.t) != Double.doubleToLongBits(other.t)) {
            return false;
        }
        if (Double.doubleToLongBits(this.x) != Double.doubleToLongBits(other.x)) {
            return false;
        }
        if (Double.doubleToLongBits(this.y) != Double.doubleToLongBits(other.y)) {
            return false;
        }
        if (this.link != other.link && (this.link == null || !this.link.equals(other.link))) {
            return false;
        }
        if (this.c1Gaussian != other.c1Gaussian && (this.c1Gaussian == null || !this.c1Gaussian.equals(other.c1Gaussian))) {
            return false;
        }
        if (this.c2Gaussian != other.c2Gaussian && (this.c2Gaussian == null || !this.c2Gaussian.equals(other.c2Gaussian))) {
            return false;
        }
        return true;
    }
}
