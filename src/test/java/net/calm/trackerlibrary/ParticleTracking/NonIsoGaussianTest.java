package net.calm.trackerlibrary.ParticleTracking;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class NonIsoGaussianTest {

    private static final double TOL = 1e-12;

    @Test
    void evaluateReturnsMagnitudeAtCentreForAxisAlignedGaussian() {
        NonIsoGaussian g = new NonIsoGaussian(1.0, 2.0, 5.0, 0.5, 0.5, 0.0, 0.9);
        assertEquals(5.0, g.evaluate(1.0, 2.0), TOL);
    }

    @Test
    void evaluateDecaysAlongMajorAxis() {
        NonIsoGaussian g = new NonIsoGaussian(0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.9);
        assertEquals(1.0, g.evaluate(0.0, 0.0), TOL);
        // At one sigma along x (theta=0), value should be exp(-0.5).
        assertEquals(Math.exp(-0.5), g.evaluate(1.0, 0.0), TOL);
    }

    @Test
    void rotatedGaussianPreservesPeakValue() {
        NonIsoGaussian g = new NonIsoGaussian(0.0, 0.0, 3.0, 0.5, 0.5, Math.PI / 4, 0.5);
        assertEquals(3.0, g.evaluate(0.0, 0.0), TOL);
    }

    @Test
    void accessorsExposeConstructionParameters() {
        NonIsoGaussian g = new NonIsoGaussian(10.0, -4.0, 2.5, 1.25, 0.75, 0.3, 0.42);
        assertEquals(10.0, g.getX0(), TOL);
        assertEquals(-4.0, g.getY0(), TOL);
        assertEquals(1.25, g.getxSigma(), TOL);
        assertEquals(0.75, g.getySigma(), TOL);
        assertEquals(0.3, g.getTheta(), TOL);
    }
}
