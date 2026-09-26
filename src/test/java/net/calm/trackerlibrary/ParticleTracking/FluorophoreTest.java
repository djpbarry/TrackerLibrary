package net.calm.trackerlibrary.ParticleTracking;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class FluorophoreTest {

    private static final double TOL = 1e-12;

    @Test
    void constructorInitialisesPositionMagnitudeAndThreshold() {
        Fluorophore f = new Fluorophore(3.0, -1.5, 100.0, 0.05);
        assertEquals(3.0, f.getX(), TOL);
        assertEquals(-1.5, f.getY(), TOL);
        assertEquals(100.0, f.getInitialMag(), TOL);
        assertEquals(100.0, f.getCurrentMag(), TOL);
    }

    @Test
    void updateMagWithExplicitValueOverridesCurrentMagnitude() {
        Fluorophore f = new Fluorophore(0.0, 0.0, 100.0, 0.05);
        f.updateMag(42.0);
        assertEquals(42.0, f.getCurrentMag(), TOL);
    }

    @Test
    void subclassDecayRateIsExposed() {
        DecayingFluorophore f = new DecayingFluorophore(0.0, 0.0, 100.0, 0.1);
        assertEquals(0.1, f.getDecayRate(), TOL);
    }

    @Test
    void subclassOnProbabilityIsExposed() {
        BlinkingFluorophore f = new BlinkingFluorophore(0.0, 0.0, 100.0, 0.7);
        assertEquals(0.7, f.getOnProb(), TOL);
    }
}
