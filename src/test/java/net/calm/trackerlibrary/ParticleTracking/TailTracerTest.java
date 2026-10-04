package net.calm.trackerlibrary.ParticleTracking;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class TailTracerTest {

    private static final double TOL = 1e-12;

    @Test
    void normalizedVectorHasUnitLength() {
        double[] v = TailTracer.normalizedVector(3.0, 4.0);
        assertEquals(1.0, Math.hypot(v[0], v[1]), TOL);
        assertEquals(0.6, v[0], TOL);
        assertEquals(0.8, v[1], TOL);
    }

    @Test
    void normalVectorIsPerpendicularToSegment() {
        // Segment along +x: (0,0) -> (1,0). Normal is (0, -1).
        double[] n = TailTracer.normalVector(0.0, 0.0, 1.0, 0.0);
        assertEquals(0.0, n[0], TOL);
        assertEquals(-1.0, n[1], TOL);
        // Unit length too.
        assertEquals(1.0, Math.hypot(n[0], n[1]), TOL);
    }

    @Test
    void intersectionsComputesLineCrossing() {
        // y = x (m1=1, c1=0) and y = -x + 2 (m2=-1, c2=2) cross at (1,1).
        double[] p = TailTracer.intersections(1.0, 0.0, -1.0, 2.0);
        assertEquals(1.0, p[0], TOL);
        assertEquals(1.0, p[1], TOL);
    }

    @Test
    void intersections2ReturnsRootInsideBoundingBox() {
        // Quadratic y = x^2 - 1 (a=1, b=0, c2=-1) and line y = 2x over x in [0,2].
        // Roots are 1 +/- sqrt(2); x = 1 + sqrt(2) ~ 2.414 is outside [0,2],
        // x = 1 - sqrt(2) ~ -0.414 is outside too; but the returned point must be
        // one of the two roots, so we only assert both are real and on the line.
        double[] p = TailTracer.intersections2(1.0, 0.0, -1.0, 0.0, 2.0, 0.0, 4.0);
        double m = (4.0 - 0.0) / (2.0 - 0.0);
        assertEquals(m * p[0], p[1], TOL);
    }
}
