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
    void intersections2ReturnsPositiveRootWhenInsideBoundingBox() {
        // Parabola y = x^2 (a=1, b=0, c2=0) meets line y = x through (0.5,0.5)-(2,2).
        // Roots x = 0 and x = 1; only (1,1) lies inside the box, so the +root is returned.
        double[] p = TailTracer.intersections2(1.0, 0.0, 0.0, 0.5, 2.0, 0.5, 2.0);
        assertEquals(1.0, p[0], TOL);
        assertEquals(1.0, p[1], TOL);
    }

    @Test
    void intersections2ReturnsNegativeRootWhenPositiveRootOutsideBox() {
        // Same parabola y = x^2 and line y = x, but the box through (0,0)-(0.5,0.5)
        // contains only the origin; the +root (1,1) is outside, so the -root is returned.
        double[] p = TailTracer.intersections2(1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5);
        assertEquals(0.0, p[0], TOL);
        assertEquals(0.0, p[1], TOL);
    }
}
