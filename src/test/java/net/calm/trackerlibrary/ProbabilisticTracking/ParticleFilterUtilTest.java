package net.calm.trackerlibrary.ProbabilisticTracking;

import java.util.Vector;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;

class ParticleFilterUtilTest {

    @Test
    void copyStateVectorDeepCopies() {
        Vector<float[]> orig = new Vector<float[]>();
        orig.add(new float[]{1f, 2f, 3f});
        Vector<float[]> copy = ParticleFilterUtil.copyStateVector(orig);
        assertNotSame(orig.get(0), copy.get(0));
        assertEquals(1f, copy.get(0)[0]);
        assertEquals(2f, copy.get(0)[1]);
        assertEquals(3f, copy.get(0)[2]);
    }

    @Test
    void copyParticleVectorDeepCopiesNestedArrays() {
        Vector<Vector<float[]>> orig = new Vector<Vector<float[]>>();
        Vector<float[]> inner = new Vector<float[]>();
        inner.add(new float[]{4f, 5f});
        orig.add(inner);
        Vector<Vector<float[]>> copy = ParticleFilterUtil.copyParticleVector(orig);
        assertNotSame(orig.get(0).get(0), copy.get(0).get(0));
        assertEquals(5f, copy.get(0).get(0)[1]);
    }

    @Test
    void addImageAccumulatesElementwise() {
        float[][] a = new float[][]{{1f, 2f}, {3f, 4f}};
        float[][] b = new float[][]{{10f, 20f}, {30f, 40f}};
        ParticleFilterUtil.addImage(a, b);
        assertEquals(11f, a[0][0]);
        assertEquals(22f, a[0][1]);
        assertEquals(33f, a[1][0]);
        assertEquals(44f, a[1][1]);
    }

    @Test
    void initArrayToValueFillsEveryCell() {
        float[][] a = new float[2][3];
        ParticleFilterUtil.initArrayToValue(a, 7f);
        for (float[] row : a) {
            for (float v : row) {
                assertEquals(7f, v);
            }
        }
    }
}
