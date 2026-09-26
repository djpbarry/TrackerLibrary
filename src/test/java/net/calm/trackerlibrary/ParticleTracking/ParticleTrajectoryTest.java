package net.calm.trackerlibrary.ParticleTracking;

import net.calm.iaclasslibrary.Particle.Particle;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

class ParticleTrajectoryTest {

    private static final double TOL = 1e-9;

    private static Particle particle(int frame, double x, double y, double mag) {
        return new Particle(frame, x, y, mag);
    }

    private static ParticleTrajectory trajectory() {
        return new ParticleTrajectory(1.0, 1.0);
    }

    @Test
    void emptyTrajectoryHasNoStart() {
        assertNull(trajectory().getStart());
        assertEquals(0, trajectory().getSize());
        assertEquals(0, trajectory().getNumberOfFrames());
    }

    @Test
    void addPointAppendsAndLinksToStart() {
        ParticleTrajectory traj = trajectory();
        traj.addPoint(particle(0, 0.0, 0.0, 1.0));
        traj.addPoint(particle(1, 3.0, 4.0, 2.0));
        traj.addPoint(particle(2, 6.0, 8.0, 3.0));

        assertEquals(3, traj.getSize());
        assertEquals(particle(2, 6.0, 8.0, 3.0).getFrameNumber(), traj.getEnd().getFrameNumber());
        assertEquals(0, traj.getStart().getFrameNumber());
    }

    @Test
    void getNumberOfFramesReturnsFrameSpan() {
        ParticleTrajectory traj = trajectory();
        traj.addPoint(particle(5, 0.0, 0.0, 1.0));
        traj.addPoint(particle(9, 1.0, 1.0, 1.0));
        assertEquals(4, traj.getNumberOfFrames());
    }

    @Test
    void getDisplacementSumsEuclideanDistances() {
        ParticleTrajectory traj = trajectory();
        // 3-4-5 triangle: from the newest particle, one link step is 5.0.
        traj.addPoint(particle(0, 0.0, 0.0, 1.0));
        traj.addPoint(particle(1, 3.0, 4.0, 1.0));
        assertNull(traj.getStart().getLink());
        // The chain flows newest -> oldest, so walk from the end.
        assertEquals(5.0, traj.getDisplacement(traj.getEnd(), 1), TOL);
    }

    @Test
    void getTypeClassifiesColocalisationByThreshold() {
        ParticleTrajectory traj = trajectory();
        Particle a = particle(0, 0.0, 0.0, 1.0);
        Particle b = particle(1, 1.0, 1.0, 1.0);
        b.setColocalisedParticle(particle(1, 1.0, 1.0, 1.0));
        traj.addPoint(a);
        traj.addPoint(b);

        // calcDualScore counts linked particles with a colocalised partner (b).
        // With threshold < 0.5 the trajectory (1 of 2 colocal) is NOT colocal;
        // with threshold > 0.5 it is.
        assertEquals(ParticleTrajectory.NON_COLOCAL, traj.getType(0.6));
        assertEquals(ParticleTrajectory.COLOCAL, traj.getType(0.4));
    }

    @Test
    void addTrajectorySplicesPointsInOrder() {
        ParticleTrajectory a = trajectory();
        a.addPoint(particle(0, 0.0, 0.0, 1.0));

        ParticleTrajectory b = trajectory();
        b.addPoint(particle(1, 1.0, 1.0, 1.0));
        b.addPoint(particle(2, 2.0, 2.0, 1.0));

        a.addTrajectory(b);
        assertEquals(3, a.getSize());
        assertEquals(0, a.getStart().getFrameNumber());
        assertEquals(2, a.getEnd().getFrameNumber());
    }
}
