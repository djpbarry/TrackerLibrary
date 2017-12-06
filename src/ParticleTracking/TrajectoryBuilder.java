/*
 * Copyright (C) 2017 David Barry <david.barry at crick.ac.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package ParticleTracking;

import Particle.Particle;
import Particle.ParticleArray;
import IAClasses.ProgressDialog;
import IAClasses.Region;
import java.util.ArrayList;
import org.apache.commons.math3.linear.ArrayRealVector;

public class TrajectoryBuilder {

    private final static int TRAJ_MAX_STEP = 3;

    public static void updateTrajectories(ParticleArray objects, double timeRes, double minStepTol, double spatialRes, double magNormFactor, ArrayList<ParticleTrajectory> trajectories, boolean morph) {
        double mw, vw, pw;
        if (UserVariables.getMotionModel() == UserVariables.RANDOM) {
            mw = 0.0;
            vw = 0.0;
            pw = 1.0;
        } else {
            mw = 0.0;
            vw = 1.0;
            pw = 1.0;
        }
        if (objects == null) {
            return;
        }
        int depth = objects.getDepth();
        ParticleTrajectory traj = null;
        double maxScore;
        ProgressDialog progress = new ProgressDialog(null, "Building Trajectories...", false, "Trajectory Builder", false);
        progress.setVisible(true);
        for (int m = 0; m < depth; m++) {
            progress.updateProgress(m, depth);
            for (int k = m; (k < depth) && (((k - m)) < TRAJ_MAX_STEP); k++) {
                int size = trajectories.size();
                ArrayList<Particle> detections = objects.getLevel(k);
                for (int j = 0; j < detections.size(); j++) {
                    Particle currentParticle = detections.get(j);
                    if (currentParticle != null) {
                        /*
                         * If no trajectories have yet been built, start a new
                         * one:
                         */
                        if (k == m) {
                            traj = new ParticleTrajectory(timeRes, spatialRes);
                            /*
                             * Particles need to be cloned as they are set to
                             * null once inserted into trajectories.
                             */
                            traj.addPoint(currentParticle.makeCopy());
                            trajectories.add(traj);
                            /*
                             * Otherwise, determine whether the current particle
                             * belongs to a pre-existing trajectory:
                             */
                        } else {
                            int i;
                            int minScoreIndex;
                            for (minScoreIndex = -1, maxScore = -Double.MAX_VALUE, i = 0; i < size; i++) {
                                traj = (ParticleTrajectory) trajectories.get(i);
                                Particle last = traj.getEnd();
                                if ((last != null) && (last.getFrameNumber() == m) && k != m) {
                                    Region currentRegion = currentParticle.getRegion();
                                    Region lastRegion = last.getRegion();
                                    double morphScore = 0.0;
                                    if (morph) {
                                        ArrayRealVector morphvector1 = currentRegion.getMorphMeasures();
                                        ArrayRealVector morphvector2 = lastRegion.getMorphMeasures();
                                        morphScore = 1.0 - morphvector1.getDistance(morphvector2) / morphvector1.getL1Norm();
                                    }
                                    double x = currentParticle.getX();
                                    double y = currentParticle.getY();
                                    ArrayRealVector vector1 = new ArrayRealVector(new double[]{x, y});
                                    ArrayRealVector vector2 = new ArrayRealVector(new double[]{last.getX(), last.getY()});
                                    double posScore = 1.0 - vector1.getDistance(vector2) / vector1.getL1Norm();
                                    double projScore = 1.0;
                                    if (UserVariables.getMotionModel() != UserVariables.RANDOM) {
                                        double deltaT = currentParticle.getFrameNumber() * UserVariables.getTimeRes() - last.getFrameNumber() * UserVariables.getTimeRes();
                                        ArrayRealVector vector3 = new ArrayRealVector(new double[]{x, y});
                                        ArrayRealVector vector4 = new ArrayRealVector(new double[]{last.getX() + traj.getXVelocity() * deltaT, last.getY() + traj.getYVelocity() * deltaT});
                                        projScore = 1.0 - vector3.getDistance(vector4) / vector3.getL1Norm();
                                    }
                                    double totScore = (mw * morphScore + vw * projScore + pw * posScore) / (mw + vw + pw);
                                    if (totScore > maxScore) {
                                        maxScore = totScore;
                                        minScoreIndex = i;
                                    }
                                }
                            }
                            /*
                             * If an acceptably low score has been evaluated,
                             * the particle is temporarily assigned to the
                             * "winning" trajectory:
                             */
                            if (traj != null) {
                                if (minScoreIndex > -1) {
                                    traj = (ParticleTrajectory) trajectories.get(minScoreIndex);
                                }
                                if ((maxScore > minStepTol) && (maxScore > traj.getTempScore())) {
                                    traj.addTempPoint(currentParticle.makeCopy(), maxScore, j, k);
                                }
                            }
                        }
                    }
                }
            }
            for (ParticleTrajectory trajectory : trajectories) {
                traj = (ParticleTrajectory) trajectory;
                Particle temp = traj.getTemp();
                if (temp != null) {
                    int row = traj.getTempRow();
                    int col = traj.getTempColumn();
                    if (col <= m + 1) {
                        traj.checkDetections(temp, 0.0);
                        objects.nullifyDetection(col, row);
                    }
                }
            }
        }
        progress.dispose();
    }

}
