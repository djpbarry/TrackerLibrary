/*
 * Copyright (C) 2015 David Barry <david.barry at cancer.org.uk>
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

import IAClasses.IsoGaussian;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import java.util.ArrayList;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class TrajectoryBuilder {

    public static void updateTrajectories(ParticleArray objects, double timeRes, double trajMaxStep, double spatialRes, boolean projectPos, double magNormFactor, ArrayList<ParticleTrajectory> trajectories) {
        if (objects == null) {
            return;
        }
        int depth = objects.getDepth();
        ParticleTrajectory traj = null;
        double x;
        double y;
        double score;
        double minScore;
        ProgressDialog progress = new ProgressDialog(null, "Building Trajectories...", false, "Trajectory Builder", false);
        progress.setVisible(true);
        for (int m = 0; m < depth; m++) {
            progress.updateProgress(m, depth);
            for (int k = m; (k < depth) && (((k - m)) < trajMaxStep); k++) {
                int size = trajectories.size();
                ArrayList<Particle> detections = objects.getLevel(k);
                for (int j = 0; j < detections.size(); j++) {
                    Particle currentParticle = detections.get(j);
                    if (currentParticle != null) {
                        IsoGaussian ch1G = currentParticle.getC1Gaussian();
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
                            traj.addPoint((Particle) currentParticle.clone());
                            trajectories.add(traj);
                            /*
                             * Otherwise, determine whether the current particle
                             * belongs to a pre-existing trajectory:
                             */
                        } else {
                            int i;
                            int minScoreIndex;
                            for (minScoreIndex = -1, minScore = Double.MAX_VALUE, i = 0; i < size; i++) {
                                traj = (ParticleTrajectory) trajectories.get(i);
                                Particle last = traj.getEnd();
                                if ((last != null) && (last.getTimePoint() == m) && k != m) {
                                    /*
                                     * Evaluate the probability that the current
                                     * particle belongs to the current
                                     * trajectory, based on the particle's
                                     * distance from last point on the current
                                     * trajectory and the number of frames
                                     * between the current particle, the last
                                     * point of the current trajectory and
                                     * differences in respective intensity
                                     * levels:
                                     */
                                    x = ch1G.getX();
                                    y = ch1G.getY();
                                    double[] vector1 = {x, y, currentParticle.getTimePoint(), ch1G.getMagnitude() / magNormFactor};
                                    double[] vector2 = {last.getX(), last.getY(), last.getTimePoint(), last.getC1Intensity() / magNormFactor};
                                    score = Utils.calcEuclidDist(vector1, vector2);
                                    if (projectPos) {
                                        double[] vector3 = {x, y};
                                        double[] vector4 = {last.getX() + traj.getXVelocity(), last.getY() + traj.getYVelocity()};
                                        score += Utils.calcEuclidDist(vector3, vector4);
                                    }
                                    if (score < minScore) {
                                        minScore = score;
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
                                if ((minScore < trajMaxStep) && (minScore < traj.getTempScore())) {
                                    traj.addTempPoint((Particle) currentParticle.clone(), minScore, j, k);
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
