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
package net.calm.trackerlibrary.ParticleTracking;

import java.util.ArrayList;
import java.util.Arrays;

import net.calm.iaclasslibrary.IAClasses.ProgressDialog;
import net.calm.iaclasslibrary.IAClasses.Region;
import net.calm.iaclasslibrary.Particle.Particle;
import net.calm.iaclasslibrary.Particle.ParticleArray;
import org.apache.commons.math3.linear.ArrayRealVector;

public class TrajectoryBuilder {

    /**
     * Links detections across successive frames into {@link ParticleTrajectory}s
     * using a greedy nearest-neighbour score of position, projected velocity,
     * and morphology. Each detection that cannot be linked starts a new
     * single-point trajectory.
     *
     * @param objects the detections organised by frame.
     * @param timeRes the time resolution (seconds per frame).
     * @param minStepTol the minimum step tolerance (reserved; see
     * {@code UserVariables#getTrajMaxStep()} for the active threshold).
     * @param spatialRes the spatial resolution (pixels per unit).
     * @param magNormFactor the magnitude normalisation factor.
     * @param trajectories the output list, appended to in place.
     * @param morph whether to include a morphology term in the score.
     */
    public static void updateTrajectories(ParticleArray objects, double timeRes, double minStepTol, double spatialRes, double magNormFactor, ArrayList<ParticleTrajectory> trajectories, boolean morph) {
        double mw, vw, pw;
        if (UserVariables.getInstance().getMotionModel() == UserVariables.RANDOM) {
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
        ProgressDialog progress = new ProgressDialog(null, "Building Trajectories...", false, "net.calm.trackerlibrary.Trajectory Builder", false);
        progress.setVisible(true);
        for (int m = 0; m < depth; m++) {
            progress.updateProgress(m, depth);
            ArrayList<Particle> detections = objects.getLevel(m);
            for (Particle currentParticle : detections) {
                if (currentParticle != null) {
                    ParticleTrajectory traj = new ParticleTrajectory(timeRes, spatialRes);
                    traj.addPoint(currentParticle.makeCopy());
                    trajectories.add(traj);
                }
            }
            if (m >= depth - 1) {
                continue;
            }
            int k = m + 1;
            int tSize = trajectories.size();
            detections = objects.getLevel(k);
            int dSize = detections.size();
            int[] terminatedTrajMap = getTerminatedTrajMap(trajectories, m);
            int ttSize = terminatedTrajMap.length;
            double[][] scores = new double[ttSize][dSize];
            for (int t = 0; t < ttSize; t++) {
                Arrays.fill(scores[t], Double.MAX_VALUE);
            }
            for (int j = 0; j < dSize; j++) {
                Particle currentParticle = detections.get(j);
                if (currentParticle != null) {
                    double minScore = Double.MAX_VALUE;
                    int minIndex = -1;
                    for (int i = 0; i < tSize; i++) {
                        ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(i);
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
                            double posScore = vector1.getDistance(vector2);
                            double projScore = 1.0;
                            if (UserVariables.getInstance().getMotionModel() != UserVariables.RANDOM) {
                                double deltaT = currentParticle.getFrameNumber() * UserVariables.getInstance().getTimeRes() - last.getFrameNumber() * UserVariables.getInstance().getTimeRes();
                                ArrayRealVector vector3 = new ArrayRealVector(new double[]{x, y});
                                ArrayRealVector vector4 = new ArrayRealVector(new double[]{last.getX() + traj.getXVelocity() * deltaT, last.getY() + traj.getYVelocity() * deltaT});
                                projScore = 1.0 - vector3.getDistance(vector4) / vector3.getL1Norm();
                            }
                            double score = (mw * morphScore + vw * projScore + pw * posScore) / (mw + vw + pw);
                            if (score < minScore) {
                                minScore = score;
                                minIndex = i;
                            }
                        }
                    }
                    if (minIndex > -1) {
                        ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(minIndex);
                        if ((minScore < UserVariables.getInstance().getTrajMaxStep()) && (minScore < traj.getTempScore())) {
                            traj.addTempPoint(currentParticle.makeCopy(), minScore, j, k);
                        }
                    }
                }
            }
            for (ParticleTrajectory trajectory : trajectories) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectory;
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

    private static int[] getTerminatedTrajMap(ArrayList<ParticleTrajectory> trajectories, int frame) {
        int count = 0;
        for (ParticleTrajectory traj : trajectories) {
            if (traj.getEnd().getFrameNumber() == frame) {
                count++;
            }
        }
        int[] result = new int[count];
        count = 0;
        for (int i = 0; i < trajectories.size(); i++) {
            ParticleTrajectory traj = trajectories.get(i);
            if (traj.getEnd().getFrameNumber() == frame) {
                result[count++] = i;
            }
        }
        return result;
    }
}
