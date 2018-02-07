/*
 * Copyright (C) 2018 David Barry <david dot barry at crick dot ac dot uk>
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
package Trajectory;

import IAClasses.ProgressDialog;
import IAClasses.Region;
import Particle.Particle;
import ParticleTracking.ParticleTrajectory;
import ParticleTracking.UserVariables;
import java.util.ArrayList;
import org.apache.commons.math3.linear.ArrayRealVector;

public class TrajectoryBridger {

    public static void bridgeTrajectories(ArrayList<ParticleTrajectory> trajectories, double[] scoreWeightings, int maxStep) {
        int size = trajectories.size();
        ProgressDialog progress = new ProgressDialog(null, "Processing Trajectories...", false, "Trajectory Builder", false);
        progress.setVisible(true);
        for (int m = 0; m < size; m++) {
            progress.updateProgress(m, size);
            ParticleTrajectory traj1 = trajectories.get(m);
            Particle traj1End = traj1.getEnd();
            if (traj1End == null) {
                continue;
            }
            double minScore = Double.MAX_VALUE;
            int minIndex = -1;
            for (int n = m + 1; n < size; n++) {
                ParticleTrajectory traj2 = trajectories.get(n);
                Particle traj2Start = traj2.getStart();
                if (traj2Start != null) {
                    int stepSize = traj2Start.getFrameNumber() - traj1End.getFrameNumber();
                    if (stepSize > 0 && stepSize <= maxStep) {
                        double morphScore = 0.0;
                        Region currentRegion = traj2Start.getRegion();
                        Region lastRegion = traj1End.getRegion();
                        if (currentRegion != null && lastRegion != null) {
                            ArrayRealVector morphvector1 = currentRegion.getMorphMeasures();
                            ArrayRealVector morphvector2 = lastRegion.getMorphMeasures();
                            morphScore = 1.0 - morphvector1.getDistance(morphvector2) / morphvector1.getL1Norm();
                        }
                        double x = traj2Start.getX();
                        double y = traj2Start.getY();
                        ArrayRealVector vector1 = new ArrayRealVector(new double[]{x, y});
                        ArrayRealVector vector2 = new ArrayRealVector(new double[]{traj1End.getX(), traj1End.getY()});
                        double posScore = vector1.getDistance(vector2);
                        double projScore = 1.0;
                        if (UserVariables.getMotionModel() != UserVariables.RANDOM) {
                            double deltaT = traj2Start.getFrameNumber() * UserVariables.getTimeRes() - traj1End.getFrameNumber() * UserVariables.getTimeRes();
                            ArrayRealVector vector3 = new ArrayRealVector(new double[]{x, y});
                            ArrayRealVector vector4 = new ArrayRealVector(new double[]{traj1End.getX() + traj1.getXVelocity() * deltaT, traj1End.getY() + traj1.getYVelocity() * deltaT});
                            projScore = 1.0 - vector3.getDistance(vector4) / vector3.getL1Norm();
                        }
                        double score = (scoreWeightings[0] * morphScore + scoreWeightings[1] * projScore + scoreWeightings[2] * posScore) / (scoreWeightings[0] + scoreWeightings[1] + scoreWeightings[2]);
                        if (score < minScore) {
                            minScore = score;
                            minIndex = n;
                        }
                    }
                }
            }
            if (minIndex > -1) {
                ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(minIndex);
                if (minScore < UserVariables.getTrajMaxStep()) {
                    traj1.addTrajectory(traj);
                    trajectories.remove(minIndex);
                    size--;
                    m--;
                }
            }
        }
        progress.dispose();
    }
}
