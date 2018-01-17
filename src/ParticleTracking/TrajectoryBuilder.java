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
import java.util.Arrays;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class TrajectoryBuilder {

//    private final static int TRAJ_MAX_STEP = 3;
//    private final static double END_TRAJ = 0.5;

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
        ProgressDialog progress = new ProgressDialog(null, "Building Trajectories...", false, "Trajectory Builder", false);
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
//                    if (m >= 190) {
//                        System.out.println(String.format("particle: %d X: %f Y: %f", j, currentParticle.getX(), currentParticle.getY()));
//                    }
                    double minScore = Double.MAX_VALUE;
                    int minIndex = -1;
                    for (int i = 0; i < tSize; i++) {
//                        System.out.println(String.format("m: %d k: %d i: %d j: %d tSize: %d dSize: %d", m, k, i, j, tSize, dSize));
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
                            if (UserVariables.getMotionModel() != UserVariables.RANDOM) {
                                double deltaT = currentParticle.getFrameNumber() * UserVariables.getTimeRes() - last.getFrameNumber() * UserVariables.getTimeRes();
                                ArrayRealVector vector3 = new ArrayRealVector(new double[]{x, y});
                                ArrayRealVector vector4 = new ArrayRealVector(new double[]{last.getX() + traj.getXVelocity() * deltaT, last.getY() + traj.getYVelocity() * deltaT});
                                projScore = 1.0 - vector3.getDistance(vector4) / vector3.getL1Norm();
                            }
                            double score = (mw * morphScore + vw * projScore + pw * posScore) / (mw + vw + pw);
                            if (score < minScore) {
                                minScore = score;
                                minIndex = i;
                            }
//                            if (m >= 190) {
//                                System.out.println(String.format("traj: %d X: %f Y: %f score: %f", i, last.getX(), last.getY(), scores[s - 1][j]));
//                            }
                        }
//                        if (minIndex > -1 && minScore < END_TRAJ) {
////                    if (minScores[t] < 1.0 - minStepTol) {
////                        Particle currentParticle = objects.getLevel(k).get(minIndices[t]);
//                            traj.addTempPoint(currentParticle.makeCopy(), minScore, j, k);
//                            System.out.println(String.format("t: %d traj: %d particleX: %f particleY: %f score: %f", m, i, currentParticle.getX(), currentParticle.getY(), minScore));
////                    }
//                        }
                    }
                    if (minIndex > -1) {
                        ParticleTrajectory traj = (ParticleTrajectory) trajectories.get(minIndex);
                        if ((minScore < UserVariables.getTrajMaxStep()) && (minScore < traj.getTempScore())) {
                            traj.addTempPoint(currentParticle.makeCopy(), minScore, j, k);
                        }
                    }
                }
            }
//            if (scores.length > 0) {
//                int[] minIndices = getMinScoreIndices(scores);
////                double[] minScores = getMinScores(scores, minIndices);
//                for (int t = 0; t < scores.length; t++) {
//                    ParticleTrajectory traj = trajectories.get(terminatedTrajMap[t]);
//                    if (minIndices[t] < dSize) {
////                    if (minScores[t] < 1.0 - minStepTol) {
//                        Particle currentParticle = objects.getLevel(k).get(minIndices[t]);
//                        traj.addTempPoint(currentParticle.makeCopy(), scores[t][minIndices[t]], minIndices[t], k);
//                        System.out.println(String.format("t: %d traj: %d particleX: %f particleY: %f score: %f", m, t, currentParticle.getX(), currentParticle.getY(), scores[t][minIndices[t]]));
////                    }
//                    }
//                }
//            }
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

    private static double[] getMinScores(double[][] scores, int[] indices) {
        DescriptiveStatistics ds = new DescriptiveStatistics();
        double[] result = new double[scores.length];
        for (int i = 0; i < scores.length; i++) {
            if (indices[i] >= scores[i].length) {
                ds.addValue(UserVariables.getTrajMaxStep());
            } else {
                ds.addValue(scores[i][indices[i]]);
            }
        }
        ds.addValue(UserVariables.getTrajMaxStep()); //ensures sum is non-zero
        double sum = ds.getSum();
        for (int i = 0; i < scores.length; i++) {
            if (indices[i] >= scores[i].length) {
                result[i] = UserVariables.getTrajMaxStep() / sum;
            } else {
                result[i] = scores[i][indices[i]] / sum;
            }
        }
        return result;
    }

    private static int[] getMinScoreIndices(double[][] scores) {
        int nT = scores.length;
        int nD = scores[0].length;
        int[] result = new int[nT];
        int[] currentCombo = new int[nT];
        Arrays.fill(currentCombo, 0);
        double minResult = Double.MAX_VALUE;
//        int N = getNCombs(nD, nT);
//        int[][] allCombs = calcAllPossibleCombs(nT, nD);
//        int winningComb = -1;
//        for (int i = 0; i < N; i++) {
        while (currentCombo != null) {
//            if (allCombs[i] != null) {
            if (allUnique(currentCombo, nD)) {
                double currentScore = calcScore(scores, currentCombo);
//                System.out.println(String.format("currentScore: %f currentCombo: %d %d %d currentScore: %f minScore: %f", currentScore, currentCombo[0], currentCombo[1], currentCombo[2], currentScore, minResult));
                if (currentScore < minResult) {
                    minResult = currentScore;
                    System.arraycopy(currentCombo, 0, result, 0, currentCombo.length);
                }
            }
//            }
            currentCombo = increment(currentCombo, currentCombo.length - 1, nD);
        }
        return result;
    }

    private static double calcScore(double[][] scores, int[] indices) {
        double score = 0.0;
        for (int i = 0; i < scores.length; i++) {
            if (indices[i] == scores[i].length) {
                score += UserVariables.getTrajMaxStep();
            } else {
                score += scores[i][indices[i]];
            }
        }
        return score;
    }

    private static int getNCombs(int nD, int nT) {
        int result = 1;
        for (int d = 0; d < nT; d++) {
            result *= nD + 1;
        }
        return result;
    }

//    private static int[][] calcAllPossibleCombs(int nT, int nD) {
//        int unassigned = nD;
//        int N = getNCombs(nD, nT);
//        ArrayList<int[]> result = new ArrayList();
//        int n = 1;
//        for (int[] r : result) {
//            Arrays.fill(r, 0);
//        }
//        int[] last = new int[nT];
//        System.arraycopy(result[0], 0, last, 0, last.length);
//        while (n < N) {
//            int[] current = new int[nT];
//            System.arraycopy(last, 0, current, 0, current.length);
//            current = increment(current, current.length - 1, unassigned);
////            System.out.print(String.format("n:%d %d %d %d ", n, current[0], current[1], current[2]));
//            if (allUnique(current, unassigned)) {
////                System.out.println(true);
//                System.arraycopy(current, 0, result[n], 0, current.length);
//            } else {
////                System.out.println(false);
//                result[n] = null;
//            }
//            System.arraycopy(current, 0, last, 0, last.length);
//            n++;
//        }
//        return result;
//    }
    private static int[] increment(int[] indices, int index, int unassigned) {
        indices[index]++;
        if (indices[index] > unassigned) {
            indices[index] = 0;
            if (index > 0) {
                return increment(indices, index - 1, unassigned);
            } else {
                return null;
            }
        }
        return indices;
    }

    private static int[] getFirstResult(int nT, int nD) {
        int[] result = new int[nT];
        Arrays.fill(result, nD);
        for (int n = 0; n < nT && n < nD; n++) {
            result[n] = n;
        }
        return result;
    }

    private static boolean allUnique(int[] indices, int unassigned) {
        for (int i = 0; i < indices.length; i++) {
            for (int j = i + 1; j < indices.length; j++) {
                if (indices[i] != unassigned && indices[i] == indices[j]) {
                    return false;
                }
            }
        }
        return true;
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
