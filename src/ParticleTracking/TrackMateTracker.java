/*
 * Copyright (C) 2019 David Barry <david.barry at crick dot ac dot uk>
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

import Particle.Blob;
import Particle.IsoGaussian;
import Particle.Particle;
import Particle.Point;
import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.TrackModel;
import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTracker;
import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory;
import ij.IJ;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import java.util.Set;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

/**
 *
 * @author David Barry <david.barry at crick dot ac dot uk>
 */
public class TrackMateTracker {

    private TrackModel tm;

    public TrackMateTracker() {

    }

    public void track(SpotCollection spots, Map<String, Object> settings) {
            SparseLAPTracker tracker = (SparseLAPTracker) (new SparseLAPTrackerFactory()).create(spots, settings);
            if (!tracker.process()) {
                IJ.log("Tracking failed.");
                return;
            }
        SimpleWeightedGraph<Spot, DefaultWeightedEdge> graph = tracker.getResult();
        Model model = new Model();
        model.setTracks(graph, true);
        tm = model.getTrackModel();
    }

    public void updateTrajectories(ArrayList<ParticleTrajectory> trajectories) {
        Set<Integer> trackIds = tm.trackIDs(false);
        for (Integer id : trackIds) {
            ParticleTrajectory traj = new ParticleTrajectory();
            Set<Spot> spots = tm.trackSpots(id);
            final Comparator< Spot> comparator = Spot.frameComparator;
            final ArrayList< Spot> sorted = new ArrayList< Spot>(spots);
            Collections.sort(sorted, comparator);
            for (Spot s : sorted) {
                Particle p = null;
                double x = s.getFeature(Spot.POSITION_X);
                double y = s.getFeature(Spot.POSITION_Y);
                double t = s.getFeature(Spot.FRAME);
                if (s instanceof Point) {
                    p = new Point((int) t, x, y, 0.0);
                } else if (s instanceof Blob) {
                    p = new Blob((int) t, x, y, 0.0);
                } else if (s instanceof IsoGaussian) {
                    p = new IsoGaussian((int) t, x, y, 0.0, s.getFeature(Spot.RADIUS), s.getFeature(Spot.RADIUS), s.getFeature(Spot.QUALITY), null, 0, null);
                }
                p.putFeature(Spot.FRAME, s.getFeature(Spot.FRAME));
                traj.addPoint(p);
            }
            trajectories.add(traj);
        }
    }
}
