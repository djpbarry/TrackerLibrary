package ParticleTracking;

import IAClasses.IsoGaussian;
import java.util.ArrayList;

/**
 *
 * @author barry05
 */
public class ParticleArray {

    private int depth;
    private ArrayList detections[];

    public ParticleArray(int depth) {
        this.depth = depth;
        if (depth > 0) {
            detections = new ArrayList[depth];
            for (int i = 0; i < depth; i++) {
                detections[i] = new ArrayList();
            }
        } else {
            this.depth = 0;
            detections = null;
        }
    }

    public boolean addDetection(int level, IsoGaussian g1, IsoGaussian g2) {
        if (detections == null || level > depth - 1) {
            return false;
        }
        IsoGaussian detection[] = {g1, g2};
        detections[level].add(detection);
        return true;
    }

    public boolean nullifyDetection(int level, int index) {
        if (detections == null || level > depth - 1) {
            return false;
        }
        IsoGaussian detection[] = {null, null};
        detections[level].set(index, detection);
        return true;
    }

    public int getDepth() {
        return depth;
    }

    public ArrayList getLevel(int level) {
        if (detections != null) {
            return detections[level];
        } else {
            return null;
        }
    }
}
