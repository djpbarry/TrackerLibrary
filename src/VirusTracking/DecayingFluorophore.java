package VirusTracking;

import java.util.Random;

public class DecayingFluorophore extends Fluorophore {

    private double decayRate;
    private double noise = 0.01;
    private Random rand = new Random();

    public DecayingFluorophore(double x, double y, double mag, double decayRate) {
        super(x, y, mag, 0.05);
        this.decayRate = decayRate;
    }

    public void updateMag() {
        currentMag = currentMag * (1.0 - decayRate + noise * rand.nextGaussian());
    }

    public void updateXMag(int t) {
        if (x > t) {
            currentMag = 255.0 - 2.0 * (x - t) + noise * rand.nextGaussian();
        } else {
            currentMag = 255.0 - 0.5 * (x - t) + noise * rand.nextGaussian();
        }
        if (currentMag < 0.0) {
            currentMag = 0.0;
        }
    }

    public double getDecayRate() {
        return decayRate;
    }
}
