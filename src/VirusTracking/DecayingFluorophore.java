package VirusTracking;

import java.util.Random;

public class DecayingFluorophore {

    private double x;
    private double y;
    private double initialMag;
    private double currentMag;
    private double decayRate;
    private double noise = 0.01;
    private Random rand = new Random();

    public DecayingFluorophore(double x, double y, double mag, double decayRate) {
        this.x = x;
        this.y = y;
        this.initialMag = mag;
        this.currentMag = mag;
        this.decayRate = decayRate;
    }

    public void updateMag() {
        currentMag = currentMag * (1.0 - decayRate + noise * rand.nextGaussian());
    }
    
    public void updateMag(double newMag) {
        currentMag = newMag;
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

    public double getCurrentMag() {
        return currentMag;
    }

    public double getDecayRate() {
        return decayRate;
    }

    public double getInitialMag() {
        return initialMag;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }
}
