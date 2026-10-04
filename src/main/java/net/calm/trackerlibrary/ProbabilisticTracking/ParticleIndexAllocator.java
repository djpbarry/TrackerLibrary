package net.calm.trackerlibrary.ProbabilisticTracking;

import java.util.concurrent.atomic.AtomicInteger;

/**
 * Thread-safe allocator that hands out each index in {@code [0, n)} exactly
 * once across all callers, then returns {@code -1} forever.
 */
public class ParticleIndexAllocator {

    private final AtomicInteger nextIndex;

    public ParticleIndexAllocator() {
        this.nextIndex = new AtomicInteger(0);
    }

    /**
     * @param n the (exclusive) upper bound on allocated indices
     * @return the next index in {@code [0, n)}, or {@code -1} once exhausted
     */
    public int next(int n) {
        int index = nextIndex.getAndIncrement();
        return index < n ? index : -1;
    }
}
