package net.calm.trackerlibrary.ProbabilisticTracking;

import org.junit.jupiter.api.Test;

import java.util.concurrent.atomic.AtomicIntegerArray;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ParticleIndexAllocatorTest {

    @Test
    void nextReturnsEachIndexExactlyOnceThenMinusOne() {
        ParticleIndexAllocator allocator = new ParticleIndexAllocator();
        int n = 100;
        for (int i = 0; i < n; i++) {
            assertEquals(i, allocator.next(n));
        }
        assertEquals(-1, allocator.next(n));
        assertEquals(-1, allocator.next(n));
    }

    @Test
    void nextIsThreadSafeUnderConcurrentClaiming() throws InterruptedException {
        int n = 1000;
        int nThreads = 8;
        ParticleIndexAllocator allocator = new ParticleIndexAllocator();
        AtomicIntegerArray seen = new AtomicIntegerArray(n);

        Thread[] threads = new Thread[nThreads];
        for (int t = 0; t < nThreads; t++) {
            threads[t] = new Thread(() -> {
                int index;
                while ((index = allocator.next(n)) != -1) {
                    seen.incrementAndGet(index);
                }
            });
        }
        for (Thread thread : threads) {
            thread.start();
        }
        for (Thread thread : threads) {
            thread.join();
        }

        for (int i = 0; i < n; i++) {
            assertEquals(1, seen.get(i), "index " + i + " should be claimed exactly once");
        }
    }
}
