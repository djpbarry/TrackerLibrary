[![Build](https://github.com/djpbarry/TrackerLibrary/actions/workflows/maven.yml/badge.svg)](https://github.com/djpbarry/TrackerLibrary/actions/workflows/maven.yml) [![Javadoc](https://img.shields.io/badge/docs-Javadoc-blue)](https://djpbarry.github.io/TrackerLibrary/) [![JitPack](https://jitpack.io/v/djpbarry/TrackerLibrary.svg)](https://jitpack.io/#djpbarry/TrackerLibrary) ![Commit activity](https://img.shields.io/github/commit-activity/y/djpbarry/TrackerLibrary?style=plastic) ![License](https://img.shields.io/github/license/djpbarry/TrackerLibrary?color=green&style=plastic)

# TrackerLibrary

A Java library for particle tracking in fluorescence microscopy. It is designed
to run inside [ImageJ/Fiji](https://imagej.net/software/fiji/) and builds on the
[TrackMate](https://imagej.net/plugins/trackmate/) framework.

TrackerLibrary provides two independent tracking approaches:

- **`ParticleTracking`** — deterministic nearest-neighbour linking plus a
  TrackMate wrapper. Detections are greedily linked frame-by-frame into
  `ParticleTrajectory`s, which also compute analytics (MSD / diffusion
  coefficient, directionality, curvature, angle/step spread, fluorophore
  ratio).
- **`ProbabilisticTracking`** — a sequential Monte Carlo (particle filter)
  tracker ported from ETH Zurich (Janick Cardinale, 2008).

It also ships fluorophore/motion simulation classes for generating synthetic
sequences.

## Building

This is a **library**, not an executable (`main-class` is `None`). Downstream
Fiji plugins consume its classes.

JDK 21 + Maven (a wrapper is committed):

```bash
./mvnw verify
```

On Windows use `mvnw.cmd`. Requires an ImageJ/TrackMate runtime for anything past
compilation/tests (there is no `main` to run).

## Dependencies

- `sc.fiji:TrackMate` (version parent-managed; currently resolves to 8.0.0)
- `com.github.djpbarry:IAClassLibrary` (from JitPack) — sibling library providing
  `net.calm.iaclasslibrary` (`Particle`, `IsoGaussian`, `Utils`, `Region`, …)
- `org.apache.commons.math3`, `org.jgrapht` (transitive)

## Package layout

```
net.calm.trackerlibrary
├── ParticleTracking/        # deterministic linking + analytics + simulation
│   ├── ParticleTrajectory        # linked-list trajectory + analytics
│   ├── TrajectoryBuilder         # nearest-neighbour linking
│   ├── TrackMateTracker          # wraps TrackMate's SparseLAPTracker
│   ├── TrajectoryBridger         # (in Trajectory/) re-links broken segments
│   ├── UserVariables             # runtime configuration (instance-based)
│   ├── Fluorophore + subclasses  # synthetic fluorophore/motion simulators
│   ├── NonIsoGaussian, RandomDistribution
│   └── TailTracer                # traces a filament/tail along an image
├── ProbabilisticTracking/   # particle-filter (sequential Monte Carlo) tracker
│   ├── PFTracking3D              # abstract PlugInFilter base
│   ├── FPTracker3D               # random-walk feature point tracker
│   ├── LinearMovementFPTracker3D # inertia-constrained tracker
│   ├── ProbabilisticTracker      # 7-dim state, intensity-aware tracker
│   └── ParticleFilterUtil        # shared static helpers
└── Trajectory/
    └── TrajectoryBridger
```

## Consuming

The artifact is published via JitPack:

- Group: `net.calm`
- Artifact: `trackerlibrary`

See the [IAClassLibrary Javadoc](https://djpbarry.github.io/IAClassLibrary/) for
the authoritative `net.calm.iaclasslibrary.*` API surface that this library
extends.

## License

GPL-3.0-or-later (see [`LICENSE`](LICENSE)).
