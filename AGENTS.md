# AGENTS.md

## Overview

TrackerLibrary is a **Java library** (not an executable) for particle tracking in
fluorescence microscopy. It is written to run **inside ImageJ/Fiji** and leans
heavily on ImageJ's `ij.*` API and the TrackMate tracking framework. There is no
`main` class (`main-class` is `None` in the POM); downstream Fiji plugins
consume these classes.

Two independent tracking approaches coexist here:

- `ParticleTracking` — deterministic nearest-neighbour linking plus TrackMate.
- `ProbabilisticTracking` — a particle-filter (sequential Monte Carlo) tracker
  ported from ETH Zurich (Janick Cardinale).

## Build

Maven project, JDK 21 (via `scijava.jvm.version`). Parent POM is
`org.scijava:pom-scijava:45.1.0`. Build/toolchain is pinned by a committed
Maven wrapper (`mvnw` / `mvnw.cmd`) at Maven 3.9.9, matching IAClassLibrary.

```bash
./mvnw verify
```

The canonical CI command (`.github/workflows/maven.yml`) is:

```bash
./mvnw --batch-mode --update-snapshots verify
```

Notes:

- CI builds on JDK 21 with Maven dependency caching; the wrapper is used
  directly. The old private-repo `mvn_settings.xml` / `-Dinternal.repo.password`
  flow is no longer used (both dependencies — TrackMate via SciJava and
  IAClassLibrary via JitPack — are public). `mvn_settings.xml` still exists but
  is not referenced by CI.
- A `.gitattributes` forces `mvnw` to LF and `mvnw.cmd` to CRLF so the wrapper
  runs on both Linux CI and Windows.
- The explicit `central` repository was added to `pom.xml` because
  `pom-scijava:45.1.0` drops the implicit Maven Central; without it JitPack
  returned empty artifacts for some transitive deps.

## Dependencies

Declared in `pom.xml`:

- `sc.fiji:TrackMate` (Fiji's TrackMate, version parent-managed; currently
  resolves to 8.0.0)
- `com.github.djpbarry:IAClassLibrary:37a1be016a` (from **JitPack**, with
  TrackMate excluded) — a sibling library by the same author providing
  `net.calm.iaclasslibrary` (`Particle.Particle`, `Particle.IsoGaussian`,
  `IAClasses.Utils`, `IAClasses.Region`, `IAClasses.ProgressDialog`,
  `Math.Optimisation.NonIsoGaussianFitter`, etc.).
- `org.apache.commons.math3` (transitive) for statistics, linear algebra, and
  distributions.
- `org.jgrapht` (transitive, via TrackMate) for the tracking graph.

Extra repositories: `maven.scijava.org` (scijava.public) and `jitpack.io`.

### IAClassLibrary API reference (authoritative)

The current IAClassLibrary API is browsable at
**https://djpbarry.github.io/IAClassLibrary/** (documents `2.0.0-SNAPSHOT`).
Treat this Javadoc as the **authoritative source** for the
`net.calm.iaclasslibrary.*` surface, not the raw IAClassLibrary source.

Note the pin (`37a1be016a`) predates `2.0.0-SNAPSHOT`, so the Javadoc reflects
the *target* API after the planned re-pin, not what the current build compiles
against. Two classes this repo still uses are `@Deprecated` there:

- `net.calm.iaclasslibrary.IAClasses.DataStatistics` — used in
  `ParticleTrajectory.java` (3 sites). Prefer
  `org.apache.commons.math3.stat.descriptive.DescriptiveStatistics`, which this
  repo already imports elsewhere.
- `net.calm.iaclasslibrary.IAClasses.ProgressDialog` — used in
  `TrajectoryBuilder.java` and `TrajectoryBridger.java`.

Both must be migrated off before/when re-pinning to `v2.0.0`.

## Package layout

All source lives under `src/main/java/net/calm/trackerlibrary/`.

**Note the unconventional capitalized sub-package names** — `ParticleTracking`,
`ProbabilisticTracking`, `Trajectory`. Follow this naming if you add packages.

```
net.calm.trackerlibrary
├── ParticleTracking/
│   ├── ParticleTrajectory.java      # linked-list trajectory + analytics
│   ├── TrajectoryBuilder.java       # nearest-neighbour linking of detections
│   ├── TrackMateTracker.java        # wraps TrackMate's SparseLAPTracker
│   ├── UserVariables.java           # global static configuration
│   ├── Fluorophore.java             # simulation base class
│   ├── BlinkingFluorophore.java     # simulation (on/off blinking)
│   ├── DecayingFluorophore.java     # simulation (exponential decay)
│   ├── MotileGaussian.java          # simulation (persistent/Brownian motion)
│   ├── ConfinedGaussian.java        # simulation (motion bounded by an ROI)
│   ├── NonIsoGaussian.java          # simulation (anisotropic Gaussian)
│   ├── RandomDistribution.java
│   └── TailTracer.java              # traces filament/tail along image
├── ProbabilisticTracking/
│   ├── PFTracking3D.java            # abstract particle-filter base (PlugInFilter)
│   ├── FPTracker3D.java             # random-walk feature point tracker
│   ├── LinearMovementFPTracker3D.java # inertia-constrained tracker
│   └── ProbabilisticTracker.java    # 7-dim state, intensity-aware tracker
└── Trajectory/
    └── TrajectoryBridger.java       # re-links broken trajectory segments
```

## Architecture and data flow

**`ParticleTrajectory` is the central data structure.** It is a singly linked
list of `Particle` objects (from IAClassLibrary), always appending at the `end`
reference:

- `addPoint(Particle)` calls `particle.setLink(end)` then sets `end = particle`.
- `getStart()` walks the `getLink()` chain to the head.
- `addTrajectory(other)` splices another trajectory's points into this one (used
  by `TrajectoryBridger`).
- `Particle` is external (`net.calm.iaclasslibrary.Particle.Particle`); the
  common surface used here: `getX()`, `getY()`, `getFrameNumber()`,
  `getMagnitude()`, `getRegion()`, `getColocalisedParticle()`,
  `setColocalisedParticle()`, `getFeature()/putFeature()`, `makeCopy()`,
  `setLink()/getLink()`.

**Tracking pipeline (deterministic):**

1. `TrajectoryBuilder.updateTrajectories(...)` iterates detection frames and
   links each detection to the trajectory with the lowest weighted score of
   `position + projected-velocity + morphology` (`mw`, `vw`, `pw` weighting).
   Each detection becomes a new single-point trajectory on first sight.
2. `TrackMateTracker.track(...)` runs TrackMate's `SparseLAPTracker` and
   `updateTrajectories(...)` converts the resulting spots into
   `ParticleTrajectory`s.
3. `TrajectoryBridger.bridgeTrajectories(...)` joins segments whose end/start
   fall within `maxStep` frames, using the same scoring.
4. `ParticleTrajectory` then computes analytics: MSD / diffusion coefficient
   (`calcMSD`), directionality (`calcDirectionality`), curvature
   (`smooth`/`calcSpec`), angle spread, step spread, fluorophore ratio.

**Probabilistic pipeline:** `PFTracking3D` (abstract) implements an ImageJ
`PlugInFilter` particle filter. Subclasses override the abstract contract:

- `generateIdealImage_3D(...)` — render the expected image from a state vector
- `drawFromProposalDistribution(...)` — the dynamics/proposal model
- `paintOnCanvas(...)` / `paintParticleOnCanvas(...)` — visualization
- `getMSigmaOfRandomWalk()`, `getMDimensionsDescription()`,
  `getMDoPrecisionOptimization()`

State-vector layout differs per subclass (e.g. `FPTracker3D` uses
`{x,y,z,intensity}`, `LinearMovementFPTracker3D` and `ProbabilisticTracker` use
`{x,y,z,vx,vy,vz,intensity}`). Read the subclass fields before extending.

## Conventions and gotchas

### Global mutable static state

- `UserVariables` is a **global static settings holder** (spatial/time
  resolution, thresholds, motion model, detection mode, etc.). All access is via
  static getters/setters. There is no per-instance configuration; code reads
  `UserVariables.getX()` directly.
- `ParticleTrajectory.scale` is a `public static` field overwritten in the
  constructor from `spatialRes` — a hidden global side effect. The static
  `msdPlot` / `globalMSD` fields accumulate state across instances until
  `resetMSDPlot()`.

Be careful: these make classes hard to unit-test in isolation and can leak state
between operations.

### Magic-int constants, not enums

- Motion models: `UserVariables.RANDOM = 6`, `DIRECTED = 7`.
- Detection modes: `MAXIMA = 3`, `BLOBS = 4`, `GAUSS = 5`.
- Colocalization flags in `ParticleTrajectory`: `NON_COLOCAL = 0`, `UNKNOWN = 1`,
  `COLOCAL = 2`.

### ImageJ UI is woven throughout

`IJ.log`, `IJ.showMessage`, `GenericDialog`, `ProgressDialog`, `Plot`,
`TextWindow`, and `Roi` appear in the "business logic", not just UI layers.
This means most code **cannot run headless**; it expects an ImageJ runtime. There
is no test harness to work around this.

### Dead / commented-out code

`TrajectoryBuilder` contains large blocks of commented-out legacy scoring code
(the old `getMinScoreIndices` / `getMinScores` combinatorial approach). Don't
assume commented code is unused because of a bug; it was superseded by the
greedy `addTempPoint` approach. `UserVariables` also has many commented-out
fields/methods.

### Licensing is inconsistent — verify before relying on it

- `pom.xml` declares **GPL-3.0-or-later** (`license.licenseName=gpl_v3`), fixed
  in M1 (was BSD-2).
- `LICENSE` file is **GPLv3**.
- Many source file headers still say **GPLv2** (source-header cleanup is a
  deferred, non-blocking item — see `DEVELOPMENT_PLAN.md` B1).

If licensing matters, flag this discrepancy rather than asserting a single
license.

### Style

- Indentation is 4 spaces; tabs appear in the older ETH-ported
  `LinearMovementFPTracker3D.java` (mixed — preserve locally).
- Class headers vary (NetBeans auto-template, GPL boilerplate, or none). Do not
  normalize them unless asked.
- `@author` tags reference `barry05` / `David Barry` / `Dave Barry` / Janick
  Cardinale (ETH).
