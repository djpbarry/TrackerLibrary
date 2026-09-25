# TrackerLibrary Development Plan

This plan outlines a multi-phase effort to make `TrackerLibrary` more robust,
maintainable, and consumable as a library, and to bring it into lockstep with
the already-modernised sibling project, [`IAClassLibrary`](https://github.com/djpbarry/IAClassLibrary)
(which lists `TrackerLibrary` as one of its downstream consumers). It is grounded
in a full review of this repository as of the plan's writing.

## Overarching aim

**Modernise a codebase whose core is over a decade old, and align it with the
sibling modernisation already underway in IAClassLibrary.** As with
IAClassLibrary, modernisation takes precedence over backward compatibility with
downstream consumers: ADAPT and `AdaptDataProcessing` will themselves be updated
in the same coordinated pass. Where a choice is between a cleaner, more
maintainable design and preserving a legacy API, prefer the cleaner design.

- Deprecate legacy symbols rather than preserving them indefinitely.
- Prefer deletion of dead/experimental code over keeping it "just in case"
  (it is recoverable from git).
- Record every decision, mistake, and lesson in
  [`REVISION_LOG.md`](REVISION_LOG.md) so the remaining sibling
  (`AdaptDataProcessing`, ADAPT) does not repeat them.

## Current state (context for the plan)

- `TrackerLibrary` (`net.calm.trackerlibrary`) is a **Java library** (not a
  runnable plugin; `main-class` is `None`) in the ImageJ/Fiji ecosystem. It
  provides two independent particle-tracking approaches — deterministic
  nearest-neighbour linking and a probabilistic particle filter — plus
  fluorophore/motion simulation classes for generating synthetic sequences.
- **Build:** Maven, parent `org.scijava:pom-scijava:35.0.0`, group `net.calm`,
  artifact `trackerlibrary`, version `3.0.10`. Declared license **BSD-2**
  (`license.licenseName=bsd_2`). No Maven wrapper.
- **CI:** `.github/workflows/maven.yml` runs a single `mvn verify` on **JDK 11**
  (AdoptOpenJDK) with a `PAT` secret and `--settings mvn_settings.xml` (adds a
  GitHub Packages repo for `djpbarry`).
- **Tests:** none — `src/test/java/` and `src/main/resources/` are empty.
  No lint/format tooling.
- **License:** a **GPL-3.0** `LICENSE` file exists, but `pom.xml` declares
  **BSD-2**, and source headers are a three-way mix of **GPL-2** boilerplate,
  NetBeans "change this header" stubs, and files with no header at all.
- **Dependencies:** `sc.fiji:TrackMate:7.10.0` (explicitly pinned) and
  `com.github.djpbarry:IAClassLibrary:37a1be016a` (JitPack commit hash, with
  TrackMate excluded). The IAClassLibrary pin predates that project's `v2.0.0`
  modernization and its Java 21 / TrackMate 8.0.0 move.
- **Structure:** 17 Java files across three subpackages under
  `net.calm.trackerlibrary` — `ParticleTracking` (11 files), `ProbabilisticTracking`
  (4 files), and `Trajectory` (1 file). Package names are capitalised, which is
  unconventional in Java.

### Key architectural facts

- **Two independent tracking paradigms coexist.**
  - `ParticleTracking` — deterministic linking. `TrajectoryBuilder` greedily
    links detections frame-by-frame; `TrackMateTracker` wraps TrackMate's
    `SparseLAPTracker` and converts its spots into trajectories;
    `TrajectoryBridger` re-joins gap-broken segments.
  - `ProbabilisticTracking` — a sequential Monte Carlo (particle filter) tracker
    ported from ETH Zurich (Janick Cardinale, 2008). `PFTracking3D` is the
    abstract `PlugInFilter` base; `FPTracker3D`, `LinearMovementFPTracker3D`,
    and `ProbabilisticTracker` are concrete subclasses.
- **`ParticleTrajectory` is the central data structure:** a singly linked list of
  `Particle` objects (from IAClassLibrary, itself extending TrackMate `Spot`),
  appended at the `end` reference via `particle.setLink(end)`. It doubles as the
  analytics engine (MSD/diffusion coefficient, directionality, curvature, step
  and angle spread, fluorophore ratio).
- **`UserVariables` is a global static configuration holder.** All runtime
  settings (spatial/time resolution, thresholds, motion model, detection mode)
  are static fields accessed through static getters/setters. There is no
  instance-based configuration.
- **The `Particle` type is external** (`net.calm.iaclasslibrary.Particle.Particle`);
  the common surface used here is `getX()`, `getY()`, `getFrameNumber()`,
  `getMagnitude()`, `getRegion()`, `getColocalisedParticle()`,
  `setColocalisedParticle()`, `getFeature()/putFeature()`, `makeCopy()`,
  `setLink()/getLink()`.
- **Simulation classes** (`Fluorophore` hierarchy, `MotileGaussian`,
  `ConfinedGaussian`, `NonIsoGaussian`, `RandomDistribution`) generate synthetic
  particles. They are part of the public API and may be consumed by downstream
  projects, not just internally.

### Verified gotchas (from the code review)

- **Three-way license inconsistency** — `pom.xml` says BSD-2, `LICENSE` is
  GPL-3.0, and source headers are a mix of GPL-2 (`TrajectoryBuilder`,
  `TrackMateTracker`, `ConfinedGaussian`, `TrajectoryBridger`), NetBeans stubs
  (`Fluorophore`, `BlinkingFluorophore`, `MotileGaussian`, `NonIsoGaussian`,
  `RandomDistribution`, `TailTracer`, `UserVariables`, `ProbabilisticTracker`),
  and no header at all (`ParticleTrajectory`, `DecayingFluorophore`,
  `PFTracking3D`, `FPTracker3D`, `LinearMovementFPTracker3D`).
- **Hidden global mutable state.** `ParticleTrajectory.scale` is a `public
  static` field overwritten in the constructor from `spatialRes`; `msdPlot` and
  `globalMSD` are static and accumulate across instances until `resetMSDPlot()`.
- **`UserVariables` is a static singleton** — cross-cutting, untestable global
  state that every class reads directly.
- **Magic-int constants instead of enums** — motion models (`RANDOM=6`,
  `DIRECTED=7`), detection modes (`MAXIMA=3`, `BLOBS=4`, `GAUSS=5`), and
  colocalisation flags (`NON_COLOCAL=0`, `UNKNOWN=1`, `COLOCAL=2`).
- **ImageJ UI woven into business logic.** `IJ.log`, `IJ.showMessage`,
  `GenericDialog`, `ProgressDialog`, `Plot`, `TextWindow`, `Roi`, and
  `IJ.showProgress` appear throughout "core" code. Most classes cannot run
  headless, which is why there are no tests.
- **Large commented-out dead code in `TrajectoryBuilder`.** The old combinatorial
  scoring approach (`getMinScoreIndices`/`getMinScores`/`calcAllPossibleCombs`)
  is left commented out alongside the live greedy `addTempPoint` path.
- **Copy-pasted private methods across the three particle-filter subclasses.**
  `addFeaturePointTo3DImage`, `addBackgroundToImage`, and (in `FPTracker3D` and
  `LinearMovementFPTracker3D`) `calculateExpectedZPositionAt` are duplicated
  verbatim rather than inherited from `PFTracking3D`.
- **Inconsistent `@Override` usage.** `FPTracker3D` and
  `LinearMovementFPTracker3D` annotate their overrides; `ProbabilisticTracker`
  overrides the same abstract methods with no `@Override` at all.
- **Ineffective/empty exception handling.** `ParticleTrajectory.projectVelocity`
  swallows `Exception` with a no-op `e.toString()`; `PFTracking3D` has an
  `catch (NullPointerException) { /*do nothing*/ }` in `paint()`, and ~8
  `printStackTrace()` calls in its file-I/O methods.
- **Active debug output in `main`.** `PFTracking3D` still has live
  `System.out.println("NAN at ...")` and `System.err.println("wrong argument...")`
  statements, plus a commented-out `main` with a hardcoded
  `C:\Users\barry05\...` path in `ProbabilisticTracker`.
- **Stale `.gitignore`.** It still lists `/nbproject/private/` even though the
  NetBeans files were deleted in commit `23f51e9`; `.idea/`, `.junie/`, and
  `.crush/` are untracked IDE/agent dirs not covered by the ignore rules.
- **`PFTracking3D` is ~1900 lines** — a single `PlugInFilter` base class holding
  the state machine, particle-filter core, file I/O, and GUI wiring in one file.

---

## Phase A — Build & CI modernization

### A1. Add a Maven wrapper

Pin a reproducible toolchain by committing `mvnw`/`mvnw.cmd` +
`.mvn/wrapper/`, mirroring IAClassLibrary's A1. Currently the build depends on a
system Maven.

### A2. Harden CI

`.github/workflows/maven.yml` currently runs a single `mvn verify` on JDK 11.

1. Use the wrapper (`./mvnw verify`) instead of a system `mvn`.
2. Cache Maven dependencies (`actions/cache` on `~/.m2/repository`).
3. Pin the Java target to 21 and build on JDK 21 (see Decision 3), matching
   IAClassLibrary.
4. Split `build` and `test` jobs once Phase C lands.
5. Re-evaluate the `PAT` / `mvn_settings.xml` GitHub Packages flow — the only
   declared dependencies are TrackMate (scijava) and IAClassLibrary (JitPack),
   both public, so the private-repo settings may be vestigial.

### A3. Expand `.gitignore`

Current `.gitignore` only excludes `/nbproject/private/`, `/build/`, `/dist/`,
`/target/`. Add `.idea/`, `*.iml`, `.junie/`, `.crush/`, and OS files. Remove the
now-stale `/nbproject/private/` entry (NetBeans files deleted in `23f51e9`).

---

## Phase B — License, metadata & versioning (follows IAClassLibrary's B)

IAClassLibrary has already resolved the license (GPL-3.0-or-later), the Java
target (21), the TrackMate version (8.0.0), and the semver tag policy. This repo
must adopt the same resolutions rather than re-litigate them.

### B1. Resolve the license

`TrackerLibrary` has the same unresolved contradiction IAClassLibrary fixed:
BSD-2 in `pom.xml`, a GPL-3.0 `LICENSE`, and mixed GPL-2 / NetBeans / missing
source headers. Adopt **GPL-3.0-or-later** (Decision 1, inherited):

1. Correct `pom.xml` (`<licenses>`, `license.licenseName=gpl_v3`,
   `license.copyrightOwners` = Francis Crick Institute / David Barry).
2. Confirm the existing root `LICENSE` matches GPL-3.0-or-later.
3. **Deferred (future tidy-up).** Source-header cleanup (removing the NetBeans
   stubs and normalising the GPL-2 vs no-header split) is a separate, non-blocking
   pass.

### B2. Versioning & tagging

The repo has **no tags** and sits at version `3.0.10`. Adopt semver `vX.Y.Z`
tags (Decision 2, inherited). Decide with the maintainer whether to jump to a new
major (`4.0.0-SNAPSHOT`) to reflect the scale of the modernization, mirroring
IAClassLibrary's jump to `2.0.0`. Consider wiring `maven-release-plugin` with
`tagNameFormat=v@{project.version}` as IAClassLibrary did.

### B3. TrackMate version web (blocks the coordinated move)

This repo pins `sc.fiji:TrackMate:7.10.0` explicitly, while IAClassLibrary now
resolves to **8.0.0** via its parent. Bump to 8.0.0 (or drop the explicit version
and let `pom-scijava:45.1.0` manage it) in lockstep with the Java 21 move.

### B4. Re-pin the IAClassLibrary dependency

`IAClassLibrary` is pinned to JitPack commit `37a1be016a`, which predates its
`v2.0.0` release. Re-pin to the `v2.0.0` tag once it ships, so downstream
consumers of `TrackerLibrary` resolve a tagged, modernised IAClassLibrary.

The authoritative API reference is the Javadoc at
**https://djpbarry.github.io/IAClassLibrary/** (documents `2.0.0-SNAPSHOT`).
Cross-checking this repo's IAClassLibrary usage against it reveals two classes
this repo still uses that are now `@Deprecated`:

- `IAClasses.DataStatistics` (used 3× in `ParticleTrajectory`) — migrate to
  `org.apache.commons.math3.stat.descriptive.DescriptiveStatistics`, already
  imported elsewhere in the same file.
- `IAClasses.ProgressDialog` (used in `TrajectoryBuilder` and
  `TrajectoryBridger`) — migrate to ImageJ's native progress mechanism.

These migrations are prerequisites for the re-pin, not merely cleanup.

---

## Phase C — Introduce tests (the biggest maintainability win)

There is no test infrastructure. Mirror IAClassLibrary's Phase C.

1. **Add a JUnit 5 (Jupiter) harness** (parent-managed `junit-jupiter-api` /
   `-engine`).
2. **Start with pure-logic classes that need no ImageJ runtime:**
   - `DecayingFluorophore.updateMag` / `updateXMag`, `BlinkingFluorophore.updateMag`
     (deterministic with a seeded `Random`).
   - `NonIsoGaussian.evaluate(...)` and its `a`/`b`/`c` construction math.
   - `TailTracer`'s pure geometry helpers (`normalizedVector`, `normalVector`,
     `intersections`, `candidatesIntensityCentre`).
   - `ParticleTrajectory`'s ImageJ-free methods: `getDisplacement`,
     `getNumberOfFrames`, `calcDualScore`, `getFluorRatio`, `getType`, `getStart`.
   - `TrajectoryBuilder`/`TrajectoryBridger` scoring math, if extracted (see D2).
3. **Extract pure logic before testing it.** The MSD/directionality/curvature
   analytics in `ParticleTrajectory` are entangled with `ij.gui.Plot`; separate
   the numeric computation from the rendering so it can be unit-tested (see D2).
4. **Note:** full end-to-end tests are blocked by the ImageJ runtime dependency;
   target headless-safe units first, as IAClassLibrary did.

---

## Phase D — Code hygiene & refactoring

### D1. Remove dead/experimental code (do first)

- Strip the large commented-out combinatorial-scoring blocks in
  `TrajectoryBuilder` (`getMinScoreIndices`, `getMinScores`,
  `calcAllPossibleCombs`, and the surrounding `//` debug `println` blocks). They
  are superseded by the greedy `addTempPoint` path and recoverable from git.
- Remove the commented-out `main` in `ProbabilisticTracker` (hardcoded
  `C:\Users\barry05\...` path) and the commented-out debug `System.out.println`
  blocks in `PFTracking3D`.
- Sweep for other commented-out `IJ.saveAs(...)`/`System.out.println` debug lines
  with machine-specific paths.

### D2. Decompose the largest, most intertwined units

- **`PFTracking3D` (~1900 lines)** is the largest file: split the particle-filter
  core, the file I/O (`writeInitFile`/`readInitFile`/`writeResultFile`/
  `readResultFile`), and the GUI (`ImageCanvas` inner class, button handlers) into
  separate responsibilities.
- **`ParticleTrajectory` (~688 lines)** mixes the linked-list structure with
  analytics and plotting. Extract the pure-math analytics (MSD, directionality,
  curvature, step/angle spread) from the `Plot` rendering so Phase C can test them.
- **`TailTracer` (~400 lines)** is geometry-heavy with magic numbers; extract
  the pure vector/intersection math.

### D3. Deduplicate the particle-filter subclasses

`addFeaturePointTo3DImage`, `addBackgroundToImage`, and
`calculateExpectedZPositionAt` are copy-pasted across `ProbabilisticTracker`,
`FPTracker3D`, and `LinearMovementFPTracker3D`. Promote the shared implementations
into `PFTracking3D` (protected), leaving only the genuinely model-specific
overrides in each subclass.

### D4. Normalise error handling

- Replace the no-op `catch (Exception e) { e.toString(); }` in
  `ParticleTrajectory.projectVelocity` with either real handling or a comment
  explaining the intentional swallow.
- Replace `printStackTrace()` (8× in `PFTracking3D`) and the bare
  `catch (NullPointerException) { /*do nothing*/ }` in `paint()` with a
  consistent logging pattern (e.g. `IJ.log`/`IJ.handleException`), matching the
  IAClassLibrary D6 discipline.
- Remove or convert the live `System.out`/`System.err` debug output in
  `PFTracking3D` (the "NAN at ..." and "wrong argument..." lines).

### D5. Modernise the class surface

- Add missing `@Override` annotations in `ProbabilisticTracker` (its sibling
  classes already have them).
- Replace raw types / add diamond generics where missing (`new ArrayList<>()` vs
  `new ArrayList()`).
- Normalise tabs vs spaces: `LinearMovementFPTracker3D` is tab-indented while the
  rest of the repo uses 4 spaces.

### D6. Resolve the static mutable state (largest architectural decision)

`UserVariables` (global settings) and the static fields in `ParticleTrajectory`
(`scale`, `msdPlot`, `globalMSD`) are the main obstacles to testability and safe
concurrent use. Decide the approach with the maintainer:

- Option A — convert `UserVariables` to an instance passed through the tracking
  pipeline (larger change, most testable).
- Option B — keep the static holder but document it and reset it explicitly in
  tests (smaller change, defers the real problem).

This decision should be made before Phase C test work proceeds too far, since it
changes how state is injected.

---

## Phase E — Documentation

- **`README.md`** is a single `# TrackerLibrary` line. Expand it with: what the
  library is, the two tracking paradigms, how to build (`mvn verify`), how it is
  consumed (JitPack), and the package map.
- **`AGENTS.md`** already documents architecture, conventions, and gotchas. Keep
  it in sync as Phases A–D land.
- **Javadoc:** most classes have minimal or no Javadoc (only `ParticleTrajectory`
  and `PFTracking3D` have any). Add it to the public API surface, especially the
  classes a downstream consumer would touch (`ParticleTrajectory`,
  `TrajectoryBuilder`, `TrackMateTracker`, `TrajectoryBridger`, `UserVariables`).

---

## Phase F — Upstream coordination

`TrackerLibrary` is downstream of `IAClassLibrary` and (transitively) upstream of
ADAPT, so several items are prerequisites for, or follow from, the sibling plans:

| Sibling plan item | TrackerLibrary action | This plan |
|---|---|---|
| IAClassLibrary B1 (license) | adopt GPL-3.0-or-later | B1 |
| IAClassLibrary B2 (tag) | re-pin IAClassLibrary to `v2.0.0` | B4 |
| IAClassLibrary B3 (TrackMate 8 / Java 21) | bump TrackMate to 8.0.0, Java 21 | A2, B3 |
| ADAPT G2 (tag all three) | tag `v4.0.0` (or agreed version) | B2 |

Coordinate the Java-target and TrackMate-version decision with the maintainer of
`IAClassLibrary` (already resolved there) and confirm the same applies here.

---

## Proposed decisions (to confirm with the maintainer)

These follow directly from IAClassLibrary's resolved decisions and should be
confirmed, not re-derived:

1. **License — GPL-3.0-or-later.** Correct `pom.xml` from BSD-2; defer the
   source-header tidy-up (inherited from IAClassLibrary Decision 1).
2. **Versioning — semver, jump to `4.0.0`.** Reflect the scale of modernization
   with a major bump (`4.0.0-SNAPSHOT`), tag `v4.0.0` on release. Use `vX.Y.Z`
   tags (inherited from Decision 2). Confirm the exact major with the maintainer.
3. **Java target — 21, TrackMate 8.0.0.** Upgrade the parent to
   `pom-scijava:45.1.0` and drop the explicit `TrackMate:7.10.0` pin (inherited
   from Decision 3).
4. **Static mutable state — decide Options A vs B** (see D6). This is
   TrackerLibrary-specific and needs a maintainer call.
5. **Public API — "no in-repo reference" ≠ "dead code".** The simulation classes
   and `TailTracer` may be consumed externally (inherited from IAClassLibrary
   Decision 0). Removal of public symbols requires a deprecation window and an
   external-usage check.

---

## Suggested sequencing & milestones

1. **M1 — Foundations (low risk, high value):** `.gitignore`, Maven wrapper, CI
   hardening (Java 21), license + pom metadata fix, re-pin IAClassLibrary to
   `v2.0.0`, bump version and tag. (Phases A, B)
2. **M2 — Dead-code removal:** strip commented-out combinatorial scoring and
   debug output. (Phase D1)
3. **M3 — Test harness:** JUnit 5 + pure-logic unit tests. (Phase C)
4. **M4 — Refactor core:** decompose `PFTracking3D`/`ParticleTrajectory`,
   deduplicate the particle-filter subclasses, normalise error handling. (D2–D5)
5. **M5 — Static-state decision + documentation:** resolve D6, expand `README.md`
   and Javadoc. (Phase E)
6. **M6 — Upstream hand-off:** confirm Java 21 / TrackMate 8 / tag with ADAPT and
   `AdaptDataProcessing`. (Phase F)

Each milestone is independently shippable. M1 is the immediate next step and
unblocks the coordinated downstream modernization.
