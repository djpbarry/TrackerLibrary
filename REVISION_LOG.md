# Revision Log

A chronological record of significant development decisions and changes to
`TrackerLibrary`, and — more importantly — the mistakes made and lessons learned,
so that the same modernisation applied to the remaining sibling projects
(`AdaptDataProcessing`, ADAPT) does not repeat them.

Git commit messages provide the fine-grained record; this file is the distilled,
dated narrative plus the "what not to do again" notes.

## Conventions

- Entries are dated and reference commits and/or `DEVELOPMENT_PLAN.md` phases.
- "Lesson" entries state what went wrong and the rule to apply next time.
- Lessons marked *(inherited)* are carried forward from
  [`IAClassLibrary/REVISION_LOG.md`](https://github.com/djpbarry/IAClassLibrary/blob/development/REVISION_LOG.md),
  whose modernization explicitly lists `TrackerLibrary` as a downstream consumer
  and intends these rules to guide it.

---

## 2026-10-04 — G4: try-with-resources in `PFTracking3D` file I/O

Converted the four file-I/O methods to try-with-resources:

- `writeInitFile` / `writeResultFile` (`BufferedWriter`).
- `readInitFile` / `readResultFile` (`BufferedReader`).

The `read*` methods previously used two nested `try` blocks (open + read) with a
manual `finally`-close; they now open in the resource clause and add
`catch (FileNotFoundException)` (silent `return false`) before
`catch (IOException)` (logged `return false`), preserving the original error
semantics while making the resource lifecycle automatic.

26/26 tests green. Version → `4.0.4` (patch, `refactor`). Next: G8 (redundant
math — `Math.hypot`/`Math.sqrt` swaps behind `TailTracerTest`).

---

## 2026-10-04 — G3: deprecated boxing constructors removed (M7 start)

First Phase G step landed: the 20 `new Double(...)` calls in `TailTracer`
(the `Double(double)` constructor, deprecated-for-removal in Java 21) were
converted to autoboxing. Mechanical, no behaviour change.

- The 20 `Double(double)` compiler warnings are gone.
- 26/26 tests green.
- Version → `4.0.3` (patch, `refactor`).

The only remaining deprecation note is the pre-existing
`ProbabilisticTracker` "uses or overrides a deprecated API" (tracked as a G8
investigation item). Next: G4 (try-with-resources in `PFTracking3D` file I/O).

---

## 2026-10-04 — Phase G scoped (Java 21 modernisation plan)

Thorough code review to scope the next milestone (M7). Added `DEVELOPMENT_PLAN.md`
"Phase G" with a step-by-step, risk-ordered plan mirroring IAClassLibrary's
Phase G. No code changed yet.

**Key finding — a real data race in `PFTracking3D`:** the particle-likelihood
pass (`updateParticleWeights`) hand-rolls a `Thread[]` pool; the shared
work-stealing counter `mControllingParticleIndex` is an *outer* instance field
(`:1316`), but `getNewParticleIndex()` is `synchronized` on the *inner*
`ParallelizedLikelihoodCalculator` instance (one per thread), so the mutex does
not protect the shared counter. Two threads can return the same particle index.
Fix planned under G1 (convert to `Runnable` + `AtomicInteger` + `ExecutorService`).

**Other Phase G items (evidence-backed):**

- G3: 20 `new Double(...)` deprecation-for-removal sites in `TailTracer`.
- G4: manual reader/writer close in `PFTracking3D` file I/O (4 methods) → try-with-resources.
- G6: 3 `instanceof` sites in `TrackMateTracker` → pattern matching.
- G8: `Math.pow(x,0.5)`→`Math.sqrt`, `Math.pow(x²+y²,0.5)`→`Math.hypot` (small).
- G2/G5 already done (D6/D4); `msdPlot`/`globalMSD` left static by design.

---

## 2026-10-04 — Cross-repo review of IAClassLibrary (plan + docs)

Reviewed the sibling's `DEVELOPMENT_PLAN.md`/`AGENTS.md`/`REVISION_LOG.md` to
fold its newer decisions into this project's roadmap. Outcomes:

- **Namespace rename (IAClassLibrary Decision 7) added to our plan** as the main
  open M6 item — `net.calm.*` → `io.github.djpbarry.*` (package root + Maven
  `groupId`), to be executed in lockstep across IAClassLibrary, TrackerLibrary,
  ADAPT, and `AdaptDataProcessing`. This is a breaking change not previously
  tracked here.
- **Conventional Commits + version-bump-on-every-change** recorded as an open
  decision (our Decision 7) — the sibling does it (no `-SNAPSHOT`), we do plain
  semver tags.
- **Lessons L12–L14 carried forward** (marked inherited) from the sibling's
  L13–L15: field-lifetime resources → `AutoCloseable`; version-bump discipline;
  verify upstream is `public` before declaring a copy redundant.
- Confirmed our `v2.0.1` pin is current (their working tree is `2.0.21`, last
  release still `v2.0.1`) and that we are already aligned on `pom-scijava:45.1.0`,
  Java 21, TrackMate 8.0.0, `central`, `jitpack.yml` JDK 21, and `.gitattributes`.

No code change; docs only.

---

## 2026-09-28 — IAClassLibrary re-pinned to v2.0.1 (B4) + patch release 4.0.2

The IAClassLibrary dependency was pinned to a raw JitPack commit hash
(`37a1be016a`) rather than a tagged release — an error, since a commit pin is
not a reproducible semver coordinate and can drift silently. Re-pinned to the
latest JitPack tag, `v2.0.1`.

| Change | What |
|---|---|
| IAClassLibrary pin | `37a1be016a` → `v2.0.1` (latest tag; `v2.0.0` builds `Error` on JitPack, so `v2.0.1` is the first clean tag). |
| `DataStatistics` migration | The 3 `IAClasses.DataStatistics` sites in `ParticleTrajectory` moved to `org.apache.commons.math3.stat.descriptive.DescriptiveStatistics`. The std-dev sites use `Math.sqrt(getPopulationVariance())` to preserve the deprecated class's population (÷N) semantics, which differ from `DescriptiveStatistics.getStandardDeviation()` (sample, ÷(N−1)). |
| `ProgressDialog` migration | `TrajectoryBuilder` and `TrajectoryBridger` now use ImageJ's native `IJ.showProgress(m, n)` loop + `IJ.showProgress(1.0)` clear instead of `IAClasses.ProgressDialog`. |
| Version | `4.0.1` → `4.0.2` (patch — dependency re-pin only; no change to TrackerLibrary's own public API). |

`mvn verify` passes on JDK 21 (26/26 tests).

### Lesson

No new lesson — this is the B4 re-pin the plan already flagged, executed once
`v2.0.1` shipped. Reinforces the existing rule: pin dependencies to **tagged
releases**, not commit hashes (L2 in spirit), and land the deprecated-API
migrations in the same pass as the re-pin so the build gate proves them.

---

## 2026-09-27 — Release 4.0.0 (B2, M6 kick-off)

Version promoted `4.0.0-SNAPSHOT` → `4.0.0` and the malformed SCM URL
(`github.com/github.com/…`) corrected, so the project can be tagged and served by
JitPack to ADAPT/`AdaptDataProcessing`.

- **Tag:** `v4.0.0` (semver, no elided patch zero — L2 satisfied).
- **License:** already resolved to GPL-3.0-or-later before tagging — L3
  satisfied.
- **Release mechanism:** the first tag was cut manually; `maven-release-plugin`
  is configured with `tagNameFormat=v@{project.version}` for future releases.
- **Pending verification (L5):** confirm JitPack builds the tag under
  `pom-scijava:45.1.0` (the sibling IAClassLibrary builds `ok` with the same
  parent and no `jitpack.yml`, so this is expected to work). If the JitPack build
  fails on the JDK/enforcer, the likely fix is a `jitpack.yml` with
  `jdk: [openjdk21]`.

  **Follow-up (same day):** the first `v4.0.0` JitPack build **failed** — JitPack
  defaults to JDK 8 (`Java version: 1.8.0_292`), and the enforcer's
  `RequireJavaVersion` rule rejects it (needs `[21,)`). Fixed by adding a
  `jitpack.yml` with `jdk: [openjdk21]`. See **L11**.

  **Second follow-up:** the `jitpack.yml` commit (`fdae91d`) builds **ok** on
  JitPack, confirming the fix. However, force-moving the `v4.0.0` tag did not
  invalidate JitPack's cached `ref → commit` mapping — the `v4.0.0` name kept
  serving the old failed build. Resolution: since `v4.0.0` never shipped, bump
  to **`4.0.1`** and tag `v4.0.1` (a fresh tag name JitPack has never cached),
  which builds clean. *(This is a JitPack quirk worth remembering — see L11.)*

---

## 2026-09-26 — M1 (Foundations) landed

**Phases A1–A3 and B1–B3 complete.** The build now targets Java 21, the license
is resolved to GPL-3.0-or-later, and the project is shadowed by a reproducible
Maven wrapper matching IAClassLibrary. `mvn verify` passes on JDK 21.

| Change | What |
|---|---|
| A1 | Added Maven wrapper (`mvnw`, `mvnw.cmd`, `.mvn/wrapper/maven-wrapper.properties`) pinned to Maven **3.9.9** (mirrors IAClassLibrary). |
| A2 | CI hardened: JDK 21 (Temurin), `checkout@v4` / `setup-java@v4` with Maven caching, now runs `./mvnw verify` instead of a system `mvn`. Dropped the vestigial `PAT` / `mvn_settings.xml` GitHub Packages flow. |
| A3 | `.gitignore` expanded (`.idea/`, `*.iml`, `.junie/`, `.crush/`); stale `/nbproject/private/` removed. |
| B1 | `pom.xml` license → `gpl_v3` / "GNU General Public License v3.0 or later". |
| B2 | Version `3.0.10` → `4.0.0-SNAPSHOT`; added `maven-release-plugin` with `tagNameFormat=v@{project.version}`. |
| B3 | Parent `pom-scijava:35.0.0` → `45.1.0`; `scijava.jvm.version=21`; removed the explicit `TrackMate:7.10.0` pin (parent resolves 8.0.0). |

### New lessons (see `DEVELOPMENT_PLAN.md` "M1 deviations")

1. **`pom-scijava:45.1.0` drops implicit Maven Central.** JitPack then returned
   empty `xml-apis-ext` artifacts and broke the `BanDuplicateClasses` enforcer
   rule. Fix: declare `central` explicitly in `<repositories>`. Flagged for the
   sibling projects' parent bumps.
2. **TrackMate is parent-managed**, so "bump to 8.0.0" in practice meant
   *remove* the pin, not set `8.0.0`.
3. **`.gitattributes` was required** to keep `mvnw` on LF endings for Linux CI
   (not anticipated by the plan).
4. **Wrapper Maven pin is a convenience, not a floor** — the enforcer's
   `RequireMavenVersion` is the real gate. 3.9.9 was chosen purely for
   consistency with IAClassLibrary.
5. **The `mvnw` wrapper must be committed with its executable bit set.** It was
   initially committed as mode `100644` (no `+x`), so the Linux CI runner's
   `./mvnw` invocation failed instantly with "Permission denied" and the whole
   build job died in ~4 seconds. Fixed with `git update-index --chmod=+x mvnw`.
   This is a recurring failure mode — see **L10**.

### Deferred (not blocking)

- **B4 (IAClassLibrary `v2.0.0` re-pin)** — no `v2.0.0` tag on JitPack yet.
  Maintainer deprioritised it: the Javadoc at `djpbarry.github.io/IAClassLibrary/`
  is the compatibility target, not the pinned artifact. Still on `37a1be016a`;
  re-raise only if the deprecated `DataStatistics` / `ProgressDialog` migration
  actually needs the new API.

---

## 2026-09-26 — M5 docs (README + Javadoc) landing

- **`README.md`** expanded from a single heading to a full overview: the two
  tracking paradigms, build instructions (wrapper on JDK 21), dependencies,
  package map, consumption via JitPack, and license. Added Build / Javadoc /
  JitPack / commit-activity / license **badges** (mirroring IAClassLibrary).
- **Javadoc:** added class/method docs to the public API surface most likely
  touched downstream — `TrajectoryBuilder.updateTrajectories`,
  `TrackMateTracker.track`/`updateTrajectories`, `TrajectoryBridger.bridgeTrajectories`
  (plus `UserVariables`, already documented in the D6 pass).
- **Javadoc publishing:** added `.github/workflows/javadoc.yml` (mirrors the
  sibling) to publish to `djpbarry.github.io/TrackerLibrary/`, and configured
  `maven-javadoc-plugin` in `pom.xml` with `doclint=none` + `quiet=true`.
  Without this the Javadoc build **failed** — two `malformed HTML` errors from
  unescaped `<david.barry at crick.ac.uk>` in `@author` tags, plus ~100
  doclint warnings. The sibling works around the same class of issues with the
  exact same plugin configuration.

---

## 2026-09-26 — D6 resolved (Option A): static mutable state → instance state

Decision 4 of the plan was settled in favour of **Option A** (instance-based
state). This unblocks the remaining D2 decomposition work and makes the
tracking code unit-testable without global-state bleed.

| Change | What |
|---|---|
| `ParticleTrajectory.scale` | `public static double` → `protected double` instance field, defaulting to `1.0`. The no-arg constructor (used by `TrackMateTracker`) no longer inherits whatever `scale` a previous trajectory last set — a latent hidden-global bug fixed. |
| `UserVariables` | Converted from a `public static` singleton to an instance holder: instance (non-static) fields, public constructor, `getInstance()` (lazy) and `setInstance(...)` (inject/reset). In-repo consumers (`TrajectoryBuilder`, `TrajectoryBridger`, `ParticleTrajectory`) now route through `UserVariables.getInstance()`. The `public static final int` constants (`RANDOM`, `MAXIMA`, etc.) are unchanged. |
| `ParticleTrajectory.calcMSD` | Extracted the pure MSD computation into ImageJ-free `calcMSDValues(...)`; the `Plot` rendering/accumulation stays in `calcMSD` as a thin layer (D2). |
| Test | Added `UserVariablesTest` (singleton + instance-isolation). Suite now **22/22**. |

### TailTracer geometry (D2, same pass)

Made the four pure geometry helpers in `TailTracer` — `normalizedVector`,
`normalVector`, `intersections`, `intersections2` — `static` (they never touched
instance state) and added `TailTracerTest` (4 tests). The ImageJ-bound helpers
(`candidatesIntensityCentre`, `findStartingVector`, `trace`, `getX/Y/TailIntens`)
were left instance-bound. This completes the `TailTracer` leg of D2; only the
`PFTracking3D` core / file I/O / GUI split remains outstanding.

### PFTracking3D static helpers extracted (D2, same pass)

Extracted the eight `public static` helpers from `PFTracking3D`
(`copyStateVector`, `copyParticleVector`, `searchLocalMaximumIntensityWithSteepestAscent`,
`addImage`, `initArrayToValue`, `getAFrameCopy`, `getSubStackFloatCopy`,
`getSubStackFloat`) into a new `ParticleFilterUtil` class, leaving thin
`public static` delegating methods in `PFTracking3D` so the inherited-API surface
is unchanged for the subclasses and any downstream callers (Lesson L4).
`PFTracking3D` shrank from ~2206 to ~2095 lines. Added `ParticleFilterUtilTest`
(4 tests: deep copies, elementwise add, array fill). Suite now **26/26**.

The remaining `PFTracking3D` clusters (file I/O ↔ `mStateVectors`/`mFrameOfInitialization`/`mDimensionsDescription`,
and the `ImageCanvas` inner classes) are field-coupled and part of the protected
API, so extracting them would churn call sites for little headless-testability
gain — parked as a lower-value follow-up rather than forced in this pass.

### Deliberately not changed

- `msdPlot`, `plotLegend`, `globalMSD` remain static: they are the **population
  MSD accumulator** (one shared chart aggregating across a run, exposed via
  `drawGlobalMSDPlot()`/`getMsdPlot()`/`resetMSDPlot()`). That is global by
  design, and `ij.gui.Plot`-bound, so headless-testing it adds no value. The
  *testable* math was extracted (above) instead.

### Compatibility note (downstream impact)

This is a **breaking change for downstream callers**: the old `public static`
`UserVariables.getX()`/`setX(...)` methods are gone, replaced by
`UserVariables.getInstance().getX()` (and `setInstance(...)` for injection).
ADAPT / `AdaptDataProcessing` must be updated in the same coordinated pass.
If a smoother migration window is desired, re-add `public static` shim methods
delegating to `getInstance()` — deferred pending the maintainer's preference.

### Lesson

No new lesson — this followed the plan's D6 route as designed. The only
reminders: (1) a source-level breaking change to a library's public API needs an
explicit downstream-usage note (Lesson L4 in spirit); (2) toggle the
"waiting-for-decision" notes in `DEVELOPMENT_PLAN.md` as soon as a decision
lands, not on a later pass.

---

## 2026-09-26 — M2, M3, and M4 (partial) landed

### M2 — Dead-code removal (Phase D1)

Stripped commented-out/experimental code across four files:

| Commit | What |
|---|---|
| — | `TrajectoryBuilder`: removed the legacy combinatorial scoring path (`getMinScores`, `getMinScoreIndices`, `calcScore`, `getNCombs`, `getFirstResult`, `increment`, `allUnique`, commented-out `calcAllPossibleCombs`) and ~20 lines of debug `println`s; dropped the now-unused `DescriptiveStatistics` import. |
| — | `ProbabilisticTracker`: removed the commented-out `main` (hardcoded `C:\Users\barry05\...`) and `mouseReleased` override. |
| — | `PFTracking3D`: removed six commented-out debug/profiling `println` blocks. |
| — | `UserVariables`: removed commented-out `c1Index`/`c2Index`/`channels`/`c2CurveFitTol`/`prevRes`/`medianThresh` fields and accessors. |

Two deliberate deferrals, both recorded:

- The **live** `System.out`/`System.err` diagnostics in `PFTracking3D` → moved to D4.
- `UserVariables.RED`/`GREEN`/`BLUE`/`FOREGROUND` are now unreferenced in-repo but
  left as `public static final` — public API (Lesson L4).

### M3 — Test harness (Phase C)

- Added `org.junit.jupiter:junit-jupiter-api` + `-engine` (test scope, version
  managed by parent → 5.13.4). Note the parent's `junit.version` is 4.13.2
  (JUnit 4), so JUnit 5 required explicit declaration — mirror IAClassLibrary.
- 4 test classes / 16 tests, all headless-safe (no ImageJ runtime):
  `NonIsoGaussianTest`, `FluorophoreTest`, `DecayingFluorophoreTest`,
  `ParticleTrajectoryTest`.
- One test initially failed and exposed a real behavioural detail: `addPoint`
  links newest→oldest, so `getDisplacement` must walk from `getEnd()`, not
  `getStart()`.

### M4 — Refactor (D3–D5 done; D2 gated on D6)

- **D4 (error handling):** 8× `printStackTrace()` → `IJ.handleException(...)`;
  live `System.out`/`System.err` → `IJ.log(...)`; the no-op `e.toString()` and
  the `catch (NullPointerException) { /*do nothing*/ }` now carry explanatory
  comments.
- **D3 (dedupe):** `addBackgroundToImage`, `addFeaturePointTo3DImage`,
  `calculateExpectedZPositionAt` promoted to `protected` in `PFTracking3D` and
  removed from the three subclasses. Kept the defensively-bounded loop variant
  (superset, so behaviour is unchanged) and the `aGhostImage` parameter
  (call-site compatibility).
- **D5 (surface):** 9 `@Override`s added to `ProbabilisticTracker`; no raw
  types; `LinearMovementFPTracker3D` tabs → 4 spaces.
- **D2 blocked** on the D6 decision (see below).

### New lesson from M3/M4

1. **Imports can be "used" transitively without being obvious.** Removing
   `calculateExpectedZPositionAt` from `FPTracker3D`/`LinearMovementFPTracker3D`
   looked like it orphaned `import ij.ImageStack;`, but both files still use
   `ImageStack` elsewhere (`mouseReleased` → `getAFrameCopy`, `autoInitFilter`).
   Always recompile before deleting an import; `findstr` on multiple files via
   `cmd.exe` silently returns "no match" on a parse error. *(Two round-trips of
   compile failures were avoided only by the clean-test gate.)*

---

## 2026-09-25 — Review & plan kick-off (pre-M1)

Initial full review of `TrackerLibrary` and creation of the modernization plan.
This repo is the first of the three IAClassLibrary downstream consumers to be
brought into the coordinated modernization pass, so its plan adopts the
already-resolved sibling decisions (GPL-3.0, Java 21, TrackMate 8.0.0, semver)
rather than re-deriving them.

| Commit | What |
|---|---|
| — | Added `AGENTS.md` (architecture, conventions, gotchas). |
| — | Added `DEVELOPMENT_PLAN.md` (six-phase roadmap). |
| — | Added this `REVISION_LOG.md`. |

Review findings recorded for the plan (see `DEVELOPMENT_PLAN.md` "Verified
gotchas"):

- **Three-way license inconsistency:** `pom.xml` says BSD-2, `LICENSE` is
  GPL-3.0, and source headers are a mix of GPL-2, NetBeans stubs, and no header.
  Resolved by inheriting IAClassLibrary Decision 1 (GPL-3.0-or-later).
- **No tests, no Maven wrapper, no tags, no README content** — the repo is a
  bare library with only a build and CI.
- **Hidden global mutable state:** `ParticleTrajectory.scale` (public static,
  set in the constructor), `msdPlot`/`globalMSD` (static), and the `UserVariables`
  static settings holder.
- **Copy-pasted particle-filter methods** across `ProbabilisticTracker`,
  `FPTracker3D`, and `LinearMovementFPTracker3D`.
- **Inconsistent `@Override`** (missing entirely in `ProbabilisticTracker`) and
  **tabs vs spaces** (`LinearMovementFPTracker3D` uses tabs).
- **Ineffective exception handling:** a no-op `e.toString()` in
  `ParticleTrajectory.projectVelocity`, a `catch (NullPointerException) { /*do
  nothing*/ }` in `PFTracking3D.paint()`, and ~8 `printStackTrace()` calls.
- **Active debug output in `main`:** live `System.out`/`System.err` in
  `PFTracking3D` and a commented-out `main` with a hardcoded Windows path in
  `ProbabilisticTracker`.

### 2026-09-25 — IAClassLibrary API cross-check (authoritative Javadoc)

Cross-checked this repo's `net.calm.iaclasslibrary.*` usage against the
authoritative Javadoc at https://djpbarry.github.io/IAClassLibrary/
(`2.0.0-SNAPSHOT`). The `Particle`/`IsoGaussian`/`ParticleArray`/
`NonIsoGaussianFitter` surfaces match; two findings flagged:

- **`IAClasses.DataStatistics` and `IAClasses.ProgressDialog` are `@Deprecated`**
  in the current IAClassLibrary API but still used here (`ParticleTrajectory`,
  `TrajectoryBuilder`, `TrajectoryBridger`). Migration is a prerequisite for the
  `v2.0.0` re-pin (recorded in `DEVELOPMENT_PLAN.md` B4).
- **`Particle.getFrameNumber()` Javadoc is wrong** — it reads "This particle's
  z-position within an image stack" but returns the frame/time index. Flagged for
  the maintainer to review independently in the IAClassLibrary project.

---

## 2019–2020 — Mavenisation & early history

From git history (pre-modernization):

| Commit | What |
|---|---|
| `0300392` | Mavenised the project and removed the reference to the old `Revision` class. |
| `e99bee1` | Updated the POM to deploy to GitHub. |
| `51a2a25` | Added `mvn_settings.xml` (GitHub Packages repo, PAT-authenticated). |
| `e2cdcd2` | Added the GitHub Actions CI workflow (`maven.yml`). |
| `23f51e9` | Deleted old NetBeans files (`nbproject/`). |
| `0b19cdf` | Added the (one-line) `README.md`. |
| `b21f9db` | Added the GPL-3.0 `LICENSE`. |

Note: `nbproject/` was deleted in `23f51e9`, but `.gitignore` still lists
`/nbproject/private/` (stale entry, see `DEVELOPMENT_PLAN.md` A3).

---

## Lessons learned (mistakes to avoid in sibling projects)

### L1 — Decide the Java target *before* touching build/CI *(inherited)*

IAClassLibrary pinned Java 11 and reversed to Java 21 two commits later because
the coordinated move was not settled upfront.

**Rule:** resolve the cross-project Java / parent-POM / TrackMate version and the
dependency pins in writing before any build or CI edit. For `TrackerLibrary`
this means confirming Java 21 + TrackMate 8.0.0 + `pom-scijava:45.1.0` and the
`v2.0.0` IAClassLibrary re-pin in one pass (Phase A/B together).

### L2 — Use consistent semver tag names *(inherited)*

IAClassLibrary's only tag was `v1.032` (a mislabel of `v1.0.32`) while its POM
said `1.0.37`. `TrackerLibrary` has **no tags at all**, so it starts clean.

**Rule:** use `vX.Y.Z` tags from the first tag onward; never elide the patch
zero. Decide the major bump (`3.0.10` → `4.0.0`) before tagging.

### L3 — Resolve the licence before tagging *(inherited)*

IAClassLibrary shipped a BSD-2 `pom.xml` with GPL-3 sources and no `LICENSE`.
`TrackerLibrary` is worse: BSD-2 POM, GPL-3.0 `LICENSE`, and GPL-2 / NetBeans /
missing source headers.

**Rule:** add/align the root `LICENSE`, `pom.xml` (`<licenses>`,
`license.licenseName`, `license.copyrightOwners`), and treat the source-header
tidy-up as a separate, deferred item so it does not block tagging.

### L4 — "No references in this repo" ≠ "dead code" *(inherited)*

IAClassLibrary nearly deleted `LocationAgnosticBioFormatsImg` because it had no
in-repo references; it is public API. `TrackerLibrary`'s simulation classes and
`TailTracer` are public API with the same risk.

**Rule:** in a library, removing a public symbol requires an external-usage check
and a deprecation window, not just a grep of this repo.

### L5 — Check JitPack compatibility when moving the parent POM *(inherited)*

The parent `pom-scijava` move broke JitPack builds for IAClassLibrary.

**Rule:** verify the parent-POM version is consumable via JitPack before adopting
`45.1.0`; document any downgrade and its reason.

### L6 — Experimental code must not leak into `main` *(inherited)*

`TrackerLibrary` has live `System.out`/`System.err` debug in `PFTracking3D` and a
commented-out `main` with a hardcoded `C:\Users\barry05\...` path, plus a large
block of commented-out combinatorial scoring in `TrajectoryBuilder`.

**Rule:** gate experiments behind a flag or keep them on a branch; strip debug
`println`/hardcoded-path code before merge (recoverable from git).

### L7 — Clean up *all* legacy IDE artifacts, not just the obvious ones *(inherited)*

NetBeans files were deleted, but the `.gitignore` still references them and the
untracked `.idea/`, `.junie/`, and `.crush/` directories are not ignored.

**Rule:** after removing an IDE/legacy build, sweep for and remove or gitignore
the remaining config and output files.

### L8 — Pin the canonical branch explicitly *(inherited)*

Work is on `master` and `origin/HEAD` points at `master` (consistent today), but
the plan must name this explicitly.

**Rule:** name the canonical branch in the plan and keep `origin/HEAD` consistent
before tagging.

### L9 — Java LSP / static-analysis setup *(inherited, parked)*

IAClassLibrary could not get a Java LSP (JDTLS) running under Crush across drives
(`H:` project, `C:` Python/JDK). For the dead-code sweep and `@Override` pass
planned here, fall back to IntelliJ's "unused declaration" / "missing @Override"
inspections rather than relying on a Crush LSP.

### L10 — Commit the Maven wrapper with its executable bit set

After adopting the wrapper in M1, `mvnw` was committed with mode `100644` (no
`+x`). The Linux CI runner invoked `./mvnw …` and failed instantly with
"Permission denied", killing the whole build in ~4 seconds before Maven even
started. This is easy to miss on Windows where the executable bit is not part of
the filesystem (and where `mvnw.cmd` is used instead), so the problem only
surfaced on the Ubuntu runner.

**Rule:** whenever a repo gains a Maven/Gradle wrapper, verify it is tracked with
the executable bit (`git ls-files -s mvnw` must show `100755`, not `100644`). On
Windows, set it explicitly with `git update-index --chmod=+x mvnw` as
`chmod +x` is a no-op there. This applies to `AdaptDataProcessing` and ADAPT when
they adopt the wrapper: confirm `100755` *and* LF line endings (see M1 deviation
3) in the same pass. *(Recurring — has bitten this family of projects before.)*

### L11 — Pin the JDK in `jitpack.yml` for Java 21 projects

JitPack's default build image runs **JDK 8**. A Java 21 project (like
`TrackerLibrary` after the `pom-scijava:45.1.0` move) fails JitPack's build at
the enforcer's `RequireJavaVersion` rule — not because the code is wrong, but
because the runner is on `1.8.0_292`. The failure only surfaces on the first
real JitPack request, since the local `mvnw verify` runs on the developer's JDK 21.

**Rule:** any Java 21 (or otherwise non-default-JDK) project consumed via JitPack
must carry a `jitpack.yml` with `jdk: [openjdk21]` (or the matching version).
Don't assume the sibling already has one — IAClassLibrary's `development` branch
does *not*, so its older (pre-Java-21) commits build `ok` while a fresh Java 21
tag would not. Add `jitpack.yml` in the same pass as the Java 21 / parent-POM
move, before cutting the release tag.

Corollary: **JitPack caches build results by ref name.** Force-moving a tag does
not reliably invalidate the cached `ref → commit` mapping, so a "fixed" tag can
keep serving the stale failed build. If a release tag fails on JitPack, prefer
cutting a fresh tag (next patch version) over force-moving the broken one — a
never-seen tag name always triggers a clean build.

### L12 — Field-lifetime resources need `close()`/`AutoCloseable`, not try-with-resources *(inherited)*

A resource held as a field and used across many methods (e.g. a Bio-Formats
`ImageReader` opened in one method, consumed in others) cannot be wrapped in a
method-local try-with-resources. The correct fix is to make the owning class
`implements AutoCloseable` and expose a `close()`.

**Rule:** distinguish method-scoped resources (try-with-resources) from
object-lifetime resources (implement `AutoCloseable` + `close()`); never just
ignore a field-held reader/stream because it "can't be wrapped in a try".

### L13 — Version via Conventional Commits (bump on every change) *(inherited)*

IAClassLibrary adopted Conventional Commits with the `pom.xml` `<version>` bumped
on **every** code change (`fix`/`refactor`/`chore`/`docs`/`test` → patch,
`feat` → minor, breaking → major) and no `-SNAPSHOT` suffix. Leaving the pom on a
stale version across a batch of commits breaks the version↔tag mapping.

**Rule:** decide the versioning discipline up front; if adopting the sibling's
scheme, bump the version on every change and never leave a `-SNAPSHOT` or stale
version. *(Open question for TrackerLibrary — see `DEVELOPMENT_PLAN.md` Decision 7.)*

### L14 — Verify the upstream method is public before declaring a copy redundant *(inherited)*

A "copied from X" method is only a redundant reimplementation if X's equivalent is
actually `public`/callable and still behaviourally equivalent. IAClassLibrary
found two cases where a grep-based "redundant" verdict was wrong: `OverlayToRoi`
(no public `OverlayCommands.overlayToRoi` exists) and `ImageBlurrer` (already a
thin delegation).

**Rule:** before replacing a "copied from X" method, confirm X's replacement is
public/callable and equivalent — inspect the dependency's class/source (e.g.
`javap`), don't trust a grep-only redundancy verdict.
