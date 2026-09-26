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

## 2026-09-26 — D6 resolved (Option A): static mutable state → instance state

Decision 4 of the plan was settled in favour of **Option A** (instance-based
state). This unblocks the remaining D2 decomposition work and makes the
tracking code unit-testable without global-state bleed.

| Change | What |
|---|---|
| `ParticleTrajectory.scale` | `public static double` → `protected double` instance field, defaulting to `1.0`. The no-arg constructor (used by `TrackMateTracker`) no longer inherits whatever `scale` a previous trajectory last set — a latent hidden-global bug fixed. |
| `UserVariables` | Converted from a `public static` singleton to an instance holder: instance (non-static) fields, public constructor, `getInstance()` (lazy) and `setInstance(...)` (inject/reset). In-repo consumers (`TrajectoryBuilder`, `TrajectoryBridger`, `ParticleTrajectory`) now route through `UserVariables.getInstance()`. The `public static final int` constants (`RANDOM`, `MAXIMA`, etc.) are unchanged. |
| `ParticleTrajectory.calcMSD` | Extracted the pure MSD computation into ImageJ-free `calcMSDValues(...)`; the `Plot` rendering/accumulation stays in `calcMSD` as a thin layer (D2). |
| Test | Added `UserVariablesTest` (singleton + instance-isolation). Suite now **18/18**. |

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
