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
