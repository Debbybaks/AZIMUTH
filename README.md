# AZIMUTH

Astronomical Zonal Intelligence, Mapping & Unified Trajectory Hub.

Canonical source: https://github.com/Debbybaks/AZIMUTH, branch `main`.

Workspace and remote inspection on 2026-10-24 found no source files, commits, or branches and no AZIMUTH Vercel project in the connected team. This repository is the source of truth. Production must be deployed from its verified commit.

Scientific simulation software under development. Model assumptions and verification are documented with the implementation.

## Run and verify

Requires Node.js 22+ and npm. `npm ci`, `npm run dev`, `npm test`, `npm run build`. Vercel uses the Vite static build in `dist/` (ignored by Git). Fonts and Earth texture are self-hosted; the application makes no runtime API calls and needs no secrets.

## Implemented scientific workflows

- Interactive WebGL solar system: Sun, eight planets and Moon; camera orbit/pan/zoom/focus, labels, orbital paths, scale modes, and replay spacecraft tracking.
- Mission designer: approximate Kepler planetary ephemerides, prograde zero-revolution universal-variable Lambert solver, numerical propagation, departure velocity perturbation, calculated encounter/miss, telemetry, energy error, and replay.
- Hohmann calculations, orbital elements ↔ Cartesian states, ecliptic/equatorial/horizon coordinate transforms, patched-conic Earth gravity assist.
- Worker-based porkchop and bounded-grid optimization, seeded Gaussian velocity Monte Carlo, fully interacting three-body experiment, integrator comparison (Euler, Heun, RK4, velocity Verlet).
- Spherical-Earth ground tracks, satellite constellations, geometric footprints, simplified local atmospheric decay rate, planetary reference-site great-circle distances, explicitly synthetic terrain experiment.
- Local experiment notebook with complete computed trajectories, notes, duplication, renaming, and recomputation. Comparison with the previous mission. JSON/CSV export and downloadable HTML scientific report with plots.
- Mars intercept challenge: exposed score formula, actual computed tolerance, no guaranteed success.

## Scientific limitations and remaining scope

This release is a research/educational simulator, not operational flight software. See the in-app model panel and `docs/SCIENCE.md` for model assumptions and numerical validation. It does not implement operational ephemerides, finite-burn propulsion, parking-orbit escape, orbital capture, SOI event transitions, measured planetary DEM/raster ingestion, real hazard assessment, global coverage/revisit optimization, adaptive Runge–Kutta, multirevolution Lambert, general multi-leg mission design, position/model-parameter uncertainty, or cloud-synced notebooks. These are future extensions, not hidden or simulated capabilities.

## Source and release discipline

Canonical repository: `Debbybaks/AZIMUTH`, production branch `main`.
Vercel project: `azimuth` (`prj_p6HAk64JwTeXiupVM75QLhuiIOdN`) in `debbybaks-projects`.
Every release must pass tests and build, commit all source, push without rewriting history, verify remote files/head, then deploy that Git commit. Generated `dist/`, credentials, and dependency caches stay untracked. Commits use Deborah Bakare's configured Git identity and contain no assistant co-author attribution.

## Asset acknowledgments

Barlow and IBM Plex Mono fonts: Google Fonts, SIL Open Font License (copies under `public/fonts/`). Earth visual texture: Three.js example `earth_atmos_2048.jpg`, from the Three.js r160 example assets (https://github.com/mrdoob/three.js/tree/r160/examples/textures/planets). Other body textures are procedural illustrations, not scientific surface datasets. Icons: Lucide (ISC). Renderer: Three.js (MIT).
