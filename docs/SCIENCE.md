# Scientific model and verification

## Units and frames

Dynamics use AU and days internally: 1 AU = 149597870700 m and 1 day = 86400 s. Solar μ = 0.0002959122082855911 AU³/day². UI labels state where AU or km/s are used. CSV uses SI. JSON declares units. Solar coordinates are heliocentric mean J2000 ecliptic. The three-body experiment is barycentric. Earth ground tracks use a spherical Earth of mean radius 6371 km and sidereal rotation period 86164 s.

## Planetary model

Fixed J2000 approximate Kepler elements, with mean anomaly advanced using n = √(μ/a³). Planet elements are adapted from NASA/JPL's approximate positions model: https://ssd.jpl.nasa.gov/planets/approx_pos.html . Secular element rates and outer-planet correction terms are not implemented; these states must not be interpreted as precise ephemerides, particularly far from J2000. The Moon is an approximate Earth-relative Kepler ellipse superposed on Earth’s heliocentric state.

## Transfer and propagation

Universal-variable Lambert solution: Stumpff C(z), S(z), bracketed bisection in the zero-revolution domain, prograde geometry, maximum 160 iterations, explicit residual in days. Near collinear geometries reject rather than fabricating a result. Target endpoints are calculated from the approximate body states at departure and arrival epochs.

Numerical integration supports explicit Euler, improved Euler (Heun), RK4, and velocity Verlet. All use fixed steps; no adaptive tolerance is claimed. Restricted N-body propagation adds planetary accelerations and the heliocentric indirect term. Planet states are sampled at the current full-step epoch, not each internal RK stage. Forces within two planetary radii of the departure body are excluded to avoid a point singularity; this is not a launch/escape model. The fully interacting laboratory instead advances all three point masses together using velocity Verlet, with close-approach termination.

An encounter means final separation below 100000 km. It is not capture or insertion. Departure and arrival v∞ are relative velocities, not parking-orbit maneuver budgets. The actual trajectory follows the integrator; perturbations and numerical error can create a miss. Closest approach is sampled at integration steps, not a root-refined event. Sampling may underestimate temporal precision.

## Other laboratories

Hohmann: circular, coplanar, ideally phased heliocentric orbits with two instantaneous burns.
Gravity assist: Earth patched conic with fixed 29.78 km/s planet velocity, no atmosphere; equal planet-relative excess speeds, rotated by 2 asin(1/e).
Coordinate transforms: fixed J2000 obliquity 23.43928°, user-provided latitude and local sidereal angle; no precession, nutation, refraction, or observer parallax.
Ground tracks: circular two-body Earth orbit, spherical rotating Earth, zero-elevation geometric footprints. No geopotential or repeat-groundtrack scheduling.
Decay: fixed exponential atmosphere anchored at 700 km, density 3.614e-13 kg/m³, scale height 88.667 km; instantaneous rate only, not lifetime prediction.
GIS: reference NASA mission and JPL DSN site coordinates are approximate. There is no measured terrain raster. Synthetic h = 1500 sin(7φ) cos(4λ) m is explicitly marked. Great-circle distances use a sphere. This is not an actual landing-site safety assessment.
Optimization: finite bounded grid search on departure plus arrival excess-speed cost; no global convergence guarantee.
Monte Carlo: seeded LCG and Box–Muller Gaussian velocity perturbations only. Mean, empirical quantile, and empirical encounter fraction are computed from actual propagations. Small ensembles do not establish precise probabilities.

## Automated validation

`npm test` checks Hohmann Earth–Mars benchmarks (~258.87 d, ~5.594 km/s), elements/state round trips, Lambert endpoint recovery, intentional miss behavior, RK4 refinement, Verlet energy behavior, obliquity/horizon transforms, geostationary period, gravity-assist constraints, N-body barycenter conservation, restricted perturbations, and input bounds.

Browser QA covers WebGL initialization, calculation, launch, replay scrub, all labs, worker search and ensembles, experiment save, exports, and responsive layout. No test validates these simplified models against flight-grade SPICE ephemerides. Production accuracy is intentionally not claimed.
