# One-dimensional exhaust gas dynamics (`exhaust-1d`)

This branch replaces the lumped-volume exhaust of engine-sim with a
one-dimensional, compressible, unsteady gas-dynamics model of the header
and collector. Header geometry (primary length and diameter, collector
diameter and length, a megaphone cone, a tailpipe) now changes both the
engine's breathing and the sound. The old model is kept behind
`legacy_exhaust: true`.

Files: `include/pipe_1d.h`, `src/pipe_1d.cpp` (duct solver),
`include/exhaust_network.h`, `src/exhaust_network.cpp` (primaries, junction,
collector, open end, port coupling), changes in `exhaust_system.*`,
`combustion_chamber.*`, `piston_engine_simulator.*`,
`scripting/include/exhaust_system_node.h`, `es/objects/objects.mr`, tests in
`test/pipe_1d_tests.cpp`, `test/exhaust_network_tests.cpp`,
`test/exhaust_engine_tests.cpp`.

## 1. What the legacy model did

Each `exhaust_system` was one `GasSystem` volume (`collector_cross_section_area
x length`). Every cylinder had a lumped "runner + primary" volume
(`exhaust_runner_volume + exhaust_runner_cross_section_area x
(primary_tube_length + primary_length)`) that flowed through an orifice
(`primary_flow_rate`) into the collector, which flowed to the atmosphere
through `outlet_flow_rate`. Nothing propagated; `length` only set the collector
volume and an audio delay. The sound was each cylinder's runner gauge
pressure, delayed by `(primary_length + length) / 343 m/s`, scaled by
`audio_volume / (cylinders x L^2)`, summed per exhaust system and then
convolved with the impulse-response WAV in the synthesizer.

## 2. The 1-D model

### 2.1 Equations

Quasi-1-D Euler equations for an ideal gas in a duct of slowly varying area
`A(x)`, with wall friction and wall heat transfer as source terms. Per cell
the conserved variables are `q1 = rho A`, `q2 = rho u A`, `q3 = E A`,
`E = p/(gamma-1) + rho u^2/2`, plus two passive scalars `rho A x_k` carrying
the fuel and the O2 fraction of the mixture (so the exhaust O2 gauge keeps
working):

```
d(rho A)/dt     + d(rho u A)/dx           = 0
d(rho u A)/dt   + d((rho u^2 + p) A)/dx   = p dA/dx - F_w
d(E A)/dt       + d(u (E + p) A)/dx       = Q_w
F_w = 2 f rho |u| u A / D            (Fanning friction, tau_w = f rho u^2 / 2)
Q_w = 2 f rho |u| A / D  cp (T_w - T) (Reynolds analogy, St = f/2)
```

Gas model: the same ideal gas as `GasSystem` (5 degrees of freedom, so
`gamma = 1 + 2/5 = 1.4`, `R = R_u / M_air = 287.0 J/(kg K)` with
`units::AirMolecularMass`). Using the same gas makes the energy exchange with
the port volumes exact (section 2.4).

Constants: Fanning `f = 0.005` and wall temperature `T_w = 700 K` are the
values used in the ExhaustTuner `wsFactory()` reference this code is ported
from; both are script parameters (`wall_temperature`; the friction factor is
in `Pipe1D::Parameters`). `f = 0.005` is in the range Blair uses for exhaust
pipes (Design and Simulation of Four-Stroke Engines, ch. 2); 700 K is a
plausible steady header skin temperature, kept fixed (no wall thermal
inertia).

### 2.2 Scheme

MUSCL-Hancock finite-volume scheme (Toro, *Riemann Solvers and Numerical
Methods for Fluid Dynamics*, 3rd ed., ch. 14): minmod-limited slopes of the
primitive variables, a half-step predictor with the linearised primitive
equations, then an HLLC Riemann flux (Toro ch. 10.4-10.5) at every face,
and a conservative update with the sources. Passive scalars travel with the
mass flux of the side that holds the contact wave. Density is floored at
1e-6 of its reference. CFL number 0.8 on `|u| + c` over all cells of a
network; the network sub-steps inside each simulation step as often as
needed, so any `simulation_frequency` works (for the shipped engines, 1 to 6
sub-steps at 4 cm cells).

The scheme, the junction solvers and the engine coupling are a port of the
validated JavaScript reference (`wsFactory()` in the ExhaustTuner
calculator: Sod shock tube within 0.4%, linear-acoustic split at a step
within 0.3%). Two things were changed on purpose, see 2.3 (open end) and 2.4
(port coupling).

### 2.3 Boundaries

* **Port end** of a primary: the first cell receives (or gives up) the mass
  exchanged with the cylinder's port volume as a cell source; the face itself
  is a wall (section 2.4).
* **N-way junction** (primaries into the collector): method of
  characteristics, one pressure for all ends, the isentropic change of state
  along each incoming Riemann invariant, mass flows that sum to zero (Newton
  with analytic derivative, bisection fallback), inflows mixed by mass before
  they are sent down the branches that receive gas (reference `junctionN`,
  Blair's constant-pressure junction).
* **Step** between two collector pipes of different diameter: the two-pipe
  version of the same (reference `junction`). Cones are not steps: a linear
  diameter taper stays inside one pipe through `A(x)`.
* **Open end** to still air at 1 atm, 25 C: method of characteristics
  (Blair ch. 2). Outflow: boundary pressure = atmospheric, state reached
  isentropically from the end cell; a supersonic exit above atmospheric keeps
  the cell state. Inflow: the invariant of the (hot) pipe gas is matched to
  the cold inflowing air that expands isentropically from rest,
  `cp T_b + u_b^2/2 = cp T_atm`, `p_b = p_atm (T_b/T_atm)^(gamma/(gamma-1))`,
  i.e. the contact between hot pipe gas and cold air sits at the exit. The
  reference used a ghost state with the Bernoulli velocity `sqrt(2 cp (T_a -
  T_g))`; that velocity grows like `sqrt(dp)` and amplified round-off until a
  pipe at rest with one closed end started to flow (test `PipeAtRestStaysAtRest`
  caught it). The characteristic form gives the acoustic inflow
  `u ~ dp / (rho c)` for small suction and the same Bernoulli limit for
  strong suction.

### 2.4 Coupling to the cylinder

The cylinder side keeps its lumped exhaust-port `GasSystem`, now holding only
the head's `exhaust_runner_volume` (floored at one cell's volume), and the
existing cylinder -> port valve flow is untouched. In place of the legacy
port -> collector orifice, each fluid sub-step calls
`ExhaustSystem::exchangePort()`:

* compressible orifice flow between the port (stagnation `p`, `T`) and the
  pipe's first cell (stagnation state from its `p`, `T`, `u`), in the
  direction of the higher stagnation pressure (reference `valveFlow`), with
  effective area `Cd A = min(k sqrt(M_air), A_port, A_pipe)`. The identity
  `Cd A = k sqrt(M)` comes from `GasSystem::flowRate` (`n_dot = k p0
  sqrt(gamma/(R T0)) f(pr)`) against the physical orifice (`m_dot = Cd A p0
  sqrt(gamma M/(R T0)) f(pr)`), so `primary_flow_rate` keeps its meaning.
* gas leaving the port takes its internal energy plus the flow work `R T`
  per mole (the port loses enthalpy); the pipe cell receives the enthalpy and
  the jet momentum `m u_jet` with `u_jet` the isentropic velocity from the
  port's stagnation state to the cell pressure (capped at sonic). Gas
  entering the port brings its stagnation enthalpy. Total energy and mass
  are therefore conserved exactly (test `PortExchangeConservesEnergy`).
* clamps: at most 90% of the port's moles per sub-step (as
  `GasSystem::flow`) and at most 20% of the cell's mass per sub-step in the
  reverse direction.

Both sides are updated at once at the fluid sub-step rate (8 x
`simulation_frequency`), and the pipes advance once per simulation step.
The first version averaged the port flux over a step and fed it to the pipe
as a face flux; with a 50 cc port against a 40-60 cc first cell that
explicit coupling oscillated and could empty the first cell. Updating the
cell directly each sub-step (as the reference does for its chambers) is
stable at all tested cell sizes and simulation frequencies.

### 2.5 Sound

The far-field pressure of an open pipe end is a monopole,
`p(r, t) = rho_0 / (4 pi r) dQ/dt`, `Q = m_dot / rho_0` the volume velocity
through the exit. `ExhaustNetwork::radiatedPressure()` evaluates it at
`r = 1 m` from the step-averaged exit mass flow. `writeToSynthesizer()`
feeds `attenuation^3 x audio_volume x 1e5 x p(1 m)` per exhaust system (the
same rpm attenuation as before). The gain `1e5` has no physical meaning: it
brings the signal to the same order as the legacy `1600 x gauge pressure x
audio_volume / L^2` so that the synthesizer's leveler, `audio_volume` and the
impulse responses keep working. The per-cylinder delay filters are no longer
used for 1-D systems; `sound_attenuation` (a per-cylinder sound hack) has no
effect on them. Legacy systems keep the old path unchanged.

## 3. Script parameters

All new inputs are optional; existing scripts load unchanged and get a 1-D
network built from their legacy parameters.

| input | default | meaning |
|---|---|---|
| `legacy_exhaust` | `false` | `true`: lumped model and per-cylinder audio delay exactly as before |
| `primary_diameter` | 0 = derived | diameter of every primary; derived from the head's `exhaust_runner_cross_section_area` |
| `collector_diameter` | 0 = derived | derived from `collector_cross_section_area` |
| `collector_length` | -1 = `length` | collector pipe from the junction to the cone/tailpipe/exit |
| `megaphone_length`, `megaphone_exit_diameter` | 0 = none | cone after the collector |
| `tailpipe_length`, `tailpipe_diameter` | 0 = none / same diameter | straight pipe after the cone |
| `pipe_cell_size` | 4 cm | target cell size; a pipe has at least 6 and at most 200 cells |
| `wall_temperature` | 700 K | wall temperature for the heat-transfer source |

Per-cylinder primary length = `primary_tube_length` + that cylinder's
`primary_length` (both existed before). `outlet_flow_rate` and
`velocity_decay` are legacy-only.

Why the new model is the default: the purpose of the branch is that header
geometry changes the sound and the torque, old files get a sensible network
without edits, and `legacy_exhaust: true` is a one-line way back for A/B
listening and regression checks. Two consequences for old files to know
about:

* `length` was often a volume fudge (the `exhaust_system_parameters` default
  is `volume / collector_cross_section_area` = 100 L / 81 cm^2 = 12.3 m). A
  12 m collector is now a 12 m pipe. It still works (the 200-cell cap makes
  its cells 6 cm) but costs time; set `collector_length` to the real value.
* torque curves change (section 6): 4-1 headers now show their tuning
  effects, including the low-rpm hole.

## 4. Validation

Linux, clang 18, headless build of the core library and the script
interpreter (`test/` with gtest; the engine-level tests need
`engine-sim-script-test`, which runs about 8 minutes).

Unit and network tests (`engine-sim-test`):

| test | result |
|---|---|
| (a) closed-open pipe at rest, 2000 steps | pressure, velocity, temperature unchanged to 1e-6 / 1e-9; mass and energy to 1e-9 |
| (b) +2 kPa pulse at the closed end of a 1 m closed-open pipe | returns inverted (-1755 Pa at the wall) after 5.740 ms vs 2L/c = 5.778 ms (-0.66%) |
| (c) quarter-wave resonance, 0.8 m closed-open pipe | 108.0 Hz vs c/4L = 108.16 Hz (-0.15%) |
| (d) mass conservation, 4 pulsed ports into a 4-1 for 1 s | imbalance 1e-10 % of the throughput (the engine runs below give 1e-13 relative) |
| (e) 3-way junction, pulse from one branch | mass conserved to 1e-9; transmitted amplitude 0.455 vs linear acoustics 2A1/(A1+A2+A3) = 0.471 (-3.4%); the two receiving branches agree to 0.004% |
| orifice | choked mass flow matches the isentropic formula to 4e-6 (the gas constant rounding) |
| port exchange | mass to 1e-9, total energy to 1e-6 over a 2 ms blowdown |

Engine-level tests (`engine-sim-script-test`):

* (f) all 25 shipped engines (`assets/engines/*/*.mr`; the five older files
  through `set_engine(...)`) compile, build their networks and run 2 s of
  cranking plus 1.5 s on the dyno at 2500 rpm without NaN; 23 start within the
  2 s crank and make positive torque (30 to 620 N m). The radial-9 and the
  Merlin V12 need a longer crank than the test gives them, in both models
  (the radial-9 starts with 6 s of cranking and makes 778 N m at 2000 rpm;
  the Merlin did not start in 8 s with either model in this harness).
* (g) `legacy_exhaust: true` on the EJ25: no network is built, the runner
  volume equals the legacy formula, the engine runs. Bit-identity with the
  pre-change code was checked outside the test suite with the same
  comparison program built against the pre-change tree (commit 36ad6a2) and
  against this branch: 3 s runs of the EJ25 (hold 3500 rpm) and the GM LS
  (hold 4000 rpm) printed rpm, torque, cylinder and runner pressures and the
  cumulative synthesizer input every 0.1 s with 17 significant digits;
  the two outputs are identical.
* physics: see section 6.

Pre-existing: upstream's own `GasSystemTests.PressureEquilibriumMaxFlow*`
and `FunctionTests.FunctionGaussianTest` fail, and
`SynthesizerTests.SynthesizerSystemTestSingleThread` hangs, identically on
the untouched tree. Not caused by this branch.

## 5. Performance

Test machine: 2-core Xeon at 2.8 GHz (shared; CPU time of the process is
reported, wall time was up to 3x higher under load). GM LS V8, two 4-1
systems, 0.74 m primaries, 2.5 m and 4.4 m collectors, `simulation_frequency`
10000, dyno hold 4000 rpm, 3 s, no audio:

| model | cells | sub-steps / step | CPU s per simulated s | relative |
|---|---|---|---|---|
| legacy | - | - | 0.87 | 1.0 |
| 1-D, 4 cm (default) | 340 | 3.6 | 1.85 | 2.1 |
| 1-D, 3 cm | 456 | 4.5 | ~2.7 | ~3.1 |
| 1-D, 2 cm | 663 | 6.6 | ~4.9 | ~5.6 |

Per cell update the solver costs about 65 ns in isolation and 80 ns inside
the engine (cache pressure, junctions). The cost scales with `cells x
sub-steps`, i.e. with `1 / cell_size^2`, which is why 2 cm (the first
target) was about 5x the legacy cost and 4 cm was chosen as the default:
it keeps a V8 near 2x legacy. Trade-off: the shortest wavelength carried
cleanly is about 8-10 cells, so 4 cm resolves up to ~1-1.5 kHz in hot gas
(c ~ 600-800 m/s) and 2 cm up to ~2-3 kHz; the firing harmonics of every
shipped engine are well inside the 4 cm range, the hiss and crack above
2 kHz are left to the impulse response. Torque changed by 3-4% between 2 cm
and 4 cm on the V8 (678 vs 703 N m at 4000 rpm). Other shipped engines at
their own `simulation_frequency`: 0.7-2.2 CPU s per simulated second with
the 1-D model (the radial-9 with two default 12 m collectors, 978 cells,
is the heaviest).

Things that were tried and kept: one half-step predictor per cell (not per
face), no divisions in the hot loop, analytic Newton derivatives and one
`pow` per branch in the junctions, cached orifice constants, the port jet
speed removed from the CFL estimate (it is a cell source, not a face flux:
15% fewer sub-steps, same results). Not done: single precision (the mass
accumulators would drift over long runs), SIMD across cells (the HLLC
branches need masking; large effort, also for MSVC), threading one exhaust
system per thread (synchronisation 10000 times a second). On a slower
machine use `pipe_cell_size: 6 * units.cm` or `legacy_exhaust: true`.

## 6. Physics check: primary length

2.0 L inline-4, 4-1 header, 721 mm x 44 mm primaries, 488 mm x 79.5 mm
collector (`test/data/exhaust_tuner_i4.mr`, rev limit raised to 9000 for the
7800 rpm point), wide open, dyno hold, mean torque over the last 2 s of a 4 s
run. Primaries scaled by 0.9 / 1.0 / 1.1. Cells 4 cm.

4 cm cells (default):

| model | primaries | 4000 rpm | 5000 rpm | 6000 rpm | 7000 rpm | 7800 rpm |
|---|---|---|---|---|---|---|
| legacy | 649 mm (x0.9) | 143.8 | 173.7 | 193.1 | 184.0 | 166.0 |
| legacy | 721 mm (x1.0) | 145.6 | 178.1 | 190.7 | 183.0 | 163.0 |
| legacy | 793 mm (x1.1) | 148.8 | 178.3 | 196.3 | 180.9 | 160.3 |
| 1-D | 649 mm (x0.9) | 54.6 | 138.0 | 180.7 | 183.1 | 186.6 |
| 1-D | 721 mm (x1.0) | 61.1 | 153.2 | 148.5 | 191.5 | 187.7 |
| 1-D | 793 mm (x1.1) | 97.4 | 133.6 | 157.2 | 197.4 | 187.7 |

2 cm cells:

| model | primaries | 4000 rpm | 5000 rpm | 6000 rpm | 7000 rpm | 7800 rpm |
|---|---|---|---|---|---|---|
| 1-D | 649 mm (x0.9) | 66.7 | 128.2 | 170.6 | 183.5 | 163.5 |
| 1-D | 721 mm (x1.0) | 69.6 | 141.8 | 183.5 | 183.1 | 163.5 |
| 1-D | 793 mm (x1.1) | 92.7 | 165.9 | 185.6 | 183.9 | 162.4 |

(torque in N m; the legacy rows do not depend on the cell size)

The legacy model changes by 1-3% (the longer primary is a bigger lumped
volume); the 1-D model changes by 10-40% at 4000-6000 rpm, in both grids.
The low-rpm hole (4000-5000 rpm) is the 4-1 interference: the neighbour's
blowdown pulse, transmitted through the junction and doubled at the closed
valve, arrives at the port during valve overlap at 1.5-1.6 bar and pushes
exhaust gas into the cylinder and the intake runner. The effect is real for
4-1 headers but its size here is large; see the limits below.

The two grids disagree above 6000 rpm (e.g. 721 mm at 6000 rpm: 148 N m at
4 cm vs 184 N m at 2 cm; at 7800 rpm 188 vs 164). The 2 cm results are the
smoother ones (torque rises monotonically with primary length at 5000-6000
rpm and the three lengths converge at 7800 rpm, as expected when the tuned
pulse falls outside the overlap window); at 4 cm the blowdown pulse front
(a few cm wide at these speeds) is under-resolved and its arrival time at
the port shifts. Where the exact torque matters, run with
`pipe_cell_size: 2 * units.cm` (section 5 gives the cost); the sound and
the existence of the tuning effect are the same at both sizes. The V8 at
4000 rpm (where the pulses are longer) moved only 3-4% between the grids.

## 7. Sound

WAV previews rendered through the harness (same synthesizer settings as the
app, impulse responses from `es/sound-library`), held rpm, last 1.5 s
analysed with a Hann window:

| run | RMS | clipped samples | strongest lines (Hz : relative) | energy 150-400 / 400-1000 / 1000-3000 Hz |
|---|---|---|---|---|
| I4 6000 rpm, legacy | 27559 | 54% | 200:1.00 600:0.15 400:0.07 800:0.06 | 97% / 3% / 0.1% |
| I4 6000 rpm, 1-D 4 cm | 12945 | 2.5% | 400:1.00 200:0.77 800:0.57 1600:0.31 600:0.29 1200:0.27 | 59% / 29% / 11% |
| V8 (LS) 4000 rpm, legacy | 17530 | 3.0% | 433:1.00 267:0.36 400:0.25 167:0.23 | 20% / 78% / 0.4% |
| V8 (LS) 4000 rpm, 1-D 4 cm | 14943 | 2.1% | 433:1.00 933:0.60 800:0.56 700:0.49 633:0.46 | 7% / 76% / 17% |

The firing fundamental (200 Hz for the I4 at 6000 rpm, 267 Hz for the V8 at
4000) and its harmonics come through in the 1-D model; the legacy I4 is a
nearly pure 200 Hz tone that the leveler clips more than half the time,
while the 1-D signal has 10-20% of its energy above 1 kHz (blowdown pulse
shape) and clips little. Files: `i4_legacy_6000.wav`, `i4_1d_6000.wav`,
`ls_legacy_4000.wav`, `ls_1d_4000.wav` in the harness output (not committed).

## 8. Known limits

* Pulse amplitudes are probably on the strong side: no junction loss
  coefficients (constant-pressure junction), a fixed wall temperature, no
  pulse attenuation other than friction; and the intake side is still a
  lumped runner, so exhaust gas pushed into it stays for the whole intake
  stroke instead of being swept out by an intake wave. Together these make
  the 4-1 low-rpm hole deeper than a dyno would show.
* The port -> primary orifice (`primary_flow_rate`) is inherited from the
  legacy structure. Physically the port and primary are one duct; the
  orifice keeps the old parameter meaningful but adds a sudden-expansion
  loss at the pipe inlet.
* Cones are limited to a linear taper inside one pipe; steps need two
  pipes; no X/H crossovers between exhaust systems, no mufflers or chambers
  (the impulse response still does that), no catalyst.
* `sound_attenuation` has no effect on 1-D systems.
* Primaries are straight pipes of one diameter; the head's port volume is
  lumped.
* Cell size trades CPU for bandwidth (section 5). The 200-cell cap coarsens
  very long pipes silently.
* The radiated sound is a monopole at 1 m; no directivity, no ground.
