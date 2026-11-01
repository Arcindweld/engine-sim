# One-dimensional exhaust gas dynamics (`exhaust-1d`)

This branch adds a one-dimensional, compressible, unsteady gas-dynamics
model of the header and collector next to the lumped-volume exhaust of
engine-sim. Header geometry (primary length and diameter, collector diameter
and length, a megaphone cone, a tailpipe) then changes both the engine's
breathing and the sound. Selection is per exhaust system with
`exhaust_model`: `"auto"` (the default) runs the 1-D model only for files
that give 1-D geometry, so the shipped engines keep their behaviour
bit-for-bit and calculator exports get the new model; `"1d"` and `"legacy"`
force it.

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
network. The port boundaries are solved after every fluid sub-step of the
engine (8 per simulation step); the pipes advance at their own CFL step. When
the fluid sub-step is the shorter of the two (4 and 2 cm cells at 20 kHz:
6 us against about 20-45 us), the pipes are held while the ports run on
against them, the face fluxes the ports give are summed, and the pipes take
one step over the time the ports have run ahead, with the averaged fluxes,
before a sub-step would put the ports more than one CFL step ahead of them
and at the end of every simulation step. When the CFL step is the shorter one
(from about 1 cm down at 20 kHz), ports and pipes advance together in CFL
steps inside each fluid sub-step (section 2.4).

The scheme, the junction solvers and the engine coupling are a port of the
validated JavaScript reference (`wsFactory()` in the ExhaustTuner
calculator: Sod shock tube within 0.4%, linear-acoustic split at a step
within 0.3%; this C++ solver is checked against the exact Riemann solution
itself, section 4). Three things were changed on purpose, see 2.3 (open end,
the receiving side of junctions) and 2.4 (port coupling).

### 2.3 Boundaries

* **Port end** of a primary: a valve-type characteristic boundary between
  the cylinder's port volume and the pipe's inlet face, solved every fluid
  sub-step (section 2.4).
* **N-way junction** (primaries into the collector): method of
  characteristics, one pressure for all ends, the isentropic change of state
  along each incoming Riemann invariant, mass flows that sum to zero, inflows
  mixed by mass before they are sent down the branches that receive gas
  (reference `junctionN`, Blair's constant-pressure junction). A receiving
  end takes the mixed inflow's state at the face, as the open end takes the
  outside air: velocity from the pipe's own characteristic at the junction
  pressure, temperature from the mixed stagnation enthalpy, density `p / (R
  T)`. (The reference, and the first version here, gave it the receiving
  pipe's own density; the mass flux then disagreed with the energy flux and
  every temperature front put a pressure step at the junction, 6.5-10 kPa for
  900 K gas meeting 300 K gas at 50 m/s, test
  `TemperatureFrontCrossesJunctionsUndisturbed`.) Newton with an analytic
  derivative with the ends' own densities starts the secant solve of the
  mass balance with the mixed state (bisection fallback); the receiving ends
  are then scaled to the inflow so the fluxes sum to zero to round-off. Each
  evaluation takes one `log` of the junction pressure and one `exp` per end
  (`(p_j / p_k)^((gamma - 1) / (2 gamma))`), the density ratio being that to
  the integer power `dof`.
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
the head's `exhaust_runner_volume` (floored at 2 cm of the runner's
cross-section), and the existing cylinder -> port valve flow is untouched. In
place of the legacy port -> collector orifice, the port volume feeds its
primary through a boundary condition at the primary's inlet face
(`ExhaustNetwork::portBoundary`), the valve boundary of the method of
characteristics (Benson, *The Thermodynamics and Gas Dynamics of
Internal-Combustion Engines*, vol. 1, Oxford 1982; Winterbone & Pearson,
*Theory of Engine Manifold Design*, 2000):

* After every fluid sub-step it finds the face pressure `p_f` at which the
  compressible orifice flow between the port and the face equals the flow
  the pipe takes along its incoming characteristic, `u_f = J + 2 c_f /
  (gamma - 1)` with `J = u - 2 c / (gamma - 1)` of the first cell and `c_f`
  isentropic from it:
  * port -> pipe (port pressure above the face's closed-end pressure): the
    gas at the face is port gas, `T_f = (h0_port - u_f^2/2) / cp`, `rho_f =
    p_f / (R T_f)` (the contact between port gas and pipe gas sits at the
    face, as at the open end), the orifice runs from the port's `p`, `T` to
    `p_f`, and `u_f` is capped at sonic;
  * pipe -> port: the gas at the face is pipe gas, isentropic from the first
    cell to `p_f`; the orifice runs from its stagnation state to the port
    pressure, and an outflow is capped at sonic.

  The root is bracketed between the closed-end pressure and the port
  pressure. The unknown is `s = (p_f / p_c)^((gamma - 1) / (2 gamma))`: with
  `gamma = 1 + 2 / dof` (the gas model's degrees of freedom), `u_f = J + dof
  c_c s`, the first cell's isentrope gives `p_f = p_c s^(dof + 2)`, `rho_f =
  rho_c s^dof` and `T_f = T_c s^2`, and the orifice function becomes `psi =
  t^dof sqrt(2 gamma / (gamma - 1) (1 - t^2))` with `t = (p_d / p_0)^((gamma -
  1) / (2 gamma))`, a product of `s` and constants. The search needs no
  `pow` (two per face to set it up): port -> pipe by Newton with the analytic
  derivative, started from the previous sub-step's root and kept inside the
  bracket; pipe -> port by secant steps from that root, with the Illinois
  method on the bracket as the fallback. Searching in `p_f` instead took
  about 15 `pow` calls per face; with two faces per port per sub-step that
  was a quarter of the run time of a V8 (section 5). For a full-area opening
  the boundary reduces to the usual isentropic inflow from a plenum and, for
  outflow, to a face at the port pressure.
* The face flux (mass, momentum with `p_f A`, stagnation enthalpy, species)
  goes to the pipe, and the port volume gives up or takes the same mass,
  total energy and species over the sub-step. Gas leaving the port takes its
  internal energy, the flow work `R T` per mole and its share of the port's
  bulk kinetic energy (the port's momentum is scaled with the gas that
  leaves, as `GasSystem::flow` does); gas entering it brings its stagnation
  enthalpy, and the bulk motion it slows becomes heat. Mass and total energy
  are conserved to round-off (test `PortExchangeConservesEnergy`, with the
  port at rest and moving, emptying and filling).
* The port volume is integrated with the explicit midpoint rule: the flux is
  taken from the port's state half a sub-step on. One explicit step made the
  port relax about 7% too fast at 20 kHz, a first-order error that left the
  torque drifting with the frequency (134.1 / 135.3 / 140.6 N m at 15 / 20 /
  40 kHz on 2 cm cells for the test I4 at 6000 rpm; 137.8 / 138.5 / 139.5
  with the midpoint rule).
* Effective area `Cd A = min(k sqrt(M_air), A_port, A_pipe)`. The identity
  `Cd A = k sqrt(M)` comes from `GasSystem::flowRate` (`n_dot = k p0
  sqrt(gamma/(R T0)) f(pr)`) against the physical orifice (`m_dot = Cd A p0
  sqrt(gamma M/(R T0)) f(pr)`), so `primary_flow_rate` keeps its meaning.
  At most 90% of the port's moles leave in one sub-step (as `GasSystem::flow`).
* Timing: the port boundaries are solved after every fluid sub-step of the
  engine (8 x `simulation_frequency`), right after the cylinders have flowed
  into their ports, against the pipes' held state; the pipes advance at their
  own CFL step with the face fluxes the ports gave them meanwhile (section
  2.2; `ExhaustNetwork::Parameters::holdPipes`, on by default), and at the end
  of every simulation step, where the outlet flow (the sound source) is
  closed. The port and the pipe exchange exactly the same mass, energy and
  species; until the pipes advance, that gas is counted as in transit at the
  inlet face. Holding the pipes is the method of characteristics' own
  assumption over a step (the incoming characteristic held), and the hold
  never exceeds one CFL step, so it is not the first form below, which
  averaged an explicit exchange over a whole simulation step.

For a small pressure difference the flow is then limited by the pipe's wave
impedance, `m ~ A dp / c`, and the port relaxes on `V / (c A)`, about
0.05 ms for a 50 cc port on a 44 mm primary, against a sub-step of 6 us at
20 kHz. Two earlier forms did not converge:

1. Averaging the port flux over a simulation step and feeding it to the pipe
   as a face flux: unstable with a 50 cc port against a 40-60 cc first cell.
2. The orifice flow computed from the first cell's own state at every fluid
   sub-step and added to the cell as a source, with the pipes advanced once
   per simulation step. An explicit square-root orifice law between two small
   volumes overshoots: its slope is unbounded at zero pressure difference, so
   the exchange settles into a two-step oscillation of about `rho (Cd A dt c^2
   (1/V_port + 1/V_cell))^2 / 2` (kPa at 20 kHz) and the flow reversed on
   nearly every sub-step; with the pipes frozen through eight exchanges the
   I4 at 6000 rpm made 148.7 / 183.6 / 118.5 / -6.7 N m on 4 / 2 / 1 / 0.5 cm
   cells. Advancing the pipes after every exchange removed the freezing but
   not the oscillation (131.8 / 134.5 / 141.1 / 143.4 N m). Section 6 has the
   numbers with the boundary form.

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
| `exhaust_model` | `"auto"` | `"auto"`: 1-D when any of the geometry inputs below (`primary_diameter` ... `pipe_cell_size`) is given, legacy otherwise; `"1d"` / `"legacy"`: force |
| `legacy_exhaust` | `false` | `true`: same as `exhaust_model: "legacy"` (lumped model and per-cylinder audio delay exactly as before) |
| `primary_diameter` | 0 = derived | diameter of every primary; derived from the head's `exhaust_runner_cross_section_area` |
| `collector_diameter` | 0 = derived | derived from `collector_cross_section_area` |
| `collector_length` | -1 = `length` | collector pipe from the junction to the cone/tailpipe/exit |
| `megaphone_length`, `megaphone_exit_diameter` | 0 = none | cone after the collector |
| `tailpipe_length`, `tailpipe_diameter` | 0 = none / same diameter | straight pipe after the cone |
| `pipe_cell_size` | 0 = 4 cm | target cell size; a pipe has at least 6 and at most 200 cells (giving it, even 4 cm, selects the 1-D model in `"auto"`) |
| `wall_temperature` | 700 K | wall temperature for the heat-transfer source |

Per-cylinder primary length = `primary_tube_length` + that cylinder's
`primary_length` (both existed before). `outlet_flow_rate` is legacy-only.
`velocity_decay` acts in both models: besides the legacy collector volume it
sets how fast the bulk flow of each cylinder's exhaust-port volume decays
(`CombustionChamber::flow`), and that port feeds the 1-D primary.

Why `"auto"`: the shipped engines and every existing user file keep running
exactly as before (section 4, test g), while a file that describes its
header in 1-D terms (the calculator exports do) gets the new model without a
switch. `exhaust_model: "1d"` turns it on for an old file; it then builds its
network from the legacy parameters. Two things to know when doing that:

* `length` was often a volume fudge (the `exhaust_system_parameters` default
  is `volume / collector_cross_section_area` = 100 L / 81 cm^2 = 12.3 m). A
  12 m collector becomes a 12 m pipe. It still works (the 200-cell cap makes
  its cells 6 cm) but costs time; set `collector_length` to the real value.
* torque curves change (section 6): 4-1 headers now show their tuning
  effects, including the low-rpm hole.

## 4. Validation

Linux, clang 18, headless build of the core library and the script
interpreter (`test/` with gtest; the engine-level tests need
`engine-sim-script-test`, which runs about 5 minutes). The Windows CI builds
`engine-sim-test` and runs the 1-D unit and network tests, and runs the test
engine of `test/data` on the 1-D model through the headless runner.

Unit and network tests (`engine-sim-test`):

| test | result |
|---|---|
| (a) closed-open pipe at rest, 2000 steps | pressure, velocity, temperature unchanged to 1e-6 / 1e-9; mass and energy to 1e-9 |
| (b) +2 kPa pulse at the closed end of a 1 m closed-open pipe | returns inverted (-1755 Pa at the wall) after 5.740 ms vs 2L/c = 5.778 ms (-0.66%) |
| (c) quarter-wave resonance, 0.8 m closed-open pipe | 108.0 Hz vs c/4L = 108.16 Hz (-0.15%) |
| (d) mass conservation, 4 pulsed ports into a 4-1 for 1 s, driven like the engine | imbalance 5e-11 % of the throughput (round-off) |
| (e) 3-way junction, pulse from one branch | mass conserved to 1e-9; transmitted amplitude 0.455 vs linear acoustics 2A1/(A1+A2+A3) = 0.471 (-3.4%); the two receiving branches agree to 0.004% |
| orifice | choked and subsonic mass flow match the isentropic nozzle formulas to 1e-9, the two meet at the critical pressure ratio, the flow falls monotonically to zero as the pressures equalise |
| (h) shock tubes against the exact Riemann solution | Sod (1 / 0.125 kg/m3, 1e5 / 1e4 Pa): L1 density error 1.09% on 100 cells, 0.34% on 400, pressure 0.23%, star-region pressure 30305 Pa vs exact 30313 Pa; exhaust blowdown (5 bar 1200 K into 1 bar 300 K): 2.24% / 0.81%, 278517 vs 278552 Pa; mass to 1e-15 |
| (i) temperature front (900 K gas behind 300 K gas, same p and u = 50 m/s) through junctions | two-pipe step, two-ended and three-ended junction, both flow directions: pressure stays within 1e-4 Pa and velocity within 1e-6 m/s of uniform (the receiving-end density of the first version gave 6.5-10 kPa and 23-45 m/s) |
| port boundary, energy | mass and total energy (port internal + bulk kinetic, pipes) to 1e-16 relative over 1.5 ms: port at rest and moving at 38 / 20 m/s, emptying (3 atm to 1 atm) and filling (0.5 atm to 1 atm) |
| port boundary, steady flow | a port held at 1.3 bar, 900 K into a frictionless pipe open at its end: 0.0593928 kg/s through the boundary and 0.0593931 kg/s out of the pipe vs the orifice formula's 0.0593936 kg/s |
| coupled ports: 4 ports with a valve-like feed into a 4-1, driven like the engine | the port flow changes direction 23-28 times per cycle (the waves; with the cell-source exchange it reversed on 90% of the 3200 sub-steps per cycle). RMS port gauge pressure on 1 cm cells at 10 / 20 / 40 kHz: 6977 / 6975 / 6976 Pa (7135 / 7016 / 6997 with one explicit port step instead of the midpoint rule, before the pipes were held); at 20 kHz on 4 / 2 / 1 cm cells: 5748 / 6839 / 6975 Pa, 7030 on 0.5 cm (the port's own ring near c A / (2 pi V), 2 kHz here, is above what 4 cm cells carry) |
| held pipes against lock-step (`HeldPipesConvergeToTheFineGrid`) | the same coupled ports at 20 kHz. Held, the pipes take 2 / 3 steps per simulation step on 4 / 2 cm cells instead of 8 and land closer to the 1 cm result: RMS 5748 / 6839 Pa held, 5402 / 6422 lock-step, 6977 on 1 cm; forward flow 46.3 / 53.6 against 43.3 / 51.0 and 55.0 mg per cycle; back-flow 24.8 / 35.3% against 19.6 / 31.8% and 37.0%. Lock-step runs the pipes at a quarter to a third of their CFL step, where the upwind scheme is most diffusive. On 1 cm cells at 40 kHz, where both resolve the waves, held and lock-step agree: 6976 / 6991 Pa, 55.2 / 55.0 mg |

Engine-level tests (`engine-sim-script-test`):

* (f) all 25 shipped engines (`assets/engines/*/*.mr`; the five older files
  through `set_engine(...)`), forced to the 1-D model with
  `exhaust_model: "1d"`, compile, build their networks and run 2 s of
  cranking plus 1.5 s on the dyno at 2500 rpm without NaN; 23 start within the
  2 s crank and make positive torque (27 to 612 N m), with the mass account
  closed to 4e-13 of the throughput. The radial-9 and the
  Merlin V12 need a longer crank than the test gives them, in both models
  (the radial-9 starts with 6 s of cranking and makes 778 N m at 2000 rpm;
  the Merlin did not start in 8 s with either model in this harness).
* (g) `legacy_exhaust: true` on the EJ25: no network is built, the runner
  volume equals the legacy formula, the engine runs. `AutoModelSelection`:
  the EJ25 as shipped runs legacy, the calculator export in `test/data` runs
  1-D, `"1d"`, `"legacy"`, `legacy_exhaust: true` and a lone `pipe_cell_size`
  select as documented. Bit-identity with the pre-change code was checked
  outside the test suite with the same comparison program built against the
  pre-change tree (commit 36ad6a2) and against this branch: 3 s runs of the
  EJ25 (hold 3500 rpm) and the GM LS (hold 4000 rpm), both as shipped (auto)
  and with `legacy_exhaust: true`, printed rpm, torque, cylinder and runner
  pressures and the cumulative synthesizer input every 0.1 s with 17
  significant digits; the outputs are identical.
* (h) primary length (`PrimaryLengthTunesOnlyTheOneDModel`): the test I4
  held at 6000 rpm with 649 / 793 mm primaries, each torque the mean of three
  random-burn seeds: 4 cm 166.8 / 154.0 N m, 2 cm 161.9 / 152.3 N m (changes
  -12.8 and -9.6 N m, 7.7% and 5.9%), legacy 192.5 / 192.4 N m (0.07%). The
  test asks for more than 5% on both grids and under 3% legacy, the two grids
  within 8% and agreeing on the sign of the change and on its size within a
  half: one seed's change carries about 3 N m of noise, three seeds' about 2.
* physics: see section 6.
* Determinism: `Synthesizer::m_latency` was never initialized (upstream
  code), and `Simulator::startFrame` reads it before the audio side first
  sets it, so the first frame's step count came from whatever the heap held
  and runs inside one process depended on what had run before them (the
  primary-length test's torques moved by up to 2% between running it alone
  and after (f), enough to fail its 4 cm / 2 cm agreement check). It is now
  initialized; the engine-level tests give the same values alone and in the
  suite, valgrind memcheck finds no uninitialized reads in a legacy or a 1-D
  run, and fresh headless runs are unchanged (the test I4 at 6000 rpm:
  141.205 N m on the 1-D model, 190.688 N m legacy, before and after).

Pre-existing: upstream's own `GasSystemTests.PressureEquilibriumMaxFlow*`
and `FunctionTests.FunctionGaussianTest` fail, and
`SynthesizerTests.SynthesizerSystemTestSingleThread` hangs, identically on
the untouched tree; `SynthesizerTests.SynthesizerSystemTest` fails its sample
checks and then crashes, here and on the exhaust-1d head before this change
(the synthesizer is upstream's code apart from two includes). Not caused by
this branch.

## 5. Performance

Test machine: 2-core Xeon at 2.8 GHz, otherwise idle (CPU time of the
process is reported). Dyno hold, 4 s with 2 s settling, no audio, each
engine at its own `simulation_frequency`. "Lock-step" is the pipes advanced
with every fluid sub-step (`holdPipes = false`), "before" the build before
the `pow`-free port face search of section 2.4; both were run in the same
session:

| engine | sim. frequency | model | cells | pipe steps per simulation step (lock-step) | CPU s per simulated s (lock-step / before) | relative to legacy |
|---|---|---|---|---|---|---|
| GM LS V8 (two 4-1), 4000 rpm | 10 kHz | legacy | - | - | 0.72 | 1.0 |
| GM LS V8 (two 4-1), 4000 rpm | 10 kHz | 1-D, 4 cm | 340 | 4.1 (8) | 2.01 (2.81 / 3.77) | 2.8 |
| GM LS V8 (two 4-1), 4000 rpm | 10 kHz | 1-D, 3 cm | 456 | 6.3 | 2.83 | 3.9 |
| GM LS V8 (two 4-1), 4000 rpm | 10 kHz | 1-D, 2 cm | 663 | 8.2 (9.3) | 4.37 (4.91 / 5.89) | 6.1 |
| Subaru EJ25 (two systems), 3500 rpm | 20 kHz | legacy | - | - | 0.72 | 1.0 |
| Subaru EJ25 (two systems), 3500 rpm | 20 kHz | 1-D, 4 cm | 424 | 2.0 (8) | 2.01 (5.10 / 6.25) | 2.8 |
| Lexus LFA V10, 5000 rpm | 6.5 kHz | legacy | - | - | 0.73 | 1.0 |
| Lexus LFA V10, 5000 rpm | 6.5 kHz | 1-D, 4 cm | 454 | 7.4 (8) | 2.41 (2.58 / 3.19) | 3.3 |
| test inline-4, 6000 rpm | 20 kHz | legacy | - | - | 0.70 | 1.0 |
| test inline-4, 6000 rpm | 20 kHz | 1-D, 4 cm | 84 | 2.2 (8) | 1.54 (2.12 / 3.48) | 2.2 |
| test inline-4, 6000 rpm | 20 kHz | 1-D, 2 cm | 168 | 5.0 (8) | 2.35 (3.05 / 4.21) | 3.4 |
| test inline-4, 6000 rpm | 20 kHz | 1-D, 1 cm | 337 | 10.1 (10.8) | 5.49-5.67 (5.94 / 7.56) | 8 |

Torque on the held pipes, N m (lock-step): LS 721.2 / 720.5 / 719.3 on
4 / 3 / 2 cm (719.8 / 719.3 / 718.7), legacy 714.0; EJ25 199.7 (199.7),
legacy 207.1; LFA 150.5 (153.1), legacy 197.5; the inline-4 in section 6.
The LFA's 6.5 kHz is 4.6 degrees of crank per step at 5000 rpm; neither
model is converged there (section 6), and small changes move it by the
scatter.

Where the time goes: each simulation step solves every port boundary twice
per fluid sub-step (the midpoint rule, 16 times per step) and advances the
pipes and junctions once per CFL step. The port face search takes no `pow`
(section 2.4; searching in `p_f` took about 15 per face and a quarter of the
run time of the V8), the junctions one `log` per pressure and one `exp` per
end. Holding the pipes cut their steps to about 2 per simulation step at
20 kHz and 4 at 10 kHz on 4 cm cells, so the gain is largest at high
simulation frequencies (the EJ25 at 20 kHz: 2.5x); on cells whose CFL step
is shorter than the fluid sub-step (1 cm at 20 kHz, 2 cm at 10 kHz) the pipes
set the pace and holding changes little. The ports' 16 face solves per step
are now a floor under the cost of every cell size. The earlier, loosely
coupled scheme cost 1.85 CPU s for the LS at 10 kHz with 4 cm cells but did
not converge. Bandwidth: the shortest wavelength carried cleanly is about
8-10 cells, so 4 cm resolves up to ~1-1.5 kHz in hot gas (c ~ 600-800 m/s)
and 2 cm up to ~2-3 kHz; the firing harmonics of every shipped engine are
inside the 4 cm range, the hiss and crack above 2 kHz are left to the impulse
response.

Kept: one half-step predictor per cell (not per face), no divisions in the
hot loop, cached orifice constants. Not done: single precision, SIMD across
cells, one exhaust system per thread. On a slow machine use
`legacy_exhaust: true`.

## 6. Physics check: primary length

2.0 L inline-4, 4-1 header, 721 mm x 44 mm primaries, 488 mm x 79.5 mm
collector (`test/data/exhaust_tuner_i4.mr`, 20 kHz, rev limit raised to 9000
for the 7800 rpm point), wide open, dyno hold, mean torque over the last 2 s
of a 4 s run. Primaries scaled by 0.9 / 1.0 / 1.1. Torque in N m:

| model | primaries | 4000 rpm | 5000 rpm | 6000 rpm | 7000 rpm | 7800 rpm |
|---|---|---|---|---|---|---|
| legacy | 649 mm (x0.9) | 144.0 | 173.7 | 193.0 | 184.1 | 166.1 |
| legacy | 721 mm (x1.0) | 145.6 | 178.1 | 190.7 | 183.0 | 163.0 |
| legacy | 793 mm (x1.1) | 148.9 | 178.5 | 190.5 | 181.4 | 159.9 |
| 1-D, 4 cm | 649 mm (x0.9) | 53.7 | 137.2 | 168.9 | 184.4 | 187.3 |
| 1-D, 4 cm | 721 mm (x1.0) | 61.5 | 159.3 | 140.7 | 193.7 | 188.2 |
| 1-D, 4 cm | 793 mm (x1.1) | 97.3 | 125.9 | 153.9 | 198.6 | 188.2 |
| 1-D, 2 cm | 649 mm (x0.9) | 52.4 | 136.6 | 161.5 | 184.2 | 187.4 |
| 1-D, 2 cm | 721 mm (x1.0) | 60.6 | 164.3 | 138.1 | 193.2 | 188.2 |
| 1-D, 2 cm | 793 mm (x1.1) | 96.7 | 123.5 | 152.5 | 198.0 | 188.0 |

The legacy model changes by 1-4% with the length (the longer primary is a
bigger lumped volume); the 1-D model by 17-33% at 5000-6000 rpm and 80-85%
at 4000 rpm, and the three lengths come together at 7800 rpm, where the tuned
pulse falls outside the overlap. 4 cm and 2 cm cells give the same picture:
every length change from 4000 to 7000 rpm goes the same way on both, and
single runs agree within 4.6% (the largest difference is the 649 mm
primaries at 6000 rpm; the random burn alone moves a run by up to 2.4%,
below). The pipes are held between their CFL steps (section 2.2); with them
in lock-step with the fluid sub-steps the same grid agreed within 3.7%, the
largest difference at 721 mm and 5000 rpm.

Convergence, 721 mm primaries at 6000 rpm: 140.7 / 138.1 / 137.5 / 136.9 N m
on 4 / 2 / 1 / 0.5 cm cells (20 kHz); 140.2 / 138.1 / 139.4 N m at 15 / 20 /
40 kHz on 2 cm cells, 140.7 / 141.5 at 20 / 40 kHz on 4 cm. The combustion
model scales each burn's efficiency by a random number
(`CombustionChamber::ignite`, `std::rand`, seeded once per run): four seeds
put the 4 cm result between 140.7 and 143.5 N m (mean 141.6) and the 2 cm
result between 138.1 and 141.4 N m (mean 139.2), so differences of a few N m
between runs are that scatter. At 4000 rpm: 61.5 / 60.6 / 60.1 N m on 4 / 2 /
1 cm. With the pipes in lock-step the same points gave 141.2 / 138.5 / 137.5 /
136.9 and 137.8 / 138.5 / 139.5 N m, seeds 141.0-143.9 (mean 141.9) and
138.5-142.0 N m (mean 140.0): the same within the scatter, the held coarse
grids a little nearer the fine ones (section 4). Below 15 kHz the torque falls
away in both models (at 10 kHz, 3.6 degrees of crank per step: 97.5 against
138.1 N m for the 1-D model on 2 cm cells, 131.5 against 190.7 N m for the
legacy one): engine-sim's own cylinder,
valve and port flows need a step of no more than about 2 degrees at the
highest rpm used (the 2.4 degrees of 15 kHz at 6000 rpm, or of the LS V8 at
4000 rpm and 10 kHz, were still converged; legacy 194.3 / 190.7 / 193.2 N m
at 15 / 20 / 40 kHz). Before
the port boundary of section 2.4 the same 721 mm, 6000 rpm point gave
148.7 / 183.6 / 118.5 / -6.7 N m on 4 / 2 / 1 / 0.5 cm, and the two grids
disagreed on the direction of the length change.

The low-rpm hole (4000 rpm, 721 mm) is the 4-1 interference. Traced over a
cycle of cylinder 1 (2 cm cells): the cylinder that fires 180 degrees later
opens its exhaust at 300 degrees of cylinder 1's cycle; its blowdown pulse
runs 0.72 m to the junction and 0.72 m back up cylinder 1's primary, and
reaches cylinder 1's port at overlap TDC: the port pressure peaks at 1.52 bar
at 367 degrees (both valves are open past 0.25 mm from 311 to 411 degrees;
the legacy model's port stays at 1.02 bar). It pushes exhaust gas back into
the cylinder (0.0021 mol per cycle, 14% of what left) and on through the
open intake valve into the lumped intake runner (0.0019 mol during overlap).
The intake stroke then draws 0.0231 mol instead of 0.0432, and at intake
closing the cylinder holds 53% of the legacy model's fuel; torque 60 against
146 N m. With 793 mm primaries the pulse peaks later, at 1.78 bar at 390
degrees with the exhaust valve almost shut (0.8 mm): 67% of the legacy
model's fuel and 97 N m. A real 4-1 header shows this interference; its depth here is larger
than a dyno would show because the junction has no losses, the wall
temperature is fixed, and the intake runner is one lumped volume, so the
exhaust gas pushed into it stays in the next charge instead of being swept
back by an intake wave.

## 7. Sound

WAV previews rendered through the headless runner (same synthesizer settings
as the app, impulse responses from `es/sound-library`), held rpm, the last
1.5 s analysed with a Hann window; line strengths are amplitude ratios to the
strongest line, band edges kept off the firing harmonics (the I4 fires at
200 Hz at 6000 rpm, the V8 at 267 Hz at 4000 rpm):

| run | RMS | clipped samples | strongest lines (Hz : relative) | energy 150-350 / 350-1000 / 1000-3000 Hz |
|---|---|---|---|---|
| I4 6000 rpm, legacy | 27558 | 54.1% | 200:1.00 600:0.15 400:0.07 800:0.06 1600:0.02 1400:0.01 | 97% / 3% / 0% |
| I4 6000 rpm, 1-D 4 cm | 13692 | 2.6% | 400:1.00 200:0.64 800:0.63 600:0.33 1600:0.28 1200:0.21 | 20% / 72% / 8% |
| V8 (LS) 4000 rpm, legacy | 17527 | 3.0% | 433:1.00 267:0.36 400:0.28 167:0.23 300:0.16 100:0.14 | 16% / 82% / 0% |
| V8 (LS) 4000 rpm, 1-D 4 cm | 15131 | 1.9% | 633:1.00 433:0.69 1000:0.59 400:0.50 567:0.39 700:0.35 | 4% / 76% / 19% |

The legacy I4 is a nearly pure 200 Hz tone that the leveler clips more than
half the time; the 1-D I4 spreads its energy over the firing harmonics (200,
400, 600, 800 Hz and on) and clips little. On the V8 the 1-D model moves energy
from the 433 Hz line into the lines from 400 Hz to 1 kHz and carries 19% of it
above 1 kHz (the blowdown pulse shape). With the pipes held between their CFL
steps (section 2.2) the 1-D renders keep a little more above 1 kHz than in
lock-step (I4 8% against 4%, V8 19% against 17%; the lines up to 1 kHz
move by at most 0.06 of the strongest, the I4's 1600 Hz line rises from 0.18
to 0.28), as the less diffused waves of section 4 would have it. Files: `i4_legacy_6000.wav`,
`i4_1d_6000.wav`, `ls_legacy_4000.wav`, `ls_1d_4000.wav` in the harness output
(not committed). The legacy renders are the same as before this branch (the
I4 render was redone: identical statistics).

## 8. Known limits

* Pulse amplitudes are probably on the strong side: no junction loss
  coefficients (constant-pressure junction), a fixed wall temperature, no
  pulse attenuation other than friction; and the intake side is still a
  lumped runner, so exhaust gas pushed into it stays for the whole intake
  stroke instead of being swept out by an intake wave. Together these make
  the 4-1 low-rpm hole deeper than a dyno would show.
* The port -> primary orifice (`primary_flow_rate`) is inherited from the
  legacy structure. Physically the port and primary are one duct; the
  orifice keeps the old parameter meaningful, and where its `Cd A` is below
  the pipe area the boundary adds a sudden-expansion loss at the pipe inlet
  (the jet's excess kinetic energy is mixed out at the face pressure). The
  calculator's exports carry the donor engine's rating: `k_carb(200)` gives
  `Cd A` = 10.8 cm^2, 70% of the test engine's 44 mm primary (15.5 cm^2).
* engine-sim itself (cylinder, valve and port flows, both exhaust models)
  needs a simulation step of no more than about 2 degrees of crank at the
  highest rpm used; at 3.6 degrees per step the test inline-4 lost a third of
  its torque in both models (section 6). Several shipped engines run at
  5-6.5 kHz, i.e. 4-6 degrees per step at 5000 rpm; their torque is not
  converged in either model, and the 1-D model departs further from the
  legacy one there.
* The port boundaries are solved at the engine's fluid sub-step (8 x
  `simulation_frequency`) whatever the cell size, so the ports set a floor
  on the cost (section 5).
* Cones are limited to a linear taper inside one pipe; steps need two
  pipes; no X/H crossovers between exhaust systems, no mufflers or chambers
  (the impulse response still does that), no catalyst.
* `sound_attenuation` has no effect on 1-D systems.
* Primaries are straight pipes of one diameter; the head's port volume is
  lumped.
* Cell size trades CPU for bandwidth and accuracy (sections 5 and 6). The
  200-cell cap coarsens very long pipes silently.
* The radiated sound is a monopole at 1 m; no directivity, no ground.
