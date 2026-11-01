# engine-sim (Arcindweld fork)

A fork of [ange-yaghi/engine-sim](https://github.com/ange-yaghi/engine-sim) (MIT), the real-time engine
simulator and sound synthesizer, imported at upstream commit `85f7c3b` with its five submodules pinned
(csv-io `2112c55`, delta-studio `b7d0a04`, direct-to-video `19f939c`, piranha `432f0b1`,
simple-2d-constraint-solver `e009f4f`). The upstream README is kept as `UPSTREAM-README.md`.

The fork exists to pair the simulator with the Arcindweld exhaust and intake design calculator, which exports
`.mr` engine files for it. Work planned here, in order:

1. **1-D exhaust model** (done, see below): primaries, collector and tailpipe as pipes with real wave
   propagation, sound taken at the tailpipe end, so header length and collector design change the sound and
   the torque curve the way they do on a car. Legacy mode stays available and unchanged.
2. **Flow coefficients from geometry**: valve, port and pipe discharge coefficients from the engine's own
   dimensions instead of per-file guesses (the exporter already ships donor orifice ratings because ideal
   orifices gave a mild 350 small-block 913 hp).
3. **Performance**: profiling and the hot loops.

## 1-D exhaust

Merged from branch `exhaust-1d`; full description in `docs/exhaust-1d.md`. Each `exhaust_system` can be
solved as a network of one-dimensional pipes (quasi-1-D Euler equations, MUSCL-Hancock + HLLC, wall friction
and heat transfer): one primary per cylinder (`primary_tube_length` + the cylinder's `primary_length`), a
constant-pressure junction, a collector with an optional megaphone cone and tailpipe, and an open end. The
sound is taken from the tailpipe exit (monopole, d(mass flow)/dt) and still goes through `audio_volume` and the
impulse response. Pipe length, diameter and collector design then change the sound and the torque curve.

Parameters on `exhaust_system` / `exhaust_system_parameters`, all optional:

* `exhaust_model`: `"auto"` (default) runs the 1-D model only when one of the geometry inputs below is given,
  so the shipped engines and existing files run the original lumped model unchanged; `"1d"` and `"legacy"`
  force one or the other. `legacy_exhaust: true` is the same as `"legacy"`.
* `primary_diameter` (default: from the head's `exhaust_runner_cross_section_area`), `collector_diameter`
  (from `collector_cross_section_area`), `collector_length` (default `length`), `megaphone_length`,
  `megaphone_exit_diameter`, `tailpipe_length`, `tailpipe_diameter`, `pipe_cell_size` (default 4 cm),
  `wall_temperature` (700 K).

Each cylinder's exhaust-port volume feeds its primary through a characteristic (valve-type) boundary solved
after every fluid sub-step of the engine, and the pipes advance at their own CFL step with the fluxes the ports
gave them, so the torque converges with the cell size and the simulation frequency: the test inline-4 at
6000 rpm makes 141 / 138 / 137 / 137 N m on 4 / 2 / 1 / 0.5 cm cells and 140 / 138 / 139 N m at 15 / 20 / 40 kHz,
against a spread of about 3 N m between seeds of the combustion model's random burn, and 4 cm and 2 cm cells
agree on every length change and within 5% at 4000-7800 rpm (`docs/exhaust-1d.md`, sections 2.2, 2.4 and 6).

Cost: a V8 with two 4-1 headers at `simulation_frequency: 10000` costs about 3x the lumped model with the
default 4 cm cells and 6x with 2 cm cells, an EJ25 at 20 kHz about 3x with 4 cm cells (2 cm resolves up to
~2-3 kHz instead of ~1-1.5 kHz). engine-sim
itself needs a simulation step of no more than about 2 degrees of crank at the highest rpm used (both exhaust
models lose torque beyond that: the inline-4 at 6000 rpm made a third less at 10 kHz, 3.6 degrees per step);
the calculator's fork exports set the frequency to 3 x redline for that reason.

Known limits: 4-1 headers show a low-rpm torque hole that is real but exaggerated (no junction losses, lumped
intake runner, fixed wall temperature; at 4000 rpm the neighbour's blowdown pulse reaches the port at 1.5 bar
during overlap and pushes exhaust through the cylinder into the lumped intake runner); the per-cylinder
`sound_attenuation` is not applied to 1-D systems; cones are linear tapers, there are no crossovers, mufflers
or catalysts (the impulse response still does that).

## Downloads

Every push to `main` builds a Windows x64 package on GitHub Actions. Open the **Actions** tab, pick the latest
green "Windows build" run, and download the `engine-sim-windows-x64` artifact. Unzip it and run
`bin/engine-sim-app.exe`; the `assets`, `es` and `basic` folders beside `bin` are what it loads. Windows
SmartScreen will ask about an unsigned program the first time.

## Building locally

The workflow in `.github/workflows/windows-build.yml` is the reference build: Visual Studio 2022, vcpkg
(`sdl2 sdl2-image boost-filesystem`, triplet x64-windows), win_flex_bison from Chocolatey, `DISCORD_ENABLED=OFF`.
One thing to know: Boost 1.84 removed `boost::filesystem::path::is_complete()`, which the pinned piranha
submodule calls once in `src/path.cpp` (line 76). Change it to `is_absolute()` (the workflow does this
automatically) or build against Boost 1.83 or older.

## License

MIT, as upstream. See `LICENSE`.
