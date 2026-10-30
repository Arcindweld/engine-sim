# engine-sim (Arcindweld fork)

A fork of [ange-yaghi/engine-sim](https://github.com/ange-yaghi/engine-sim) (MIT), the real-time engine
simulator and sound synthesizer, imported at upstream commit `85f7c3b` with its five submodules pinned
(csv-io `2112c55`, delta-studio `b7d0a04`, direct-to-video `19f939c`, piranha `432f0b1`,
simple-2d-constraint-solver `e009f4f`). The upstream README is kept as `UPSTREAM-README.md`.

The fork exists to pair the simulator with the Arcindweld exhaust and intake design calculator, which exports
`.mr` engine files for it. Work planned here, in order:

1. **1-D exhaust model** (branch `exhaust-1d`): primaries, collector and tailpipe as pipes with real wave
   propagation, sound taken at the tailpipe end, so header length and collector design change the sound and
   the torque curve the way they do on a car. Legacy mode stays available and unchanged.
2. **Flow coefficients from geometry**: valve, port and pipe discharge coefficients from the engine's own
   dimensions instead of per-file guesses (the exporter already ships donor orifice ratings because ideal
   orifices gave a mild 350 small-block 913 hp).
3. **Performance**: profiling and the hot loops.

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
