# engine-sim (Arcindweld fork)

A fork of [ange-yaghi/engine-sim](https://github.com/ange-yaghi/engine-sim) (MIT, commit 85f7c3b), the open-source
version of Engine Simulator. The goal of this fork is a more faithful exhaust: the stock model treats each exhaust
system as one lumped volume and uses the pipe length only as an audio delay, so header and pipe changes do not change
the sound. This fork is being developed alongside the Exhaust Tuner calculator and Fusion add-in.

Planned, in order:

1. A one-dimensional gas-dynamics exhaust (primaries, collector, tailpipe as pipes with wave propagation), with the
   sound taken from the pressure at the tailpipe end.
2. Flow coefficients from geometry with a discharge coefficient, and head flow tables from the valve model, instead of
   values copied from example engines.
3. Profiling and performance work on the solver and synthesizer.

## Builds

Every push builds a Windows x64 package on GitHub Actions. Open the **Actions** tab, pick the latest run and download
the `engine-sim-windows-x64` artifact. Unzip it and run `bin/engine-sim-app.exe`. Engine files from the Exhaust Tuner
calculator go in `assets/` and load with `Load Engine`.

Upstream's build notes are in [UPSTREAM-README.md](UPSTREAM-README.md). License: MIT, see [LICENSE](LICENSE).
