# LuGRE Analysis

MATLAB scripts for GNSS positioning with the public data of [LuGRE](https://etd.gsfc.nasa.gov/our-work/lunar-gnss-receiver-experiment-lugre/) (Lunar GNSS Receiver Experiment). [Mission data on Zenodo](https://zenodo.org/records/16411687).

The project uses the receiver's raw GPS/Galileo pseudoranges recorded in cislunar space to study **Direct Position Estimation (DPE)** and to compare it with least-squares baselines.

## Methods

- **Per-epoch DPE** (`mainEpoch.m`): a grid search over receiver position and clock bias that minimizes a weighted L1 pseudorange cost. It can optionally be centered on and regularized by an Earth–Moon orbit prediction. A per-epoch least-squares (LS) solution serves as the baseline.
- **Multi-epoch batch estimation** (`mainBatch.m`): estimates one orbit state (position and velocity at a reference epoch) for a whole observation window. The state is propagated through Earth–Moon dynamics, and epochs with fewer than four satellites are also used.
  - `dpeIrls`: L1 cost solved by iteratively reweighted least squares, with a spline clock model.
  - `dpeGrid`: coarse-to-fine grid over the 6-D orbit state, with the L1-optimal (weighted median) clock at each epoch.
  - `ls`: the same model with an L2 cost, used as the baseline.
- **Evaluation**: the solutions are compared with the onboard navigation solution and a reference trajectory.

## Data

- LuGRE Level-0 telemetry (`L0/TLM`), read from a local copy of the mission data
- IGS multi-GNSS broadcast ephemerides (`data/rinex/`)
- Lunar ephemeris from the JPL Horizons API (fetched online)
- ICRF/ITRF transformations based on IAU 2000/2006 with IERS Earth orientation data

## Layout

```text
mainEpoch.m, mainBatch.m   entry scripts
lugreConfig.m              shared configuration
src/                       io, gnss, frames, orbit, epoch, batch, truth, plot
tests/                     offline and integration tests
experiments/               ablation and feasibility studies for the batch methods
data/                      broadcast ephemerides and reference trajectory
```

## Usage

Requirements: MATLAB with the Aerospace Toolbox, plus network access to JPL Horizons.

1. Set `cfg.data.root` in `lugreConfig.m` to the local LuGRE data folder, and `cfg.data.window` to an operation window (for example `OP5_0`).
2. Run `mainEpoch` or `mainBatch` from the project root.
3. Run the tests with `addpath('tests'); testEpoch; testBatch`.

## License

See [LICENSE](LICENSE).
