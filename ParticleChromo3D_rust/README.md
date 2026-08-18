# ParticleChromo3D (Rust)

Rust port of `ParticleChromo3D/Ps.py`. Runs the same PSO algorithm to reconstruct 3D chromosome structure from Hi-C contact matrices, with parallelism over the alpha sweep via Rayon.

## Install from crates.io

```bash
cargo install particle_chromo3d
```

Then run directly:

```bash
particle_chromo3d <input_matrix> [OPTIONS]
```

### Example

```bash
particle_chromo3d chr22_matrix.txt -s 15 -i 30000 -o ./out/chr22
```

## Use as a library

```bash
cargo add particle_chromo3d --no-default-features
```

`--no-default-features` drops the CLI-only dependencies (`clap`, `env_logger`, `uuid`),
leaving just `log`, `rand`, and `rayon`.

```rust
use particle_chromo3d::swarm::LossFunc;
use particle_chromo3d::{PsoParams, alpha_grid, par_choice};

let params = PsoParams {
    rand_range: 1.0,
    swarm_size: 15,
    threshold: 1e-6,
    itt_count: 30000,
    loss_func: LossFunc::Rmse,
};
let alphas = alpha_grid(0.1, 2.0, 0.1);
let best = par_choice("chr22_matrix.txt", "./out/chr22", &alphas, params);

println!("spearman {} at alpha {}", best.spearman, alphas[best.alpha_idx]);
```

Lower-level pieces are available too: `helper` for matrix and PDB I/O, and `swarm`
for the `Swarm` type and the `pdist` / `pearsonr` / `spearmanr` / `loss_function`
primitives but frankly I dont see many use cases for these.

## Requirements (build from source)

- Rust 1.88+ (install via [rustup](https://rustup.rs/) or `conda install conda-forge::rust`)

## Build

From the repo root or from this directory:

```bash
cargo build --release
```

The binary is written to `target/release/particle_chromo3d`.

For a debug build (slower, more logging detail available):

```bash
cargo build
# binary at target/debug/particle_chromo3d
```

## Run

```bash
cargo run --release -- <input_matrix> [OPTIONS]
```

Or call the binary directly after building:

```bash
./target/release/particle_chromo3d <input_matrix> [OPTIONS]
```

### Example

```bash
cargo run --release -- ../exampleIfs/chr22_matrix.txt -s 15 -i 30000 -o ./out/chr22
```

### Options

| Flag | Long | Default | Description |
|------|------|---------|-------------|
| (positional) | | | Input contact matrix file (space-delimited) |
| `-s` | `--swarmSize` | `5` | Number of particles per swarm |
| `-i` | `--ittCount` | `30000` | Maximum iterations |
| `-t` | `--threshold` | `0.000001` | Early-stop threshold (min cost improvement per 1000 iters) |
| `-r` | `--randRange` | `1.0` | Initial xyz coordinate range `[-r, r]` |
| `-o` | `--outfile` | `./out/chr.pdb` | Output PDB filename prefix |
| `-l` | `--lossFunction` | `2` | Loss function: `0`=SSE, `1`=MSE, `2`=RMSE, `3`=Huber |
| | `--logLevel` | `info` | Log verbosity: `error`, `warn`, `info`, `debug`, `trace` |

## Output

- `<outfile><uuid>.pdb` — 3D bead coordinates in PDB format, scaled to [-10, 10]
- `<outfile><uuid>.log` — best alpha, cost, Spearman, and Pearson correlations

## Algorithm

Mirrors the Python implementation:

1. Strips the input matrix of extra whitespace
2. Converts the upper-triangle of the contact matrix to a contact list `(i, j, freq)`
3. Sweeps alpha from 0.1 to 1.9 in 0.1 steps; for each alpha, target distance = `1 / freq^alpha`
4. Runs PSO for each alpha in parallel (one thread per alpha value via Rayon)
5. Selects the alpha with the best Spearman correlation between reconstructed and target distances
6. Writes the best configuration to a PDB file
