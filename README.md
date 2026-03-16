# Probabilistic Allocation of Multi-Mapping Reads

This project demonstrates how to use an Expectation-Maximization (EM) algorithm to probabilistically allocate ambiguous RNA-seq reads across repetitive genomic features (for example, L1 insertions).  
The current implementation is notebook-first and focused on a synthetic simulation workflow in `em.ipynb`.

## Why This Project Exists

Reads that map to multiple repetitive loci are hard to assign with a single "best" alignment rule.  
Instead of forcing a hard assignment, this project estimates fractional assignments using:

- alignment likelihoods per read-feature pair, and
- iteratively updated feature activity priors (`theta`).

This gives a principled estimate of expected counts per feature.

## Method (EM Overview)

Let:

- `L` be the read-by-feature likelihood matrix (`N x M`)
- `theta` be the current prior weights over `M` features
- `Z` be the posterior fractional assignment matrix

E-step:

$$
Z_{ij} = \frac{L_{ij}\theta_j}{\sum_{k=1}^{M} L_{ik}\theta_k}
$$

M-step:

$$
\theta_j^{(\text{new})} = \frac{\sum_{i=1}^{N} Z_{ij}}{N}
$$

The notebook iterates these updates and visualizes convergence dynamics.

## Repository Layout

- `em.ipynb`: main notebook with formulation, simulation, EM loop, and plots
- `.devcontainer/devcontainer.json`: VS Code/Cursor devcontainer config
- `.devcontainer/Dockerfile`: container image and dependency installation
- `.devcontainer/requirements.txt`: Python dependencies
- `LICENSE`: GNU GPLv3 license

## Quickstart

### Option A: Dev Container (recommended)

1. Open the repository in Cursor or VS Code.
2. Reopen in container (uses `.devcontainer/devcontainer.json`).
3. Open `em.ipynb` and run all cells.

The container installs dependencies from `.devcontainer/requirements.txt`.

### Option B: Local Python Environment

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r .devcontainer/requirements.txt
```

## Run the Notebook

From the repository root:

```bash
jupyter notebook em.ipynb
```

or:

```bash
jupyter lab em.ipynb
```

Then run all cells in order.

## What You Should See

The notebook currently:

- builds a synthetic likelihood matrix for 1000 reads across 5 insertions,
- initializes `theta` uniformly,
- runs EM updates for 600 iterations,
- prints periodic `theta` updates and final expected counts,
- plots likelihood structure and `theta` evolution over iterations.

These outputs help show how anchor reads and mismatch penalties create informative dynamics.

## Assumptions and Caveats

- This is a prototype notebook, not a production pipeline.
- The current demonstration uses synthetic data only.
- If all features are equally likely for all reads, EM has no informative dynamics (the uniform ambiguity caveat).
- Results depend on the quality and calibration of the likelihood matrix `L`.

## Roadmap

Potential next improvements:

- refactor notebook logic into reusable Python modules,
- add automated tests for E-step/M-step correctness,
- add a CLI for batch runs,
- support external/real alignment likelihood inputs.

## License

This project is licensed under GNU GPLv3. See `LICENSE` for full terms.
