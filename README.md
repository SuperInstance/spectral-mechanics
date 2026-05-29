# spectral-mechanics

**Graphs as mechanical systems — nodes are masses, edges are springs, eigenvalues are normal modes.**

Pure Rust, zero dependencies. A library that treats graph spectral theory as Hamiltonian mechanics: the Laplacian IS the potential energy operator, and symplectic integration preserves the total energy.

## What This Gives You

- **SpringGraph** — a graph where nodes have mass and edges have spring constants
- **Energy bookkeeping** — potential, kinetic, and total energy with conservation tracking
- **Störmer-Verlet integration** — symplectic integrator that conserves energy over millions of steps
- **Normal mode analysis** — eigenvalue frequencies of the coupled spring system
- **Thermal analysis** — spectral temperature, equipartition, virial theorem
- **Conservation ratio** — CR = λ₂/λ_max of the underlying graph

## The Core Idea

The graph Laplacian `L = D - A` is the potential energy operator for a spring network. If node positions are `x` and spring constants are the edge weights, then:

```
V = ½ x^T L x
```

This is not an analogy. It's a literal spring-mass system. The eigenvalues of L are the squares of the natural frequencies: `ωᵢ = √λᵢ`. The zero eigenvalue corresponds to translation (the whole graph moves as a rigid body).

## Quick Start

```rust
use spectral_mechanics::{SpringGraph, spectral_temperature, virial_ratio};

// Two masses connected by a spring
let mut g = SpringGraph::new(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
g.positions = vec![1.0, -1.0]; // displaced from equilibrium
g.velocities = vec![0.0, 0.0];

// Energy exchange: PE → KE → PE
let initial_pe = g.potential_energy();
let report = g.integrate(0.01, 10000);

// Energy drift should be tiny (symplectic integrator)
assert!(report.max_drift < 1e-3);

// Normal mode frequencies
let modes = g.normal_modes(); // [0.0, √2] — translation + oscillation

// Conservation ratio of the coupling graph
let cr = g.cr();

// Thermal properties
let temp = spectral_temperature(&g);
let virial = virial_ratio(&g); // ~1.0 at equilibrium
```

## API Reference

### SpringGraph

| Method | Description |
|--------|-------------|
| `SpringGraph::new(adj)` | Create with unit masses |
| `.with_masses(masses)` | Set node masses |
| `.potential_energy()` | V = ½ Σ kᵢⱼ(xᵢ-xⱼ)² |
| `.kinetic_energy()` | T = ½ Σ mᵢvᵢ² |
| `.total_energy()` | H = T + V |
| `.forces()` | Fᵢ = -Σⱼ kᵢⱼ(xᵢ-xⱼ) |
| `.normal_modes()` | ωᵢ = √(λᵢ/mᵢ) |
| `.verlet_step(dt)` | One symplectic step |
| `.integrate(dt, steps)` | Full integration with EnergyReport |
| `.cr()` | Conservation ratio of graph |

### Free Functions

| Function | Description |
|----------|-------------|
| `spectral_temperature(&g)` | T = 2⟨KE⟩ / n_active_modes |
| `virial_ratio(&g)` | 2⟨T⟩/⟨x·F⟩, should be ~1.0 |
| `equipartition_check(&g)` | Energy per mode |

## How It Fits

Part of the SuperInstance spectral ecosystem:

- **[spectral-graph-core](https://github.com/SuperInstance/spectral-graph-core)** — Eigenvalues, CR, Fiedler vectors
- **spectral-mechanics** — Graphs as physics (this repo)
- **[symplectic-spin](https://github.com/SuperInstance/symplectic-spin)** — General symplectic integrators in Rust
- **[symplectic-physics](https://github.com/SuperInstance/symplectic-physics)** — Fortran 2008 symplectic integrators

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
spectral-mechanics = { git = "https://github.com/SuperInstance/spectral-mechanics" }
```

## License

MIT

Part of the [SuperInstance](https://github.com/SuperInstance) ecosystem.
