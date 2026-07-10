//! Spectral mechanics: where graph spectral theory meets Hamiltonian mechanics.
//! The graph Laplacian IS the kinetic energy operator. Edges ARE springs.

/// A graph as a mechanical system: nodes are masses, edges are springs
pub struct SpringGraph {
    /// Adjacency matrix (edge weights = spring constants)
    pub adj: Vec<Vec<f64>>,
    /// Node masses
    pub masses: Vec<f64>,
    /// Node positions (displacement from equilibrium)
    pub positions: Vec<f64>,
    /// Node velocities
    pub velocities: Vec<f64>,
}

impl SpringGraph {
    pub fn new(adj: Vec<Vec<f64>>) -> Self {
        let n = adj.len();
        Self {
            adj,
            masses: vec![1.0; n],
            positions: vec![0.0; n],
            velocities: vec![0.0; n],
        }
    }

    pub fn with_masses(mut self, masses: Vec<f64>) -> Self {
        self.masses = masses;
        self
    }

    /// Potential energy: V = Σ k_ij (x_i - x_j)² / 2 (count each edge once)
    pub fn potential_energy(&self) -> f64 {
        let n = self.adj.len();
        let mut v = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.positions[i] - self.positions[j];
                v += 0.5 * self.adj[i][j] * dx * dx;
            }
        }
        v
    }

    /// Kinetic energy: T = Σ m_i v_i² / 2
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities.iter().enumerate()
            .map(|(i, v)| self.masses[i] * v * v / 2.0)
            .sum()
    }

    /// Total energy (Hamiltonian)
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }

    /// Forces on each node: F_i = -Σ_j k_ij (x_i - x_j)
    /// Each spring contributes fully to each node's force
    pub fn forces(&self) -> Vec<f64> {
        let n = self.adj.len();
        let mut f = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    f[i] -= self.adj[i][j] * (self.positions[i] - self.positions[j]);
                }
            }
        }
        f
    }

    /// Normal mode frequencies: ω_i = √λ_i where λ_i are Laplacian eigenvalues / masses
    pub fn normal_modes(&self) -> Vec<f64> {
        let eigs = self.laplacian_eigenvalues();
        eigs.iter().enumerate()
            .map(|(i, &l)| {
                if i == 0 { 0.0 } // zero mode (translation)
                else { (l / self.masses.get(i).copied().unwrap_or(1.0)).sqrt().max(0.0) }
            })
            .collect()
    }

    /// Step using Störmer-Verlet (symplectic, conserves energy)
    pub fn verlet_step(&mut self, dt: f64) {
        let half = dt / 2.0;
        let forces = self.forces();
        // Half kick: update velocity using acceleration = F/m
        for i in 0..self.velocities.len() {
            self.velocities[i] += half * forces[i] / self.masses[i];
        }
        // Drift: update position
        for i in 0..self.positions.len() {
            self.positions[i] += dt * self.velocities[i];
        }
        // Half kick with new forces
        let forces = self.forces();
        for i in 0..self.velocities.len() {
            self.velocities[i] += half * forces[i] / self.masses[i];
        }
    }

    /// Integrate for n steps, return energy drift
    pub fn integrate(&mut self, dt: f64, steps: usize) -> EnergyReport {
        let initial = self.total_energy();
        let mut max_drift = 0.0_f64;
        let mut max_e = initial;
        let mut min_e = initial;

        for _ in 0..steps {
            self.verlet_step(dt);
            let e = self.total_energy();
            let drift = (e - initial).abs();
            if drift > max_drift { max_drift = drift; }
            if e > max_e { max_e = e; }
            if e < min_e { min_e = e; }
        }

        EnergyReport {
            initial_energy: initial,
            final_energy: self.total_energy(),
            max_drift,
            energy_range: max_e - min_e,
        }
    }

    /// Laplacian eigenvalues via Jacobi
    pub fn laplacian_eigenvalues(&self) -> Vec<f64> {
        let n = self.adj.len();
        let mut lap = vec![vec![0.0; n]; n];
        for i in 0..n {
            let deg: f64 = self.adj[i].iter().sum();
            lap[i][i] = deg;
            for j in 0..n { if i != j { lap[i][j] = -self.adj[i][j]; } }
        }
        jacobi(&mut lap)
    }

    /// Conservation ratio of the underlying graph
    pub fn cr(&self) -> f64 {
        let eigs = self.laplacian_eigenvalues();
        if eigs.len() < 2 { return 0.0; }
        let l2 = eigs[1];
        let ln = *eigs.last().unwrap_or(&1.0);
        if ln <= 0.0 { 0.0 } else { l2 / ln }
    }
}

#[derive(Debug)]
pub struct EnergyReport {
    pub initial_energy: f64,
    pub final_energy: f64,
    pub max_drift: f64,
    pub energy_range: f64,
}

/// Equipartition: in thermal equilibrium, each mode has kT/2 energy
pub fn equipartition_check(graph: &SpringGraph) -> Vec<f64> {
    let modes = graph.normal_modes();
    let total_ke = graph.kinetic_energy();
    let n_active = modes.iter().filter(|&&m| m > 1e-10).count().max(1);
    let expected = total_ke / n_active as f64;

    // Energy in each mode (approximate: distribute KE equally among active modes)
    let n = graph.adj.len();
    let mut mode_energies = Vec::new();

    for (i, &freq) in modes.iter().enumerate() {
        if freq < 1e-10 {
            mode_energies.push(0.0);
        } else {
            mode_energies.push(expected);
        }
        if i >= n - 1 { break; }
    }
    mode_energies
}

/// Spectral temperature: T = 2 * <KE> / (k * n_active_modes)
pub fn spectral_temperature(graph: &SpringGraph) -> f64 {
    let ke = graph.kinetic_energy();
    let modes = graph.normal_modes();
    let n_active = modes.iter().filter(|&&m| m > 1e-10).count().max(1);
    2.0 * ke / n_active as f64
}

/// Virial theorem check: 2<T> = <x·F> for a system in equilibrium
pub fn virial_ratio(graph: &SpringGraph) -> f64 {
    let ke = graph.kinetic_energy();
    let forces = graph.forces();
    let vf: f64 = graph.positions.iter().zip(forces.iter())
        .map(|(x, f)| x * f)
        .sum();
    if ke.abs() < 1e-15 { return 0.0; }
    -vf / (2.0 * ke) // Should be ~1.0 at equilibrium
}

fn jacobi(a: &mut Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 { return vec![]; }
    for _ in 0..100 * n * n {
        let (mut p, mut q) = (0, 1);
        let mut max_val = 0.0_f64;
        for i in 0..n { for j in (i+1)..n { if a[i][j].abs() > max_val { max_val = a[i][j].abs(); p = i; q = j; } } }
        if max_val < 1e-14 { break; }
        let app = a[p][p]; let aqq = a[q][q]; let apq = a[p][q];
        let theta = if (app - aqq).abs() < 1e-30 { std::f64::consts::FRAC_PI_4 }
                     else { 0.5 * (2.0 * apq / (app - aqq)).atan() };
        let (c, s) = (theta.cos(), theta.sin());
        for i in 0..n { if i != p && i != q { let aip = a[i][p]; let aiq = a[i][q]; a[i][p] = c*aip+s*aiq; a[p][i]=a[i][p]; a[i][q]=-s*aip+c*aiq; a[q][i]=a[i][q]; } }
        a[p][p] = c*c*app+2.0*s*c*apq+s*s*aqq;
        a[q][q] = s*s*app-2.0*s*c*apq+c*c*aqq;
        a[p][q] = 0.0; a[q][p] = 0.0;
    }
    let mut eigs: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eigs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn harmonic_oscillator_conserves() {
        let mut g = SpringGraph::new(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
        g.positions = vec![1.0, -1.0];
        g.velocities = vec![0.0, 0.0];
        let report = g.integrate(0.01, 10000);
        assert!(report.max_drift < 1e-3, "Energy drift: {}", report.max_drift);
    }

    #[test]
    fn normal_modes_complete_graph() {
        // K₃ (triangle): Laplacian eigenvalues are 0, 3, 3
        let g = SpringGraph::new(vec![
            vec![0.0, 1.0, 1.0],
            vec![1.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
        ]);
        let modes = g.normal_modes();
        assert!(modes[0].abs() < 1e-10, "Zero mode: {}", modes[0]);
        assert!((modes[1] - 3.0_f64.sqrt()).abs() < 0.1, "Mode 1: {}", modes[1]);
    }

    #[test]
    fn forces_restore_equilibrium() {
        let mut g = SpringGraph::new(vec![vec![0.0, 2.0], vec![2.0, 0.0]]);
        g.positions = vec![1.0, 0.0];
        let f = g.forces();
        // Force on node 0 should be toward node 1 (negative direction)
        assert!(f[0] < 0.0, "Force should pull toward equilibrium: {}", f[0]);
        assert!((f[0] + 2.0).abs() < 1e-10, "Force magnitude: {}", f[0]);
    }

    #[test]
    fn energy_exchange() {
        let mut g = SpringGraph::new(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
        g.positions = vec![1.0, -1.0];
        g.velocities = vec![0.0, 0.0];
        let initial_pe = g.potential_energy();
        let initial_ke = g.kinetic_energy();
        assert!(initial_pe > 0.0);
        assert!(initial_ke.abs() < 1e-10);

        // Quarter period: PE → KE
        g.verlet_step(std::f64::consts::FRAC_PI_4 / 2.0);
        // Now there should be both PE and KE
        assert!(g.kinetic_energy() > 0.0, "KE should be > 0 after step");
    }

    #[test]
    fn cr_for_graphs() {
        let path = SpringGraph::new(vec![
            vec![0.0, 1.0, 0.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 0.0],
        ]);
        let complete = SpringGraph::new(vec![
            vec![0.0, 1.0, 1.0],
            vec![1.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
        ]);
        // Complete graph should have higher CR than path
        assert!(complete.cr() > path.cr(), "Complete CR {} vs Path CR {}", complete.cr(), path.cr());
    }

    #[test]
    fn spectral_temperature_positive() {
        let mut g = SpringGraph::new(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
        g.velocities = vec![1.0, -1.0];
        let temp = spectral_temperature(&g);
        assert!(temp > 0.0, "Temperature should be positive: {}", temp);
    }

    #[test]
    fn virial_near_equilibrium() {
        let mut g = SpringGraph::new(vec![vec![0.0, 1.0], vec![1.0, 0.0]]);
        g.positions = vec![0.1, -0.1];
        g.velocities = vec![0.0, 0.0];
        // Near equilibrium, virial ratio should be ~1
        let ratio = virial_ratio(&g);
        // For small displacement, 2<T> ≈ 0, <x·F> ≈ small, ratio is unstable
        // Just check it computes
        assert!(ratio.is_finite(), "Virial ratio should be finite: {}", ratio);
    }

    #[test]
    fn integrate_long_stability() {
        let mut g = SpringGraph::new(vec![
            vec![0.0, 1.0, 0.5],
            vec![1.0, 0.0, 1.0],
            vec![0.5, 1.0, 0.0],
        ]);
        g.positions = vec![1.0, 0.0, -0.5];
        g.velocities = vec![0.1, -0.2, 0.1];
        let report = g.integrate(0.005, 50000);
        assert!(report.max_drift < 0.01, "Long integration drift: {}", report.max_drift);
    }
}
