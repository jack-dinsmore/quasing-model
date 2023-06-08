use std::fmt::Display;

use crate::spin::{Spin, SmallVec};
use rand::random;

const CLUSTER_QUEUE_MAX_SIZE: usize = 256_000;
const BURN_IN: usize = 128;
const RUNS_PER_TRIAL: usize = 16;

/// Result from a complete run of a lattice
pub struct Report {
    pub magnetization: f32,
    pub susceptibility: f32,
}

impl Display for Report {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Magnetization {}", self.magnetization)?;
        writeln!(f, "Susceptibility {}", self.susceptibility)
    }
}

pub struct Lattice<S: Spin> {
    data: Vec<S>,
    neighbors: Vec<SmallVec<(usize, f32)>>,
    cluster_queue: [usize; CLUSTER_QUEUE_MAX_SIZE],
}

impl<S: Spin> Lattice<S> {
    pub fn new(num_sites: usize, neighbor_func: &impl Fn(usize) -> SmallVec<(usize, f32)>) -> Self {
        let mut data = Vec::with_capacity(num_sites);
        let mut neighbors = Vec::with_capacity(num_sites);
        for site_index in 0..num_sites {
            data.push(S::start());
            neighbors.push(neighbor_func(site_index));
        }

        Self {
            data,
            neighbors,
            cluster_queue: [0; CLUSTER_QUEUE_MAX_SIZE],
        }
    }

    pub fn zero(&mut self) {
        for item in self.data.iter_mut() {
            *item = S::start();
        }
    }

    pub fn run(&mut self, beta: f32, n_trials: usize, checkpoint: (usize, f32)) -> Report {
        let mut num_avgs = n_trials - BURN_IN;

        self.evolve(beta, RUNS_PER_TRIAL * BURN_IN);
        let mut avg_mag = 0.0;
        let mut avg_mag2 = 0.0;
        for trial_index in BURN_IN..n_trials {
            self.evolve(beta, RUNS_PER_TRIAL);
            let mag = self.magnetization();
            avg_mag += mag;
            avg_mag2 += mag * mag;

            if trial_index == checkpoint.0 {
                if avg_mag / (trial_index - BURN_IN + 1) as f32 > checkpoint.1 {
                    num_avgs = trial_index - BURN_IN + 1;
                    break;
                }
            }
        }
        avg_mag /= num_avgs as f32;
        avg_mag2 /= num_avgs as f32;
        let sus = (avg_mag2 - avg_mag*avg_mag) * self.data.len() as f32;
        Report {
            magnetization: avg_mag,
            susceptibility: sus,
        }
    }

    fn evolve(&mut self, beta: f32, count: usize) -> f32 {
        let mut cluster_size = 0;
        for _ in 0..count {
            // 1. Choose random site
            let start_index = random::<usize>() % self.data.len();
            let vec = S::random_vec();

            // 2. Mark, and flip.
            let mut stack_pointer: i32 = 0;
            self.cluster_queue[0] = start_index;
            self.data[start_index].flip(&vec, true);
            let mut marked_sites = vec![false; self.data.len()];

            while stack_pointer >= 0 {
                // 3. Iterate through neighbors and flip
                let my_index = self.cluster_queue[stack_pointer as usize];
                let mut new_stack_pointer = stack_pointer - 1;
                for (neighbor, strength) in self.neighbors[my_index].iter_mut() {
                    // if marked_sites[*neighbor] { continue; }
                    let prob = 1. - (0.0f32.min(2. * beta * strength
                        * self.data[my_index].dot(&vec)
                        * self.data[*neighbor].dot(&vec))
                    ).exp();
                    
                    let success = !marked_sites[*neighbor] && (random::<f32>() < prob);
                    self.cluster_queue[new_stack_pointer as usize + 1] = *neighbor;
                    self.data[*neighbor].flip(&vec, success);
                    new_stack_pointer += success as i32;
                    marked_sites[*neighbor] = marked_sites[*neighbor] || success;
                    cluster_size += success as usize;
                }
                stack_pointer = new_stack_pointer;
            }
        }
        cluster_size as f32 / count as f32 / self.data.len() as f32
    }

    pub fn magnetization(&self) -> f32 {
        let mut tot = S::zero();
        for item in &self.data {
            tot += item;
        }
        tot.norm() / self.data.len() as f32
    }
}