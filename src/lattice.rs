use std::fmt::Display;

use crate::spin::{Spin, SmallVec};
use rand::random;

const CLUSTER_QUEUE_MAX_SIZE: usize = 10000;
const BURN_IN: usize = 256;
const RUNS_PER_TRIAL: usize = 16;

pub struct Report {
    pub order_parameter: f64,
    pub order_parameter_err: f64,
    pub susceptibility: f64,
    pub susceptibility_err: f64,
}

impl Display for Report {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Order parameter {} +/- {}", self.order_parameter, self.order_parameter_err)?;
        writeln!(f, "Susceptibility {} +/- {}", self.susceptibility, self.susceptibility_err)
    }
}

pub struct Lattice<S: Spin> {
    data: Vec<S>,
    neighbors: Vec<SmallVec<usize>>,
    cluster_queue: [usize; CLUSTER_QUEUE_MAX_SIZE],
}

impl<S: Spin> Lattice<S> {
    pub fn new(num_sites: usize, neighbor_func: impl Fn(usize) -> SmallVec<usize>) -> Self {
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

    pub fn run(&mut self, beta: f64, n_trials: usize) -> Report {
        let num_avgs = n_trials - BURN_IN;

        self.evolve(beta, RUNS_PER_TRIAL * BURN_IN);
        let mut avg_param = 0.0;
        let mut avg_sus = 0.0;
        let mut avg_param2 = 0.0;
        let mut avg_sus2 = 0.0;
        for _ in BURN_IN..n_trials {
            let sus = self.evolve(beta, RUNS_PER_TRIAL);
            let param = self.order_param();
            avg_sus += sus;
            avg_param += param;
            avg_sus2 += sus*sus;
            avg_param2 += param*param;
        }
        avg_param /= num_avgs as f64;
        avg_sus /= num_avgs as f64;
        avg_param2 /= num_avgs as f64;
        avg_sus2 /= num_avgs as f64;
        Report {
            order_parameter: avg_param,
            susceptibility: avg_sus,
            order_parameter_err: (avg_param2 - avg_param*avg_param).sqrt(),
            susceptibility_err: (avg_sus2 - avg_sus*avg_sus).sqrt(),
        }
    }

    fn evolve(&mut self, beta: f64, count: usize) -> f64 {
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
                for neighbor in self.neighbors[my_index].iter_mut() {
                    // if marked_sites[*neighbor] { continue; }
                    let prob = 1. - (0.0f64.min(2. * beta
                        * self.data[my_index].dot(&vec)
                        * self.data[*neighbor].dot(&vec))
                    ).exp();

                    let success = !marked_sites[*neighbor] & (random::<f64>() < prob);
                    self.cluster_queue[new_stack_pointer as usize + 1] = *neighbor;
                    self.data[*neighbor].flip(&vec, success);
                    new_stack_pointer += success as i32;
                    marked_sites[*neighbor] = success;
                    cluster_size += success as usize;
                }
                stack_pointer = new_stack_pointer;
            }
        }
        cluster_size as f64 / count as f64 / self.data.len() as f64
    }

    pub fn order_param(&self) -> f64 {
        let mut tot = S::zero();
        for item in &self.data {
            tot += item;
        }
        tot.norm() / self.data.len() as f64
    }
}