use std::collections::VecDeque;

use crate::spin::{Spin, SmallVec};
use rand::random;

pub struct Lattice<S: Spin> {
    data: Vec<S>,
    neighbors: Vec<SmallVec<usize>>,
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
        }
    }

    pub fn zero(&mut self) {
        for item in self.data.iter_mut() {
            *item = S::start();
        }
    }

    pub fn run(&mut self, beta: f64, n_trials: usize) -> f64 {
        let rolling_avg_count = 500;
        let runs_per_trial = 16;
        
        self.evolve(beta, runs_per_trial * (n_trials - rolling_avg_count));
        let mut rolling_avg = 0.0;
        for _ in 0..rolling_avg_count {
            self.evolve(beta, runs_per_trial);
            rolling_avg += self.order_param();
        }
        rolling_avg / rolling_avg_count as f64
    }

    fn evolve(&mut self, beta: f64, count: usize) {
        for _ in 0..count {
            // 1. Choose random site
            let start_index = random::<usize>() % self.data.len();
            let vec = S::random_vec();

            // 2. Mark, and flip.
            self.data[start_index].flip(&vec);
            let mut new_sites = VecDeque::new();
            new_sites.push_back(start_index);
            let mut marked_sites = vec![false; self.data.len()];

            while let Some(s) = new_sites.pop_front() {
                // 3. Iterate through neighbors and flip
                for neighbor in self.neighbors[s].iter_mut() {
                    if marked_sites[*neighbor] { continue; }
                    let prob = 1. - (0.0f64.min(2. * beta * self.data[s].dot(&vec) * self.data[*neighbor].dot(&vec))).exp();
                    if random::<f64>() < prob {
                        self.data[*neighbor].flip(&vec);
                        new_sites.push_back(*neighbor);
                        marked_sites[*neighbor] = true;
                    }
                }
            }
        }
    }

    pub fn order_param(&self) -> f64 {
        let mut tot = S::zero();
        for item in &self.data {
            tot += item;
        }
        tot.norm() / self.data.len() as f64
    }
}