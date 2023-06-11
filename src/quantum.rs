use crate::spin::{Spin, SmallVec};
use crate::lattice::Report;

use std::fmt::Display;
use rand::random;

const CLUSTER_STACK_MAX_SIZE: usize = 256_000;
const BURN_IN: usize = 4;
const RUNS_PER_TRIAL: usize = 8;
const SPIN_BC: bool = true;
const EPSILON : f32 = 1e-5;

pub struct QLattice {
    length_three: f32,
    data: Vec<Vec<f32>>,
    neighbors: Vec<SmallVec<(usize, f32)>>,
    cluster_stack: [(usize, f32); CLUSTER_STACK_MAX_SIZE],
}

fn verify_sorted(l: &Vec<f32>) -> bool {
    let mut last = None;
    for i in l {
        if let Some(las) = last {
            if *i < las {
                return false;
            }
        }
        last = Some(*i)
    }
    true
}

impl QLattice {
    /// Create a new lattice from a function that generates neighbors.
    pub fn new(num_sites: usize, neighbor_func: &impl Fn(usize) -> SmallVec<(usize, f32)>) -> Self {
        let mut data = Vec::with_capacity(num_sites);
        let mut neighbors = Vec::with_capacity(num_sites);
        let length_three = 2. * (num_sites as f32).sqrt();
        for site_index in 0..num_sites {
            data.push(Vec::new());
            neighbors.push(neighbor_func(site_index));
        }

        Self {
            length_three,
            data,
            neighbors,
            cluster_stack: [(0, 0.0); CLUSTER_STACK_MAX_SIZE],
        }
    }

    /// Zero out the data set
    pub fn zero(&mut self) {
        for item in self.data.iter_mut() {
            item.clear();
            // *item = (0..self.length_three as usize).into_iter().map(|_| { random::<f32>() * self.length_three }).collect::<Vec<_>>();
            item.sort_by(|a,b| a.partial_cmp(b).unwrap() );
        }
    }

    /// Same run function as for `Lattice`.
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
        {
            let mut num_lengths = 0.;
            for item in &self.data {
                num_lengths += item.len() as f32;
            }
            num_lengths /= self.data.len() as f32;
            println!("Mag {} num_lengths {}", avg_mag, num_lengths);
        }
        Report {
            magnetization: avg_mag,
            susceptibility: sus * beta,
        }
    }

    /// Flips exactly one cluster. This code has been optimized to make it branchless. 
    fn evolve(&mut self, beta: f32, count: usize) {
        for _ in 0..count {
            // 1. Choose random site
            let start_index = random::<usize>() % self.data.len();
            let start_height = random::<f32>() * self.length_three;
            
            // Set up the stack
            self.cluster_stack[0] = (start_index, start_height);
            let mut stack_pointer: i32 = 0;
            let mut marked_sites = vec![false; self.data.len()];

            while stack_pointer >= 0 {
                let (my_index, my_height) = self.cluster_stack[stack_pointer as usize];
                let mut new_stack_pointer = stack_pointer - 1;
                let my_column = &mut self.data[my_index];

                // 2. Get dl and dr
                let (my_spin, dl, dr, left_index, right_index) = {
                    let mut interface_iter = my_column.iter().enumerate();
                    let mut my_spin = SPIN_BC;
                    loop {
                        match interface_iter.next() {
                            Some((i, interface)) => {
                                if *interface > my_height {
                                    if i > 0 {
                                        break (my_spin,
                                            my_height - my_column[i-1],
                                            *interface - my_height,
                                            Some(i-1), Some(i)
                                        );
                                    } else {
                                        break (my_spin,
                                            my_height,
                                            *interface - my_height,
                                            None, Some(i)
                                        );
                                    }
                                }
                            },
                            None => {
                                let zero_index = if my_column.len() == 0 { None } else { Some(my_column.len() - 1)};
                                break (my_spin,
                                    my_height - match my_column.last() {
                                        Some(f) => *f,
                                        None => 0.,
                                    },
                                    self.length_three - my_height,
                                    zero_index, None
                                );
                            }
                        }
                        my_spin = !my_spin;
                    }
                };

                // 3. Flip spins in this column
                let mut lr = -random::<f32>().ln();
                let mut ll = -random::<f32>().ln();

                for (site, height) in &self.cluster_stack[0..stack_pointer as usize] {
                    // Do not allow the flipping of points in the queue
                    if *site != my_index { continue; }
                    if *height > my_height && *height < my_height + lr {
                        lr = height - my_height - EPSILON;
                    }
                    if *height < my_height && *height > my_height - ll {
                        ll = my_height - height - EPSILON;
                    }
                }

                if lr < dr {
                    // Add new interface to the right
                    match right_index {
                        Some(i) => my_column.insert(i, my_height + lr),
                        None => my_column.push(my_height + lr)
                    }
                } else {
                    // Delete right interface
                    match right_index {
                        Some(i) => {my_column.remove(i); ()},
                        None => my_column.push(self.length_three)
                    }
                }
                if ll < dl {
                    // Add new interface to the left
                    match left_index {
                        Some(i) => my_column.insert(i + 1, my_height - ll),
                        None => my_column.insert(0, my_height - ll)
                    }
                } else {
                    // Delete left interface
                    match left_index {
                        Some(i) => { my_column.remove(i); () },
                        None => my_column.insert(0, 0.)
                    }
                }
                // assert_eq!(verify_sorted(my_column), true);
                let cr = dr.min(lr);
                let cl = dl.min(ll);

                // 4. Make bridges
                for (neighbor, strength) in self.neighbors[my_index].iter_mut() {
                    if marked_sites[*neighbor] { continue; }
                    let neighbor_column = &self.data[*neighbor];
                    let mut cluster_done_so_far = 0.;
                    loop {
                        let lx = -random::<f32>().ln() / (2. * strength * beta);
                        cluster_done_so_far += lx;
                        let neighbor_height = my_height - cl + cluster_done_so_far;
                        if cluster_done_so_far > cl + cr { break; }

                        // Get the spin of the neighbor at this point
                        let mut neighbor_spin = SPIN_BC;
                        for interface in neighbor_column {
                            if *interface > neighbor_height {
                                break;
                            }
                            neighbor_spin = !neighbor_spin;
                        }
                        if neighbor_spin != my_spin { continue; }

                        if new_stack_pointer + 1 >= CLUSTER_STACK_MAX_SIZE as i32 {
                            break;
                        }

                        // Write the neighbor height into the stack
                        self.cluster_stack[(new_stack_pointer + 1) as usize] = (*neighbor, neighbor_height);
                        marked_sites[*neighbor] = true;
                        new_stack_pointer += 1;
                    }
                }
                if new_stack_pointer + 1 >= CLUSTER_STACK_MAX_SIZE as i32 {
                    println!("Stack exceeded");
                    break;
                }
                stack_pointer = new_stack_pointer;
            }
        }
    }

    /// Computes the magnetization of the crystal
    pub fn magnetization(&self) -> f32 {
        let mut frac = 0.0;
        for item in &self.data {
            let mut spin_up = SPIN_BC;
            let mut last_interface = 0.;
            for interface in item {
                if spin_up {
                    frac += interface - last_interface;
                } else {
                    frac -= interface - last_interface;
                }
                last_interface = *interface;
                spin_up = !spin_up;
            }
            if spin_up {
                frac += self.length_three - last_interface;
            } else {
                frac -= self.length_three - last_interface;
            }
        }
        frac.abs() / self.data.len() as f32 / self.length_three
    }
}