#![allow(dead_code)]
mod lattice;
mod spin;
use lattice::Lattice;
use spin::{Ising, SmallVec};

fn main() {
    let row_size = 100;
    let mut square_periodic_lattice = Lattice::<Ising>::new(row_size*row_size,
        |site: usize| -> SmallVec<usize> {
            let mut neighbors = SmallVec::new();
            if site % row_size != 0 {
                neighbors.push(site - 1)
            } else {
                neighbors.push(site + row_size - 1)
            }
            if site % row_size != row_size - 1 {
                neighbors.push(site + 1)
            } else {
                neighbors.push(site - row_size + 1)
            }
            if site / row_size != 0 {
                neighbors.push(site - row_size)
            } else {
                neighbors.push(row_size * (row_size-1) + site)
            }
            if site / row_size != row_size - 1 {
                neighbors.push(site + row_size)
            } else {
                neighbors.push(site - row_size * (site / row_size))
            }
            neighbors
        }
    );

    let n_trials = 1000;
    for beta in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0] {
        println!("Beta {}, param {}", beta, square_periodic_lattice.run(beta, n_trials));
        square_periodic_lattice.zero();
    }
}