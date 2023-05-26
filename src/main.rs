#![allow(dead_code)]
mod lattice;
mod spin;
mod funcs;
use lattice::Lattice;
use spin::Ising;
use funcs::{square, SQUARE_ROW_SIZE};

fn linspace(start: f64, end: f64, count: usize) -> Vec<f64> {
    if count == 1 {
        return vec![start];
    }
    let delta = (end - start) / (count - 1) as f64;
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        out.push(start + delta * i as f64);
    }
    out
}

fn main() {
    // let mut lattice = Lattice::<Ising>::new(
    //     SQUARE_ROW_SIZE*SQUARE_ROW_SIZE, square
    // );
    let mut lattice = Lattice::<Ising>::new(
        SQUARE_ROW_SIZE*SQUARE_ROW_SIZE, square
    );

    let n_trials = 1000;
    let mut order_params = Vec::new();
    let mut susceptibilities = Vec::new();
    let betas = linspace(0.4, 0.5, 11);
    for beta in &betas {
        let result = lattice.run(*beta, n_trials);
        order_params.push(result.order_parameter);
        susceptibilities.push(result.susceptibility);
        println!("Beta {}: {}", beta, result);
        lattice.zero();
    }

    println!("{:?}", betas);
    println!("{:?}", order_params);
    println!("{:?}", susceptibilities);
}