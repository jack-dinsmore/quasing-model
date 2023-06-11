#![allow(dead_code)]
#![allow(unused_imports)]
mod lattice;
mod quantum;
mod spin;
mod funcs;
use core::num;
use std::{fs::File, io::Write, thread, sync::Arc};

use lattice::Lattice;
use quantum::QLattice;
use spin::{Ising, XY, Heisenberg, Spin};
use funcs::{square_fn, load_penrose, load_einstein, rect_fn};

const NUM_THREADS: usize = 8;

struct Data {
    betas: Vec<f32>,
    magnetizations: Vec<f32>,
    susceptibilities: Vec<f32>,
}

impl Data {
    fn save(&self, filename: &str) {
        let mut f = File::create(format!("data/output/{}.dat", filename)).unwrap();
        for entry in &self.betas {
            write!(&mut f, "{},", entry).unwrap();
        }
        writeln!(f,"").unwrap();
        for entry in &self.magnetizations {
            write!(&mut f, "{},", entry).unwrap();
        }
        writeln!(f,"").unwrap();
        for entry in &self.susceptibilities {
            write!(&mut f, "{},", entry).unwrap();
        }
        writeln!(f,"").unwrap();
    }
}

fn linspace(start: f32, end: f32, count: usize) -> Vec<f32> {
    if count == 1 {
        return vec![start];
    }
    let delta = (end - start) / (count - 1) as f32;
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        out.push(start + delta * i as f32);
    }
    out
}

fn reciprocal_linspace(start: f32, end: f32, count: usize) -> Vec<f32> {
    linspace(end, start, count).into_iter().map(|f|{1./f}).collect::<Vec<_>>()
}

fn get_t1_from_eta(eta: f32) -> f32 {
    ((-eta + 1.) / 4. * std::f32::consts::PI).tan()
}

fn get_t2_from_eta(eta: f32) -> f32 {
    ((eta + 1.) / 4. * std::f32::consts::PI).tan()
}

fn linspace_ex(start: f32, end: f32, count: usize) -> Vec<f32> {
    if count == 1 {
        return vec![start];
    }
    let delta = (end - start) / (count + 1) as f32;
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        out.push(start + delta * (i + 1) as f32);
    }
    out
}

fn one_pass<S: Spin>(lattice: &mut Lattice<S>, start_temp: f32, end_temp: f32,
    n_trials: usize, num_betas: usize) -> Data {

    let mut susceptibilities = Vec::new();
    let mut magnetizations = Vec::new();
    let betas = reciprocal_linspace(start_temp, end_temp, num_betas);
    for beta in &betas {
        lattice.zero();
        let result = lattice.run(*beta, n_trials, (n_trials / 10, 0.5));
        magnetizations.push(result.magnetization);
        susceptibilities.push(result.susceptibility);
    }
    Data {
        betas, magnetizations, susceptibilities
    }
}

fn qone_pass(lattice: &mut QLattice, start_temp: f32, end_temp: f32,
    n_trials: usize, num_betas: usize) -> Data {

    let mut susceptibilities = Vec::new();
    let mut magnetizations = Vec::new();
    let betas = reciprocal_linspace(start_temp, end_temp, num_betas);
    for beta in &betas {
        print!("Temp {} \t", 1./beta);
        lattice.zero();
        let result = lattice.run(*beta, n_trials, (n_trials / 10, 0.5));
        magnetizations.push(result.magnetization);
        susceptibilities.push(result.susceptibility);
    }
    Data {
        betas, magnetizations, susceptibilities
    }
}

fn search<S: Spin>(lattice: &mut Lattice<S>, bottom: f32, top: f32,
    layers: usize, n_trials: usize, count_per_iteration: usize) -> Data {

    let mut betas = Vec::new();
    let mut susceptibilities = Vec::new();
    let mut magnetizations = Vec::new();
    let mut start = bottom;
    let mut end = top;
    for iter_count in 0..layers {
        let mut new_start = start;
        let mut new_end = end;
        let beta_line = if iter_count == 0 {
            linspace(start, end, count_per_iteration)
        } else {
            linspace_ex(start, end, count_per_iteration)
        };
        for beta in beta_line {
            lattice.zero();
            let result = lattice.run(beta, n_trials, (n_trials / 2, 0.5));
            magnetizations.push(result.magnetization);
            susceptibilities.push(result.susceptibility);
            betas.push(beta);

            if result.magnetization < 0.02 {
                new_start = new_start.max(beta);
            }
            if result.magnetization > 0.2 {
                new_end = new_end.min(beta);
            }

        }
        start = new_start;
        end = new_end;
    }

    Data {
        betas, magnetizations, susceptibilities
    }
}

fn main() {
    // one();
    rect::<Ising>(6.);
    einstein::<Ising>(0.8);
    rect::<XY>(1.3);
    einstein::<XY>(1.5);
    rect::<Heisenberg>(0.9);
    einstein::<Heisenberg>(0.4);
}

fn one() {
    // let (size, func) = (12*12, square_fn(12));
    // println!("{} sites", size);
    // let mut lattice = QLattice::new(size, &func);
    // // let data = qone_pass(&mut lattice, 0.01, 0.01, 1000, 1);
    // let data = qone_pass(&mut lattice, 0.1, 2., 1000, 10);
    // data.save("tim-square");
    
    // let (size, func) = (128*128, square_fn(128));
    // println!("{} sites", size);
    // let mut lattice = Lattice::<Ising>::new(size, &func);
    // let data = one_pass(&mut lattice, 0.01, 4., 10000, 50);
    // data.save("ising-square");

    // let mut lattice = Lattice::<XY>::new(size, &func);
    // let data = one_pass(&mut lattice, 0.01, 4., 10000, 50);
    // data.save("xy-square");

    // let mut lattice = Lattice::<Heisenberg>::new(size, &func);
    // let data = one_pass(&mut lattice, 0.01, 4., 10000, 50);
    // data.save("heisenberg-square");

    let (size, func) = load_penrose(9);
    println!("{} sites", size);
    let mut lattice = Lattice::<Ising>::new(size, &func);
    let data = one_pass(&mut lattice, 0.01, 4., 10000, 50);
    data.save("ising-penrose");

    let mut lattice = Lattice::<XY>::new(size, &func);
    let data = one_pass(&mut lattice, 0.01, 4., 10000, 50);
    data.save("xy-penrose");

    let mut lattice = Lattice::<Heisenberg>::new(size, &func);
    let data = one_pass(&mut lattice, 0.01, 4., 10000, 50);
    data.save("heisenberg-penrose");
}

fn rect<S: Spin>(tmax: f32) {
    println!("{} rect", S::name());
    let num_per_thread = 5;
    let size = 128;
    let mut threads = Vec::new();
    let etas = linspace(-0.99,0.99, NUM_THREADS * num_per_thread);

    for thread_index in 0..NUM_THREADS {
        let mut eta_chunk = Vec::new();
        for eta_index in 0..num_per_thread {
            let i = num_per_thread * thread_index + eta_index;
            eta_chunk.push(etas[i]);
        }

        threads.push(thread::spawn(move || {
            for eta in eta_chunk {
                println!("{}", eta);
                let func = rect_fn(size, get_t1_from_eta(eta), get_t2_from_eta(eta));
                let mut lattice = Lattice::<S>::new(size*size, &func);
                let data = one_pass(&mut lattice, 0.01, tmax, 10000, 50);
                data.save(&format!("rect-{}-{:.8}", S::name(), eta));
            }
        }));
    }
    
    for thread in threads {
        thread.join().unwrap();
    }
}

fn einstein<S: Spin>(tmax: f32) {
    println!("{} einstein", S::name());
    let num_per_thread = 5;
    let mut threads = Vec::new();
    let etas = linspace(-0.99,0.99, NUM_THREADS * num_per_thread);

    for thread_index in 0..NUM_THREADS {
        let mut eta_chunk = Vec::new();
        for eta_index in 0..num_per_thread {
            let i = num_per_thread * thread_index + eta_index;
            eta_chunk.push(etas[i]);
        }

        threads.push(thread::spawn(move || {
            for eta in eta_chunk {
                println!("{}", eta);
                let (size, func) = load_einstein("7k", get_t1_from_eta(eta), get_t2_from_eta(eta));
                let mut lattice = Lattice::<S>::new(size, &func);
                let data = one_pass(&mut lattice, 0.01, tmax, 10000, 50);
                data.save(&format!("einstein-{}-{:.8}", S::name(), eta));
            }
        }));
    }
    
    for thread in threads {
        thread.join().unwrap();
    }
}