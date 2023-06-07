#![allow(dead_code)]
#![allow(unused_imports)]
mod lattice;
mod spin;
mod funcs;
use core::num;
use std::{fs::File, io::Write, thread, sync::Arc};

use lattice::Lattice;
use spin::{Ising, Spin};
use funcs::{square_fn, load_penrose, load_einstein, rect_fn};

use crate::spin::XY;

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
        // println!("{}, {}", start, end);
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

            // println!("..{}, {}, {}", beta, result.magnetization, result.susceptibility);
        }
        start = new_start;
        end = new_end;
    }

    Data {
        betas, magnetizations, susceptibilities
    }
}

fn main() {
    xy_rect();
    xy_einstein();
    ising_rect();
    ising_einstein();
    // one();
}

fn one() {
    let (size, func) = (128*128, square_fn(128));
    // let (size, func) = load_penrose(6);
    println!("{} sites", size);
    let mut lattice = Lattice::<Ising>::new(size, func);
    let data = search(&mut lattice, 0., 2., 1, 1000, 12);
    data.save("ising-square");
}

fn ising_rect() {
    let num_per_thread = 5;
    let size = 128;
    let mut threads = Vec::new();
    let t2s = linspace(0.01,3.,NUM_THREADS * num_per_thread);

    for thread_index in 0..NUM_THREADS {
        let mut t2_chunk = Vec::new();
        for t2_index in 0..num_per_thread {
            let i = num_per_thread * thread_index + t2_index;
            t2_chunk.push(t2s[i]);
        }

        threads.push(thread::spawn(move || {
            for t2 in t2_chunk {
                println!("{}", t2);
                let func = rect_fn(size, t2);
                let mut lattice = Lattice::<Ising>::new(size*size, func);
                let data = search(&mut lattice, 0., 10., 6, 1000, 12);
                data.save(&format!("rect-ising-{:.8}", t2));
            }
        }));
    }
    
    for thread in threads {
        thread.join().unwrap();
    }
}

fn xy_rect() {
    let num_per_thread = 5;
    let size = 128;
    let mut threads = Vec::new();
    let t2s = linspace(0.01,3.,NUM_THREADS * num_per_thread);

    for thread_index in 0..NUM_THREADS {
        let mut t2_chunk = Vec::new();
        for t2_index in 0..num_per_thread {
            let i = num_per_thread * thread_index + t2_index;
            t2_chunk.push(t2s[i]);
        }

        threads.push(thread::spawn(move || {
            for t2 in t2_chunk {
                println!("{}", t2);
                let func = rect_fn(size, t2);
                let mut lattice = Lattice::<XY>::new(size*size, func);
                let data = search(&mut lattice, 0., 10., 6, 3000, 12);
                data.save(&format!("rect-xy-{:.8}", t2));
            }
        }));
    }
    
    for thread in threads {
        thread.join().unwrap();
    }
}

fn ising_einstein() {
    let num_per_thread = 5;
    let mut threads = Vec::new();
    let t2s = linspace(0.01,2.,NUM_THREADS * num_per_thread);

    for thread_index in 0..NUM_THREADS {
        let mut t2_chunk = Vec::new();
        for t2_index in 0..num_per_thread {
            let i = num_per_thread * thread_index + t2_index;
            t2_chunk.push(t2s[i]);
        }

        threads.push(thread::spawn(move || {
            for t2 in t2_chunk {
                println!("{}", t2);
                let (size, func) = load_einstein("7k", t2);
                let mut lattice = Lattice::<Ising>::new(size, func);
                let data = search(&mut lattice, 0., 10., 5, 1000, 12);
                data.save(&format!("einstein-ising-{:.8}", t2));
            }
        }));
    }
    
    for thread in threads {
        thread.join().unwrap();
    }
}

fn xy_einstein() {
    let num_per_thread = 5;
    let mut threads = Vec::new();
    let t2s = linspace(0.01,2.,NUM_THREADS * num_per_thread);

    for thread_index in 0..NUM_THREADS {
        let mut t2_chunk = Vec::new();
        for t2_index in 0..num_per_thread {
            let i = num_per_thread * thread_index + t2_index;
            t2_chunk.push(t2s[i]);
        }

        threads.push(thread::spawn(move || {
            for t2 in t2_chunk {
                println!("{}", t2);
                let (size, func) = load_einstein("7k", t2);
                let mut lattice = Lattice::<XY>::new(size, func);
                let data = search(&mut lattice, 0., 10., 5, 3000, 12);
                data.save(&format!("einstein-xy-{:.8}", t2));
            }
        }));
    }
    
    for thread in threads {
        thread.join().unwrap();
    }
}