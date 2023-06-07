use std::ops::AddAssign;
use std::fmt::Debug;

use rand::prelude::*;
use rand_distr::StandardNormal;

const MAX_NEIGHBORS: usize = 12;

#[derive(Debug)]
pub struct Ising {
    data: i32,
}

#[derive(Debug)]
pub struct XY {
    x: f32,
    y: f32,
}

impl<'a> AddAssign<&'a Self> for Ising {
    fn add_assign(&mut self, rhs: &'a Self) {
        self.data += rhs.data;
    }
}

impl<'a> AddAssign<&'a Self> for XY {
    fn add_assign(&mut self, rhs: &'a Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Spin for Ising {
    type V = u8;
    fn dot(&mut self, _vec: &Self::V) -> f32 {
        self.data as f32
    }
    fn flip(&mut self, _vec: &Self::V, success: bool) {
        self.data *= 1 - 2 * (success as i32)
    }
    fn norm(&self) -> f32 {
        (self.data as f32).abs()
    }
    fn random_vec() -> Self::V {
        0
    }
    fn start() -> Self {
        Self { data: if random::<bool>() {1} else {-1} }
    }
    fn zero() -> Self {
        Self { data: 0 }
    }
}

impl Spin for XY {
    type V = (f32, f32);
    fn dot(&mut self, vec: &Self::V) -> f32 {
        self.x * vec.0 + self.y * vec.1
    }
    fn flip(&mut self, vec: &Self::V, success: bool) {
        let dot = (vec.0 * self.x + vec.1 * self.y) * (success as i32 as f32);
        self.x -= 2. * dot * vec.0;
        self.y -= 2. * dot * vec.1;
        let norm = 1./(self.x*self.x + self.y*self.y).sqrt();
        self.x *= norm;
        self.y *= norm;
    }
    fn norm(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
    fn random_vec() -> Self::V {
        let x = 0.1 * thread_rng().sample::<f32,_>(StandardNormal);
        let y = 0.1 * thread_rng().sample::<f32,_>(StandardNormal);
        let norm = 1./(x*x + y*y).sqrt();
        (x * norm, y * norm)
    }
    fn start() -> Self {
        let (x, y) = Self::random_vec();
        Self { x, y }
    }
    fn zero() -> Self {
        Self { x: 0., y: 0. }
    }
}

#[derive(Debug)]
pub struct SmallVec<T> {
    data: [Option<T>; MAX_NEIGHBORS]
}

pub struct SmallVecIter<'a, T> {
    origin: &'a SmallVec<T>,
    index: usize,
}

impl<'a, T> Iterator for SmallVecIter<'a, T> {
    type Item = &'a T;
    fn next(&mut self) -> Option<Self::Item> {
        let val = self.origin.data[self.index].as_ref();
        self.index += 1;
        val
    }
}

impl<T: Copy> SmallVec<T> {
    pub fn new() -> Self {
        Self {
            data: [None; MAX_NEIGHBORS],
        }
    }

    pub fn push(&mut self, val: T) {
        for i in 0..MAX_NEIGHBORS {
            if let None = self.data[i] {
                self.data[i] = Some(val);
                return;
            }
        }
        panic!("Could not push because the list was full")
    }

    pub fn iter<'a>(&'a self) -> SmallVecIter<'a, T> {
        SmallVecIter {
            origin: &self,
            index: 0,
        }
    }

    pub fn iter_mut<'a>(&'a mut self) -> SmallVecIter<'a, T> {
        SmallVecIter {
            origin: self,
            index: 0,
        }
    }
}

pub trait Spin: for <'a> AddAssign<&'a Self> + Sized + Debug {
    /// Type of the seed vector used to flip
    type V;
    /// Flip the spin value
    fn flip(&mut self, vec: &Self::V, success: bool);
    /// Dot two spin values against the other
    fn dot(&mut self, vec: &Self::V) -> f32;
    /// Generate a random seed vector to flip
    fn random_vec() -> Self::V;
    /// Get a zero-valued spin for the sake of averaging
    fn zero() -> Self;
    /// Generate a random spin value to start
    fn start() -> Self;
    /// Get the norm for the sake of averaging
    fn norm(&self) -> f32;
}