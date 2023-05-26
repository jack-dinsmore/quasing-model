use std::ops::AddAssign;
use std::fmt::Debug;

use rand::random;

const MAX_NEIGHBORS: usize = 8;

#[derive(Debug)]
pub struct Ising {
    data: i32,
}

impl<'a> AddAssign<&'a Self> for Ising {
    fn add_assign(&mut self, rhs: &'a Self) {
        self.data += rhs.data;
    }
}

impl Spin for Ising {
    type V = u8;
    fn dot(&mut self, _vec: &Self::V) -> f64 {
        self.data as f64
    }
    fn flip(&mut self, _vec: &Self::V, success: bool) {
        self.data *= 1 - 2 * (success as i32)
    }
    fn norm(&self) -> f64 {
        (self.data as f64).abs()
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
    type V;
    fn flip(&mut self, vec: &Self::V, success: bool);
    fn dot(&mut self, vec: &Self::V) -> f64;
    fn random_vec() -> Self::V;
    fn zero() -> Self;
    fn start() -> Self;
    fn norm(&self) -> f64;
}