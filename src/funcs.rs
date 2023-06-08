use std::io::Read;

use npy::NpyData;

use crate::spin::SmallVec;

pub fn square_fn(row_size: usize) -> impl Fn(usize) -> SmallVec<(usize, f32)> {
    rect_fn(row_size, 1.)
}

pub fn rect_fn(row_size: usize, t2: f32) -> impl Fn(usize) -> SmallVec<(usize, f32)> {
    move |site: usize| {
        let mut neighbors = SmallVec::new();
        if site % row_size != 0 {
            neighbors.push((site - 1, 1.))
        } else {
            neighbors.push((site + row_size - 1, 1.))
        }
        if site % row_size != row_size - 1 {
            neighbors.push((site + 1, 1.))
        } else {
            neighbors.push((site - row_size + 1, 1.))
        }
        if site / row_size != 0 {
            neighbors.push((site - row_size, t2))
        } else {
            neighbors.push((row_size * (row_size-1) + site, t2))
        }
        if site / row_size != row_size - 1 {
            neighbors.push((site + row_size, t2))
        } else {
            neighbors.push((site - row_size * (site / row_size), t2))
        }
        neighbors
    }
}

pub fn load_penrose(level: usize) -> (usize, impl Fn(usize) ->SmallVec<(usize, f32)>) {
    let mut buf = vec![];
    std::fs::File::open(format!("data/penrose-{}.npy", level)).unwrap()
        .read_to_end(&mut buf).unwrap();
    let bond_indices = NpyData::<i64>::from_bytes(&buf).unwrap().to_vec();
    let mut max_bond_index = 0;
    for i in &bond_indices {
        max_bond_index = max_bond_index.max(*i);
    }
    (
        max_bond_index as usize + 1,
        move |site: usize| {
            let mut neighbors = SmallVec::new();
            for pair in bond_indices.chunks(2) {
                if pair[0] as usize == site {
                    neighbors.push((pair[1] as usize, 1.0));
                }
                if pair[1] as usize == site {
                    neighbors.push((pair[0] as usize, 1.0));
                }
            }
            neighbors
        }
    )
}

pub fn load_einstein(name: &str, t2: f32) -> (usize, impl Fn(usize) ->SmallVec<(usize, f32)>) {
    let mut small_buf = vec![];
    let mut med_buf = vec![];
    let mut long_buf = vec![];
    std::fs::File::open(format!("data/einstein-{}-short.npy", name)).unwrap()
        .read_to_end(&mut small_buf).unwrap();
    let small_indices = NpyData::<i64>::from_bytes(&small_buf).unwrap().to_vec();
    
    std::fs::File::open(format!("data/einstein-{}-medium.npy", name)).unwrap()
        .read_to_end(&mut med_buf).unwrap();
    let medium_indices = NpyData::<i64>::from_bytes(&med_buf).unwrap().to_vec();

    std::fs::File::open(format!("data/einstein-{}-long.npy", name)).unwrap()
    .read_to_end(&mut long_buf).unwrap();
    let large_indices = NpyData::<i64>::from_bytes(&long_buf).unwrap().to_vec();

    let size = small_indices.iter().max().max(
        medium_indices.iter().max().max(
            large_indices.iter().max())
    ).unwrap();

    (
        *size as usize + 1,
        move |site: usize| {
            let mut neighbors = SmallVec::new();
            for pair in small_indices.chunks(2) {
                if pair[0] as usize == site {
                    neighbors.push((pair[1] as usize, 1.0));
                }
                if pair[1] as usize == site {
                    neighbors.push((pair[0] as usize, 1.0));
                }
            }
            for pair in medium_indices.chunks(2) {
                if pair[0] as usize == site {
                    neighbors.push((pair[1] as usize, t2));
                }
                if pair[1] as usize == site {
                    neighbors.push((pair[0] as usize, t2));
                }
            }
            for pair in large_indices.chunks(2) {
                if pair[0] as usize == site {
                    neighbors.push((pair[1] as usize, 0.25));
                }
                if pair[1] as usize == site {
                    neighbors.push((pair[0] as usize, 0.25));
                }
            }
            neighbors
        }
    )
}