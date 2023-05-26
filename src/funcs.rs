use crate::spin::SmallVec;

pub const SQUARE_ROW_SIZE: usize = 128;

pub fn square(site: usize) -> SmallVec<usize> {
    let mut neighbors = SmallVec::new();
    if site % SQUARE_ROW_SIZE != 0 {
        neighbors.push(site - 1)
    } else {
        neighbors.push(site + SQUARE_ROW_SIZE - 1)
    }
    if site % SQUARE_ROW_SIZE != SQUARE_ROW_SIZE - 1 {
        neighbors.push(site + 1)
    } else {
        neighbors.push(site - SQUARE_ROW_SIZE + 1)
    }
    if site / SQUARE_ROW_SIZE != 0 {
        neighbors.push(site - SQUARE_ROW_SIZE)
    } else {
        neighbors.push(SQUARE_ROW_SIZE * (SQUARE_ROW_SIZE-1) + site)
    }
    if site / SQUARE_ROW_SIZE != SQUARE_ROW_SIZE - 1 {
        neighbors.push(site + SQUARE_ROW_SIZE)
    } else {
        neighbors.push(site - SQUARE_ROW_SIZE * (site / SQUARE_ROW_SIZE))
    }
    neighbors
}