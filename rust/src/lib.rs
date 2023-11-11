pub mod circ;

pub use linen::linen_circuit;
pub use rayon::rayon_circuit;

mod linen;
mod rayon;
