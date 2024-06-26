use std::f32::consts::{E, PI};
use rand::prelude::*;

const SAMPLE_SIZE: usize = 1_000_000_00;

fn accept(a:f32,b:f32) -> Option<f32> {
    let x = b/a;
    if a.powf(2.0) < (2.0*PI).powf(-0.5)*E.powf(-0.5*x.powf(2.0)) {
        Some(x)
    } else {
        None
    }
}

fn main() {

    let a = (2.0/E).powf(0.5)*(2.0*PI).powf(-0.25);
    let b = ((1.0/(2.0*PI).powf(0.5))).powf(0.5);

    let mut rng = rand::thread_rng();
    let arr1:Vec<f32> = (0..SAMPLE_SIZE).map(|_| rng.gen_range(0.0..b)).collect();
    let arr2:Vec<f32> = (0..SAMPLE_SIZE).map(|_| rng.gen_range(-a..a)).collect();
    let sample:Vec<f32> = arr1.iter()
        .zip(arr2.iter())
        .filter_map(|(a,b)| accept(*a,*b))
        .collect();
    println!("{:?}",sample.len());
}
