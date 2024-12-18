use extendr_api::prelude::*;
use rand::{distributions::Distribution, rngs::ThreadRng, thread_rng, distributions::WeightedIndex};
use statrs::distribution::{Normal, Uniform,Exp,DiscreteUniform, Poisson,Continuous,Gamma};


/*rextendr::document("~/Documents/simulating-random-variables/R-packages/bm.rust")
devtools::load_all("~/Documents/simulating-random-variables/R-packages/bm.rust")*/


/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

#[extendr]
fn brownian_bridge(x:f64,y:f64,s:f64,t:f64,times:Vec<f64>) -> Vec<f64> {
  let mut rng = thread_rng();
  let std_n = Normal::new(0.0, 1.0).unwrap();
  let mut bm: Vec<f64> = Vec::new();
  for tid in 0..(times.len()+1) {
    if tid==0{
      bm.push(std_n.sample(&mut rng)*(times[tid]-s).powf(0.5));
    } else if tid==times.len() {
      bm.push(std_n.sample(&mut rng)*(t-times[tid-1]).powf(0.5) + bm[tid-1]);
    } else {
      bm.push(std_n.sample(&mut rng)*(times[tid]-times[tid-1]).powf(0.5) + bm[tid-1]);
    }
  }

  for tid in 0..times.len() {
    bm[tid] = bm[tid] - (times[tid]/t)*bm[times.len()] + (1.0-(times[tid]/t))*x + (times[tid]/t)*y;
  }
  bm[0..times.len()].to_vec()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod bm_rust;
    fn hello_world;
    fn brownian_bridge;
}
