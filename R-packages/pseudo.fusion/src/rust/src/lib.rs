use extendr_api::prelude::*;

use rand::{distributions::Distribution, rngs::ThreadRng, thread_rng};
use statrs::distribution::{Cauchy, Continuous, Normal, StudentsT, Uniform};

const NUM_DIST: i32 = 3;

enum Dists {
    C(rand_distr::Cauchy<f64>,Cauchy),
    N(rand_distr::Normal<f64>,Normal),
    T(rand_distr::StudentT<f64>,StudentsT)
}

trait A { fn sample(&self,rng:&mut ThreadRng)->f64;  fn pdf(&self,x:f64)->f64;}
impl A for Dists {
    fn sample(&self,rng:&mut ThreadRng) -> f64 {
        match &self {
            Dists::C(s1,_) => s1.sample(rng),
            Dists::N(s1,_) => s1.sample(rng),
            Dists::T(s1,_) => s1.sample(rng)
        }
    }

    fn pdf(&self,x:f64) -> f64 {
        match &self {
            Dists::C(_,s1) => s1.pdf(x),
            Dists::N(_,s1) => s1.pdf(x),
            Dists::T(_,s1) => s1.pdf(x),
        }
    }
}

#[extendr]
fn mixture_fusion_example(burn:i32,ss:i32,num_is:i32) -> Vec<f64> {


    let burn = burn as usize;
    let ss = ss as usize;
    let num_is = num_is as usize;

    let mut rng = thread_rng();


    let u = Uniform::new(0.0,1.0).unwrap();

    let mut fs:Vec<Dists> = Vec::new();
    fs.push(Dists::C(rand_distr::Cauchy::new(-2.0, 1.0).unwrap(),Cauchy::new(-2.0, 1.0).unwrap()));
    fs.push(Dists::N(rand_distr::Normal::new(2.0,5.0).unwrap(),Normal::new(2.0,5.0).unwrap()));
    fs.push(Dists::T(rand_distr::StudentT::new(1.0).unwrap(),StudentsT::new(0.0,1.0,1.0).unwrap()));

    let mut samples: Vec<f64> = Vec::new();
    let mut n:usize = 0;
    let mut y:f64 = 0.0;
    let mut r:f64 = 0.000000001;


    while n<burn+ss {
        n+=1;

        let q1 = rand_distr::Normal::new(y,1.0).unwrap();
        let y_new = q1.sample(&mut rng);

        let mut new_r = 0.0;
        for _ in 0..num_is {
            let mut r_star = 1.0;
            for f in fs.iter() {
                let x = f.sample(&mut rng);
                let q2 = Normal::new(x,1.0).unwrap();
                r_star*=1.0_f64.min(f.pdf(y_new)/f.pdf(x))*q2.pdf(y_new);
            }
            new_r+=r_star/NUM_DIST as f64;
        }

        if u.sample(&mut rng) < new_r/r {
            y = y_new;
            r = new_r;
        }
        samples.push(y);

    }

    samples

}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod pseudo_fusion;
    fn mixture_fusion_example;
}
