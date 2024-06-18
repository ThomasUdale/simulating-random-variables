use std::f64::consts::{E, PI};

use array_init::array_init;
use indicatif::ProgressBar;
use rand::prelude::*;
use rand_distr::StandardNormal;

const P:f64 = 0.3;
const MU1:f64 = 1.0;
const SD1: f64 = 1.0;
const MU2:f64 = 20.0;
const SD2:f64 = 4.0;

const BURN:usize = 1_000_000;
const SS:usize = 1_000_000;
const L:usize = 5;
const ASTAR:f64 = 0.234;



fn normal_density(x:f64,mean:f64,sd:f64)->f64 {
    1.0/(2.0*PI*sd.powf(2.0)).powf(0.5)*E.powf(-0.5*((x-mean)/sd).powf(2.0))
}

fn density(x:f64,inv_temp:f64)->f64 {
    (P*normal_density(x, MU1, SD1)+(1.0-P)*normal_density(x, MU2, SD2)).powf(inv_temp)
}

fn apt(burn:usize,ss:usize)-> Vec<f64> {
    let mut rng = thread_rng();
    
    let mut rho: [f64;L-1] = [1.0;L-1];
    let mut beta: [f64;L] = array_init(|i| 0.5_f64.powf(i as f64));
    let mut x: [f64;L] = [0.0;L];
    let mut swap_ps: [f64;L-1] = [0.0;L-1];

    let mut samples: Vec<f64> = Vec::new();

    let mut n:usize = 0;

    let bar = ProgressBar::new((burn+ss) as u64);

    while n < burn+ss {
        bar.inc(1);
        n+=1;

        // Swap Probabilities:
        for i in 0..(L-1) {
            swap_ps[i] = 1.0_f64.min(
                (density(x[i], beta[i+1])*density(x[i+1], beta[i]))
                /(density(x[i], beta[i])*density(x[i+1], beta[i+1]))
            );
        }

        let swap_candidate: usize = rng.gen_range(0..L-1);
        if rng.gen::<f64>()<swap_ps[swap_candidate] {
            x.swap(swap_candidate, swap_candidate+1)
        }

        // MH
        for i in 0..L {
            let candidate:f64 = rng.sample::<f64,_>(StandardNormal) + x[i];
            if rng.gen::<f64>()<1.0_f64.min(density(candidate, beta[i])/density(x[i], beta[i])){
                x[i]=  candidate;
            }
        }

        if n>burn {
            samples.push(x[0]);
        }
        


        // Temp update
        let gamma:f64 = (n as f64).powf(-0.6);

        let mut t:f64 = 1.0;
        for l in 0..L {
            if l<L-1 {
                rho[l] = rho[l] + gamma*(swap_ps[l]-ASTAR);
            }
            if l>0 {
                t+=E.powf(rho[l-1]);
                beta[l] = 1.0/t;
            }
        }
    }

    bar.finish();

    samples
}


fn main() {

    let samples = apt(BURN,SS);
}
