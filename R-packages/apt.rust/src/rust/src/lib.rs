use extendr_api::prelude::*;
use std::f64::consts::{E, PI};

use array_init::array_init;
use rand::prelude::*;
use rand_distr::StandardNormal;

const P:f64 = 0.3;
const MU1:f64 = 1.0;
const SD1: f64 = 1.0;
const MU2:f64 = 20.0;
const SD2:f64 = 4.0;

const L:usize = 5;
const ASTAR:f64 = 0.234;



fn normal_density(x:f64,mean:f64,sd:f64)->f64 {
    1.0/(2.0*PI*sd.powf(2.0)).powf(0.5)*E.powf(-0.5*((x-mean)/sd).powf(2.0))
}

fn density(x:f64,inv_temp:f64)->f64 {
    (P*normal_density(x, MU1, SD1)+(1.0-P)*normal_density(x, MU2, SD2)).powf(inv_temp)
}

#[extendr]
fn apt(burn:i32,ss:i32,adapt:bool)-> Vec<f64> {

    let burn = burn as usize;
    let ss = ss as usize;

    let mut rng = thread_rng();

    let mut rho: [f64;L-1] = [1.0;L-1];
    let mut beta: [f64;L] = array_init(|i| 0.5_f64.powf(i as f64));
    let mut x: [f64;L] = [0.0;L];
    let mut swap_ps: [f64;L-1] = [0.0;L-1];

    let mut var: [f64;L] = [1.0;L];
    let mut var_rho: [f64;L] = [1.0;L];
    let mut mu_var: f64 = 0.0;
    let mut t_var: [f64;L] = [1.0;L];



    let mut samples: Vec<f64> = Vec::new();

    let mut n:usize = 0;

    while n < burn+ss {
        n+=1;
        let gamma:f64 = (n as f64).powf(-0.6);

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
            let candidate:f64 = var[i].powf(0.5)*rng.sample::<f64,_>(StandardNormal) + x[i];
            let accept_p: f64 = 1.0_f64.min(density(candidate, beta[i])/density(x[i], beta[i]));
            if rng.gen::<f64>()<accept_p {
                x[i]=  candidate;
            }

            // update variance
            //get mean of elements:
            let mut mean_x: f64 = 0.0;
            for y in 0..L {
              mean_x+=x[y];
            }
            mean_x = mean_x/(L as f64);
            mu_var = (1.0-gamma*0.5)*mu_var + gamma*0.5*mean_x;

            let mut var_x: f64 = 0.0;
            for y in 0..L {
              var_x+=(x[y]-mu_var).powf(2.0);
            }

            var_rho[i] = (1.0-gamma*0.5)*var_rho[i] + gamma*0.5*1.0/(L as f64)*var_x;
            t_var[i]+=gamma*(accept_p-ASTAR);
            var[i] = E.powf(t_var[i])*var_rho[i]

        }

        if n>burn {
            samples.push(x[0]);
        }


        if adapt {
            // Temp update


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

    }


    samples
}

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod apt_rust;
    fn hello_world;
    fn apt;
}
