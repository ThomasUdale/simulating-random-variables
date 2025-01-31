use extendr_api::prelude::*;

use rand::{distributions::Distribution, rngs::ThreadRng, thread_rng, distributions::WeightedIndex};
use statrs::distribution::{Normal, Uniform,Exp,DiscreteUniform, Poisson,Continuous,Gamma};

/*
rextendr::document("~/Documents/simulating-random-variables/R-packages/cts.smc")
devtools::load_all("~/Documents/simulating-random-variables/R-packages/cts.smc")
path_smc_rust(3,2,2,1)
*/

#[extendr]
fn euler_approx_rust(ss:i32,target_time:f64,n_intervals:i32) -> Vec<f64> {
  let mut rng = thread_rng();
  let std = Normal::new(0.0, 1.0).unwrap();
  let ss = ss as usize;
  let n_intervals = n_intervals as usize;
  let eps = target_time / (n_intervals as f64);
  let mut sample: Vec<f64> = Vec::new();
  let mut x: f64 = 0.0;
  for id in 0..ss {
    x=0.0;
    for n in 0..n_intervals {
      x = x + x.sin()*eps + std.sample(&mut rng)*eps.powf(0.5);
    }
    sample.push(x);
  }
  sample
}

#[extendr]
fn cts_smc_rust(ss: i32, target_time: f64) -> Vec<f64> {
  let ss: usize = ss as usize;

  let mut rng = thread_rng();

  let mut samples: Vec<f64> = vec![0.0;ss];
  let mut current_sample_time: Vec<f64> = vec![0.0;ss];
  let mut weights: Vec<f64> = vec![0.0f64;ss];
  let mut final_sample: Vec<f64> = vec![0.0f64;ss];

  let mut t: f64 = 0.0;

  let k_plus: Exp = Exp::new((ss as f64)*(5.0/8.0)).unwrap();
  let k_minus: Exp = Exp::new((ss as f64)*0.5).unwrap();
  let uniform_sample = DiscreteUniform::new(0,(ss-1).try_into().unwrap()).unwrap();
  let std = Normal::new(0.0, 1.0).unwrap();
  let u = Uniform::new(0.0,1.0).unwrap();

  let mut test_x = [0usize;2];
  while t < target_time {
    let e_plus = k_plus.sample(&mut rng);
    let e_minus = k_minus.sample(&mut rng);
    let new_t = (t + e_plus.min(e_minus)).min(target_time);
    test_x[0] = uniform_sample.sample(&mut rng)  as usize;
    samples[test_x[0]] = std.sample(&mut rng)*(new_t-current_sample_time[test_x[0]]).powf(0.5) + samples[test_x[0]];
    current_sample_time[test_x[0]] = new_t;
    let phi = 0.5*((samples[test_x[0]]).sin().powf(2.0) + (samples[test_x[0]]).cos());

    if e_plus < e_minus {
      if u.sample(&mut rng) < (phi.max(0.0)/(5.0/8.0)) {
        test_x[1] = uniform_sample.sample(&mut rng)  as usize;
        if test_x[0]!=test_x[1] {
          samples[test_x[1]] = std.sample(&mut rng)*(new_t-current_sample_time[test_x[1]]).powf(0.5) + samples[test_x[1]];
          current_sample_time[test_x[1]] = new_t;
          samples[test_x[0]] =   samples[test_x[1]];
        }
      }
    } else {
      if u.sample(&mut rng) < -phi.min(0.0)/(0.5) {
        test_x[1] = uniform_sample.sample(&mut rng)  as usize;
        if test_x[0]!=test_x[1] {
          samples[test_x[1]] = std.sample(&mut rng)*(new_t-current_sample_time[test_x[1]]).powf(0.5) + samples[test_x[1]];
          current_sample_time[test_x[1]] = new_t;
          samples[test_x[1]] =   samples[test_x[0]];
        }
      }
    }

    t=new_t;
  }




  for idx in 0..ss {
    if current_sample_time[idx] < target_time {
      samples[idx] = std.sample(&mut rng)*(target_time-current_sample_time[idx]).powf(0.5) + samples[idx];
    }
    weights[idx] = (1.0-(samples[idx]).cos()).exp()
  }



  let dist = WeightedIndex::new(&weights).unwrap();

  for idx in 0..ss {
    final_sample[idx] = samples[dist.sample(&mut rng) as usize];
  }

  final_sample
}

#[extendr]
fn brownian_bridge(rng: &mut ThreadRng,x:f64,y:&f64,s:&f64,t:&f64,times:&Vec<f64>) -> Vec<f64> {
  let std_n = Normal::new(0.0, 1.0).unwrap();
  let mut bm: Vec<f64> = Vec::new();
  for tid in 0..(times.len()+1) {
    if tid==0{
      bm.push(std_n.sample(rng)*(times[tid]-s).powf(0.5));
    } else if tid==times.len() {
      bm.push(std_n.sample(rng)*(t-times[tid-1]).powf(0.5) + bm[tid-1]);
    } else {
      bm.push(std_n.sample(rng)*(times[tid]-times[tid-1]).powf(0.5) + bm[tid-1]);
    }
  }

  for tid in 0..times.len() {
    bm[tid] = bm[tid] - (times[tid]/t)*bm[times.len()] + (1.0-(times[tid]/t))*x + (times[tid]/t)*y;
  }
  bm[0..times.len()].to_vec()
}

fn bbs(rng: &mut ThreadRng,t:&Vec<f64>,path:&Vec<f64>,times:&Vec<f64>) -> Vec<f64> {
  let mut bb: Vec<f64> = Vec::new();
  let mut cur_skel: Vec<f64> = Vec::new();
  let mut cur_index: usize = 1;
  let mut cur_tid: usize = 0;
  while cur_index<t.len() && cur_tid<times.len() {
    if times[cur_tid]<t[cur_index]  {
      cur_skel.push(times[cur_tid]);
      cur_tid+=1;
    } else {
      if cur_skel.len()>0 {
        bb.extend(brownian_bridge(
        rng,
        &path[cur_index-1],
        &path[cur_index],
        &t[cur_index-1],
        &t[cur_index],
        &cur_skel
      ));
      }
      cur_skel = Vec::new();
      cur_index+=1;
    }
  }
  if cur_skel.len()>0 {
    bb.extend(brownian_bridge(
        rng,
        &path[cur_index-1],
        &path[cur_index],
        &t[cur_index-1],
        &t[cur_index],
        &cur_skel
      ));
  }
  bb
}

fn sim_end(rng: &mut ThreadRng,t:f64) -> f64 {
  let std_n = Normal::new(0.0, 1.0).unwrap();
  let unif = Uniform::new(0.0,1.0).unwrap();
  loop {
    let u = std_n.sample(rng);
    if unif.sample(rng) < (1.0-(u).cos()).exp() {
      return u;
    }
  }
}

#[extendr]
fn path_smc_rust(n:i32,n_steps:i32,tt:f64,sample_t:f64) -> List {
  let mut paths: Vec<Vec<f64>> = Vec::new();
  let mut old_paths : Vec<Vec<f64>>;
  let mut old_times : Vec<Vec<f64>>;
  let mut times: Vec<Vec<f64>> = Vec::new();
  let mut weights: Vec<f64> = vec![1.0f64;n as usize];

  let mut rng = thread_rng();
  let std_n = Normal::new(0.0, 1.0).unwrap();
  let k = Poisson::new(tt*(9.0/8.0)).unwrap();
  let u = Uniform::new(0.0,1.0).unwrap();

  for _ in 0..n {
    times.push(vec![0.0f64,tt]);
    paths.push(vec![0.0f64,sim_end(&mut rng,tt)]);
  }

  for step in 0..n_steps {
    for p in 0..n as usize {
      let mut new_x: Vec<f64> = vec![0.0f64;1];
      let mut den_cur = 1.0;
      let mut den_prop = 1.0;
      for idx in 1..(times[p].len()) {
        new_x.push(std_n.sample(&mut rng)*0.01 + paths[p][idx]);
        den_prop *= std_n.pdf((new_x[idx]-new_x[idx-1])/(times[p][idx]-times[p][idx-1]).powf(0.5));
        den_cur *= std_n.pdf((paths[p][idx]-paths[p][idx-1])/(times[p][idx]-times[p][idx-1]).powf(0.5));
      }
      den_prop *= (1.0-(new_x[paths[p].len()-1]).cos()).exp();
      den_cur *= (1.0-(paths[p][paths[p].len()-1]).cos()).exp();
      if u.sample(&mut rng) < den_prop/den_cur {
        paths[p] = new_x;
      }
    }

    for p in 0..n as usize {
      let kp = k.sample(&mut rng);
      if (kp as usize)==0 {
        weights[p] = 1.0;
      } else {
        let mut skeleton_u: Vec<f64> = (0..(kp as usize)).map(|_| u.sample(&mut rng)*tt).collect::<Vec<f64>>();
        skeleton_u.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let bb = bbs(
          &mut rng,
          &times[p],
          &paths[p],
          &skeleton_u,
        );
        weights[p] = bb.iter().fold(
          1.0,
          |res,a| res * (9.0/8.0-0.5-0.5*((a).sin().powf(2.0) + (a).cos()))/(9.0/8.0)
        );

        if step<n_steps-1 {
          let mut update_id: usize = 0;
          let mut reference_id: usize = 1;
          while update_id<skeleton_u.len() {
            if skeleton_u[update_id]>times[p][reference_id] {
              reference_id +=1;
            } else {
              times[p].insert(reference_id,skeleton_u[update_id]);
              paths[p].insert(reference_id,bb[update_id]);
              update_id+=1;
            }
          }
        }
      }
    }

    if step<n_steps-1 {
      old_paths = paths.clone();
      old_times = times.clone();
      let dist = WeightedIndex::new(&weights).unwrap();
      for p in 0..n as usize {
        let new_id = dist.sample(&mut rng) as usize;
        times[p] = old_times[new_id].clone();
        paths[p] = old_paths[new_id].clone();
      }
    }
  }

  let mut out_sample: Vec<f64> = Vec::new();
  for p in 0..n as usize {
    out_sample.push(bbs(
      &mut rng,
      &times[p],
      &paths[p],
      &vec![sample_t]
    )[0])
  }
  list!(w=out_sample,weights=weights)
}

#[extendr]
fn brownian_motion_fpt(n:i32) -> Vec<f64> {
  let mut rng = thread_rng();
  let gam = Gamma::new(1.088870,1.2336997422).unwrap();
  let u = Uniform::new(0.0,1.0).unwrap();
  let n: usize = n as usize;
  let mut passage_times: Vec<f64> = Vec::new();
  let mut it_n: i32;
  let mut fn0: f64;
  let mut fn1: f64;
  for _id in 0..n {
    loop {
      let x:f64 = gam.sample(&mut rng);
      let y:f64 = u.sample(&mut rng) * 1.243707 * gam.pdf(x);
      let sqrtx:f64 = (2.0*std::f64::consts::PI*x.powf(3.0)).powf(0.5);
      it_n = (0.275*x).ceil().max(3.0) as i32;
      fn0 = (-it_n..=it_n).map(|z| (-1.0f64).powf(z as f64)*(1.0+2.0*(z as f64))*(-(1.0+2.0*(z as f64)).powf(2.0)/(2.0*x)).exp()).sum::<f64>()/sqrtx;
      it_n+=1;
      fn1 = fn0 + (-1.0f64).powf(it_n as f64)*((1.0-2.0*(it_n as f64))*((-(1.0-2.0*(it_n as f64)).powf(2.0)/(2.0*x)).exp())+(1.0+2.0*(it_n as f64))*((-(1.0+2.0*(it_n as f64)).powf(2.0)/(2.0*x)).exp()))/sqrtx;
      while (y-fn0)*(y-fn1) < 0.0 {
        fn0 = fn1;
        it_n+=1;
        fn1 = fn0 + (-1.0f64).powf(it_n as f64)*((1.0-2.0*(it_n as f64))*(-((1.0-2.0*(it_n as f64)).powf(2.0)/(2.0*x)).exp())+(1.0+2.0*(it_n as f64))*(-((1.0+2.0*(it_n as f64)).powf(2.0)/(2.0*x)).exp()))/sqrtx;
      }
      if y<=fn1 {
        passage_times.push(x);
        break;
      }
    }
  }
  passage_times
}

fn g(k:f64,t:f64) -> f64 {
  (-1.0f64).powf(k)*(1.0+2.0*k)/(2.0*std::f64::consts::PI*t.powf(3.0)).powf(0.5)*(-(1.0+2.0*k).powf(2.0)/(2.0*t)).exp()
}

#[extendr]
fn brownian_motion_fpt2(n:i32) -> Vec<f64> {
  let mut rng = thread_rng();
  let gam = Gamma::new(1.088870,1.2336997422).unwrap();
  let u = Uniform::new(0.0,1.0).unwrap();
  let n: usize = n as usize;
  let mut passage_times: Vec<f64> = Vec::new();
  let mut it_n: f64=0.0;
  let mut s:f64;
  for _id in 0..n {
    loop {
      let x:f64 = gam.sample(&mut rng);
      let w:f64 = u.sample(&mut rng) * 1.243707 * gam.pdf(x);
      s = g(0.0,x);
      let mut floor_n = (0.275*x).ceil().max(3.0) as i32;
      if floor_n%2==0 {
        floor_n+=1;
      }
      while (it_n as i32)<floor_n {
        s = s + g(-it_n,x)+g(it_n,x);
        it_n+=1.0;
      }
      loop {
        it_n+=2.0;
        let s_step = g(-it_n,x)+g(it_n,x)+g(-it_n+1.0,x)+g(it_n-1.0,x);
        s+=s_step;
        if (w-s).abs()>s_step {
          break;
        }
      }
      if w<=s {
        passage_times.push(x);
        break;
      }
    }
  }
  passage_times
}

#[extendr]
fn algo5(x_start:f64,x_end:f64,t_start:f64,t_end:f64,t_sim:f64) -> f64 {
  let mut rng = thread_rng();
  let u = Uniform::new(0.0,1.0).unwrap();
  let std_n = Normal::new(0.0, ((t_end-t_sim)*(t_sim-t_start)).powf(0.5)/(t_end-t_start)).unwrap();
  let theta = (x_start-x_end).abs();
  let m:f64 = if x_end<x_start {-1.0} else {1.0};
  let n0 = (((t_end-t_sim)+4.0*theta.powf(2.0)).powf(0.5)/(4.0*theta).ceil()) as usize;
  let mut it_n: usize;
  let mut pnu: f64;
  let mut pnd: f64;
  loop {
    let w = x_start + m*((t_end-t_start)*(((theta*(t_end-t_sim))/(t_end-t_start).powf(1.5)+std_n.sample(&mut rng)).powf(2.0)+std_n.sample(&mut rng).powf(2.0)+std_n.sample(&mut rng).powf(2.0))).powf(0.5);
    let u_test: f64 = u.sample(&mut rng);
    if a6(x_start,w,t_start,t_sim,theta,x_start,u_test) && a6(w,x_end,t_sim,t_end,theta,x_start,u_test) {
      return w;
    }
  }
}

/*fn b1(s:f64,t:f64,l:f64,v:f64,j:usize,x:f64,y:f64) -> f64 {
  (-2.0*((v-l)*(j as f64)+l-x)*((v-l)*(j as f64)+l-y)/(t-s)).exp()
}*/

fn b1(s:f64,t:f64,theta:f64,origin:f64,j:usize,x:f64,y:f64) -> f64 {
  (-2.0*(2.0*(j as f64)*theta - (theta+x-origin))*(2.0*theta*(j as f64) - (theta + y - origin))/(t-s)).exp()
}

/*fn b2(s:f64,t:f64,l:f64,v:f64,j:usize,x:f64,y:f64) -> f64 {
  (-2.0*(j as f64)*((v-l).powf(2.0)*(j as f64)+(v-l)*(x-y))/(t-s)).exp()
}*/

fn b2(s:f64,t:f64,theta:f64,origin:f64,j:usize,x:f64,y:f64) -> f64 {
  (-2.0*(j as f64)*(4.0*theta.powf(2.0)*(j as f64)+2.0*theta*(x-y))/s).exp()
}

#[extendr]
fn a6(x:f64,y:f64,s:f64,t:f64,theta:f64,x_start:f64,u:f64) -> bool {
  let mut k:usize = 1;
  let mut s_2k1: f64 = 0.0;
  let mut s_2k:f64 = 1.0;
  while (s_2k-u)*(u-s_2k1)>0.0 {
    s_2k = 1.0 - (1..=k).map(|j| b1(s,t,theta,x_start,j,x,y)+b1(s,t,theta,-x_start,j,-x,-y)-b2(s,t,theta,x_start,j,x,y)-b2(s,t,theta,-x_start,j,-x,-y)).sum::<f64>();
    s_2k1 = s_2k-b1(s,t,theta,x_start,k+1,x,y)-b1(s,t,theta,-x_start,k+1,-x,-y);
    k+=1;
  }
  if u<s_2k1 {
    return true;
  } else {
    return false;
  }
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod cts_smc;
    fn cts_smc_rust;
    fn path_smc_rust;
    fn euler_approx_rust;
    fn brownian_motion_fpt;
    fn brownian_motion_fpt2;
    fn brownian_bridge;
    fn algo5;
    fn a6;
}

