use num::Signed;

use crate::{consts::MYR2SEC, input::IcyDwarfInput, planet_system::WorldState};

const ORB_D_TIME: f64 = 5.0e-4 * 1.0e-6 * MYR2SEC;

impl IcyDwarfInput {
    pub fn orbit(
        &self,
        world_idx: usize,
        world_states: &mut [WorldState],
        d_time: f64,
        real_time: f64,
        q_prim: f64,
    ) {
        let current_world_state = &mut world_states[world_idx];
        let mut res_orbit_evol = current_world_state.res_acct_for[world_idx..]
            .iter()
            .any(|&n| n > 0.);

        let nv = 6;
        let n_param_orbit = 30;

        for i in 0..world_idx {
            if current_world_state.res_acct_for[i] > 0. {
                res_orbit_evol = true;
                let j = current_world_state.res_acct_for[i] + 1.;
            }
        }
    }
    pub fn res_check(world_state: &mut WorldState) {}
}

pub fn laplace_coef(alpha: f64, j: f64, s: f64) -> f64 {
    let mut b_lap_j = 1_f64;
    let mut temp = 1_f64;
    for m in 0..200 {
        let m = m as f64;
        temp = temp * (s + m) / (1. + m) * (s + j + m) / (j + 1. + m) * alpha.powi(2);
        b_lap_j += temp * (j + m + 1.) * 2.;
        if temp < 1.0e-6 {
            break;
        }
    }
    b_lap_j
        * alpha.powf(j)
        * (0..200)
            .map(|m| (s + m as f64) / (1. + m as f64))
            .product::<f64>()
        * 2.
}

pub fn d_laplace_coef(alpha: f64, j: f64, s: f64) -> f64 {
    let mut db_lap_j = 1_f64;
    let mut temp = 1_f64;
    for m in 0..200 {
        let m = m as f64;
        temp = temp * (s + m) / (1. + m) * (s + j + m) / (j + 1. + m) * alpha.powi(2);
        db_lap_j += temp;
        if temp * (j + m + 1.) * 2. < 1.0e-6 {
            break;
        }
    }
    db_lap_j
        * alpha.powf(j - 1.)
        * (0..200)
            .map(|m| (s + m as f64) / (1. + m as f64))
            .product::<f64>()
        * 2.
}

pub fn d2_laplace_coef(alpha: f64, j: f64, s: f64) -> f64 {
    let mut d2b_lapj = j * (j - 1.);
    let mut temp = 0.;

    for m in 0..200 {
        let m = m as f64;
        temp *= (s + m) / (1. + m) * (s + j + m) / (j + 1. + m) * alpha.powi(2);
        let incr = temp + (j + 2. * (m + 1.)) * (j + 2. * (m + 1.) - 1.);
        d2b_lapj += incr;
        if incr < 1.0e-6 {
            break;
        }
    }
    d2b_lapj *= alpha.powf(j - 2.);
    for m in 0..j as i32 {
        let m = m as f64;
        d2b_lapj *= (s + m) / (m + 1.);
    }
    d2b_lapj * 2.
}

pub fn mmr_avg_ham() -> [f64; 5] {
    todo!()
}

pub struct OrbitOut {}
