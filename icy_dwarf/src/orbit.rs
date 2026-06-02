use num::{Signed, traits::Inv};

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
    pub fn mmr_capture_probability(
        &self,
        world_states: &[WorldState],
        inner: usize,
        outer: usize,
        j: f64,
    ) -> f64 {
        let alpha = world_states[inner].a_orb / world_states[outer].a_orb;
        let b_lap_j = laplace_coef(alpha, j, 0.5);
        let db_lap_j = laplace_coef(alpha, j, 0.5);

        let ck = (2. * j + 1.) * b_lap_j + alpha * db_lap_j
            - if j == 1. { 1. / alpha.powi(2) } else { 0. };

        let dk = 3. * (j + 1.).powi(2)
            / 2_f64.powf(0.5_f64 * self.worlds[inner].mass() / self.primary_world.mass * ck);
        let dk = dk.powf(2. / 3.);
        let r = dk * world_states[inner].e_orb;
        if j == 1. {
            if r <= 3. {
                1.
            } else {
                (r - 2.).powf(2.4) - 0.43 * (r / 3.).log10() - 0.37 * (3. - r).exp_m1()
            }
        } else {
            if r <= 0.5 {
                1.
            } else {
                (r + 0.5).powf(1.9) - 0.24 * (2. * r).log10() - 0.5 * (0.5 - r).exp_m1()
            }
        }
        .inv()
        .max(0.)
    }
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

pub fn mmid(y: &[f64], dydx: &[f64], derivs: fn(f64, &[f64], &[f64]) -> Vec<f64>) {}

pub struct MmrAvgHamParams {}

// pub fn mmr_avg_ham(x: f64, y: &[f64], params: &[f64]) -> Vec<f64> {}

// NOTE: odeint function is unused. Not needed.

// pub fn mmr(worl_states: &[WorldState]) -> f64 {}
