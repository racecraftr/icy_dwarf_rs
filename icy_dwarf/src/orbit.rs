use crate::{
    consts::{GCGS, IJMAX, MYR2SEC},
    input::IcyDwarfInput,
    planet_system::WorldState,
};

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

    /// Checks for orbital mean-motion resonances between two moons by
    /// comparing the mean motion of a given moon `im` with that of all other moons.
    pub fn res_check(
        &self,
        im: usize,
        world_states: &[WorldState],
        resonance: &mut [Vec<f64>],
        p_capture: &mut [Vec<f64>],
        real_time: f64,
        a_old: &[f64],
        t_tide: &[f64],
        q_prim: f64,
    ) {
        let nmoons = world_states.len();
        for i in 0..nmoons {
            resonance[im][i] = 0.0;
        }

        let m_prim = self.primary_world.mass;
        let r_prim = self.primary_world.rad;
        let k2_prim = self.primary_world.k2;
        let prim_sign = if self.worlds[im].retrograde { -1.0 } else { 1.0 };
        let reslock = real_time > self.worlds[im].t_reslock * MYR2SEC;

        for i in 0..nmoons {
            let tzero_i = self.worlds[i].t_form * MYR2SEC;
            if i < im && real_time >= tzero_i {
                // Find index of inner and outer moon
                let (inner, outer) = if world_states[im].a_orb < world_states[i].a_orb {
                    (im, i)
                } else {
                    (i, im)
                };

                for j_idx in (1..=IJMAX).rev() {
                    let j = j_idx as f64;
                    // k=1 (l=1) handles k=1 resonances only
                    for l_idx in 1..=1 {
                        let l = l_idx as f64;

                        let is_already_captured = world_states[inner].res_acct_for[outer] > 0.0;
                        let commensurability = if is_already_captured {
                            (a_old[inner] / a_old[outer]).powf(-1.5) * j / (j + l)
                        } else {
                            (world_states[inner].a_orb / world_states[outer].a_orb).powf(-1.5) * j / (j + l)
                        };

                        let tol_lower = if is_already_captured { 0.985 } else { 0.99 };
                        let tol_upper = if is_already_captured { 1.015 } else { 1.01 };

                        if commensurability > tol_lower && commensurability < tol_upper {
                            resonance[inner][outer] = j;
                            resonance[outer][inner] = j;

                            let dnorb_dt_inner;
                            let dnorb_dt_outer;

                            if reslock {
                                dnorb_dt_inner = -1.5 * (GCGS * m_prim).sqrt()
                                    * world_states[inner].a_orb.powf(-2.5)
                                    * prim_sign * world_states[inner].a_orb / t_tide[inner];
                                dnorb_dt_outer = -1.5 * (GCGS * m_prim).sqrt()
                                    * world_states[outer].a_orb.powf(-2.5)
                                    * prim_sign * world_states[outer].a_orb / t_tide[outer];
                            } else {
                                dnorb_dt_inner = -1.5 * (GCGS * m_prim).sqrt()
                                    * world_states[inner].a_orb.powf(-2.5)
                                    * prim_sign * 3.0 * k2_prim * (GCGS / m_prim).sqrt()
                                    * r_prim.powi(5) * self.worlds[inner].mass() / q_prim
                                    * world_states[inner].a_orb.powf(-5.5);
                                dnorb_dt_outer = -1.5 * (GCGS * m_prim).sqrt()
                                    * world_states[outer].a_orb.powf(-2.5)
                                    * prim_sign * 3.0 * k2_prim * (GCGS / m_prim).sqrt()
                                    * r_prim.powi(5) * self.worlds[outer].mass() / q_prim
                                    * world_states[outer].a_orb.powf(-5.5);
                            }

                            if j * dnorb_dt_inner <= (j + l) * dnorb_dt_outer {
                                let p_cap = self.mmr_capture_probability(
                                    world_states,
                                    inner,
                                    outer,
                                    j,
                                );
                                p_capture[inner][outer] = p_cap;
                                p_capture[outer][inner] = p_cap;
                            } else {
                                p_capture[inner][outer] = 0.0;
                                p_capture[outer][inner] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Screens out multiple resonances to keep only the strongest or already in-place ones.
    pub fn resscreen(
        &self,
        resonance: &[f64],
        res_acct_for: &mut [f64],
        res_acct_for_old: &[f64],
    ) {
        let nmoons = resonance.len();
        res_acct_for.fill(0.0);

        let mut res_min = IJMAX as f64;
        let mut nbres = 0;

        // Find the min order of resonance for each moon and the number of moons
        // involved in resonances of this order
        for i in 0..nmoons {
            if resonance[i] > 0.0 && resonance[i] <= res_min {
                res_min = resonance[i];
            }
        }
        for i in 0..nmoons {
            if resonance[i] == res_min {
                nbres += 1;
            }
        }

        // Copy only the lowest-order resonances for each moon
        for i in 0..nmoons {
            if resonance[i] == res_min {
                res_acct_for[i] = resonance[i];
            }
        }

        // Zero out newer resonances if there are multiple of lowest order
        if nbres > 1 {
            for i in 0..nmoons {
                if res_acct_for[i] > 0.0 && res_acct_for_old[i] == 0.0 {
                    res_acct_for[i] = 0.0;
                }
            }
        }
    }

    /// Calculates capture probability in a resonance.
    pub fn mmr_capture_probability(
        &self,
        world_states: &[WorldState],
        inner: usize,
        outer: usize,
        j: f64,
    ) -> f64 {
        let alpha = world_states[inner].a_orb / world_states[outer].a_orb;
        let b_lap_j = laplace_coef(alpha, j, 0.5);
        let db_lap_j = d_laplace_coef(alpha, j, 0.5);

        let mut ck = (2. * j + 1.) * b_lap_j + alpha * db_lap_j;
        if j == 1. {
            ck -= 1. / (alpha * alpha);
        }

        let dk = (3. * (j + 1.).powi(2)
            / (2_f64.sqrt() * self.worlds[inner].mass() / self.primary_world.mass * ck))
            .powf(2. / 3.);
        let r = dk * world_states[inner].e_orb.powi(2);
        
        let pk = if r <= 3. {
            1.
        } else {
            1. / (r - 2.).powf(2.4)
                - 0.43 * (r / 3.).ln() / 10_f64.ln()
                - 0.37 * ((-r + 3.).exp() - 1.)
        };
        pk.max(0.)
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
    d2b_lapj
        * alpha.powf(j - 2.)
        * (0..j as i32)
            .map(|m| {
                let m = m as f64;
                (s + m) / (1. + m)
            })
            .product::<f64>()
        * 2.
}

pub fn mmid(y: &[f64], dydx: &[f64], derivs: fn(f64, &[f64], &[f64]) -> Vec<f64>) {}

pub struct MmrAvgHamParams {}
