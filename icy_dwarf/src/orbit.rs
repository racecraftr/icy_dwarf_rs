use crate::{
    consts::{GCGS, GYR2SEC, IJMAX, MIN_ECC, MYR2SEC},
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
        let nmoons = world_states.len();
        let mut resorbevol = world_states[world_idx].res_acct_for[world_idx..nmoons]
            .iter()
            .any(|&r| r > 0.0);

        for i in 0..world_idx {
            if world_states[world_idx].res_acct_for[i] > 0.0 {
                resorbevol = true;
                let j = world_states[world_idx].res_acct_for[i] + 1.0;

                let has_old = world_states[world_idx].a_old != 0.0 && world_states[i].a_old != 0.0;

                let (idx_0, idx_1) = if world_states[world_idx].a_orb < world_states[i].a_orb {
                    (world_idx, i)
                } else {
                    (i, world_idx)
                };

                let m_prim = self.primary_world.mass;
                let r_prim = self.primary_world.rad;
                let k2_prim = self.primary_world.k2;

                let m_0 = self.worlds[idx_0].mass();
                let r_0 = self.worlds[idx_0].planetary_rad;
                let w_0 = world_states[idx_0].w_tide_tot + world_states[idx_0].w_fluidtide_tot;

                let m_1 = self.worlds[idx_1].mass();
                let r_1 = self.worlds[idx_1].planetary_rad;
                let w_1 = world_states[idx_1].w_tide_tot + world_states[idx_1].w_fluidtide_tot;

                let psgn_0 = if self.worlds[idx_0].retrograde {
                    -1.0
                } else {
                    1.0
                };
                let psgn_1 = if self.worlds[idx_1].retrograde {
                    -1.0
                } else {
                    1.0
                };

                let ttide_0 =
                    self.worlds[idx_0].t_reslock * MYR2SEC * real_time / (4568.2 * MYR2SEC);
                let ttide_1 =
                    self.worlds[idx_1].t_reslock * MYR2SEC * real_time / (4568.2 * MYR2SEC);

                let (cs_ee, cs_eep, cr_e, cr_ep, cr_ee, cr_eep, cr_epep) = if has_old {
                    (
                        world_states[world_idx].cs_ee_old,
                        world_states[world_idx].cs_eep_old,
                        world_states[world_idx].cr_e_old,
                        world_states[world_idx].cr_ep_old,
                        world_states[world_idx].cr_ee_old,
                        world_states[world_idx].cr_eep_old,
                        world_states[world_idx].cr_epep_old,
                    )
                } else {
                    let alpha = ((j - 1.0) / j).powf(2.0 / 3.0);
                    let cs_ee = 0.125
                        * (2.0 * alpha * d_laplace_coef(alpha, 0.0, 0.5)
                            + alpha * alpha * d2_laplace_coef(alpha, 0.0, 0.5));
                    let cs_eep = 0.25
                        * (2.0 * laplace_coef(alpha, 1.0, 0.5)
                            - 2.0 * alpha * d_laplace_coef(alpha, 1.0, 0.5)
                            - alpha * alpha * d2_laplace_coef(alpha, 1.0, 0.5));
                    let cr_e = 0.5
                        * (-2.0 * j * laplace_coef(alpha, j, 0.5)
                            - alpha * d_laplace_coef(alpha, j, 0.5));
                    let cr_ep = 0.5
                        * ((2.0 * j - 1.0) * laplace_coef(alpha, j - 1.0, 0.5)
                            + alpha * d_laplace_coef(alpha, j - 1.0, 0.5))
                        - if j == 2.0 { 2.0 * alpha } else { 0.0 };
                    let p = 2.0 * j;
                    let cr_ee = 0.125
                        * ((-5.0 * p + 4.0 * p * p) * laplace_coef(alpha, p, 0.5)
                            + (-2.0 + 4.0 * p) * alpha * d_laplace_coef(alpha, p, 0.5)
                            + alpha * alpha * d2_laplace_coef(alpha, p, 0.5));
                    let cr_eep = 0.25
                        * ((-2.0 + 6.0 * p - 4.0 * p * p) * laplace_coef(alpha, p - 1.0, 0.5)
                            + (2.0 - 4.0 * p) * alpha * d_laplace_coef(alpha, p - 1.0, 0.5)
                            - alpha * alpha * d2_laplace_coef(alpha, p - 1.0, 0.5));
                    let cr_epep = 0.125
                        * ((2.0 - 7.0 * p + 4.0 * p * p) * laplace_coef(alpha, p - 2.0, 0.5)
                            + (-2.0 + 4.0 * p) * alpha * d_laplace_coef(alpha, p - 2.0, 0.5)
                            + alpha * alpha * d2_laplace_coef(alpha, p - 2.0, 0.5));
                    (cs_ee, cs_eep, cr_e, cr_ep, cr_ee, cr_eep, cr_epep)
                };

                let (h_0, k_0, a0, h_1, k_1, a1, w_0, w_1) = if has_old {
                    let h_0 = world_states[idx_0].h_old;
                    let k_0 = world_states[idx_0].k_old;
                    let a0 = world_states[idx_0].a_old;
                    let h_1 = world_states[idx_1].h_old;
                    let k_1 = world_states[idx_1].k_old;
                    let a1 = world_states[idx_1].a_old;

                    let factor_0 = 11.5
                        * r_0.powi(5)
                        * (GCGS * m_prim / a0.powi(3)).powf(2.5)
                        * (h_0 * h_0 + k_0 * k_0)
                        / GCGS;
                    let w_0 = if factor_0 > 0.0 { w_0 / factor_0 } else { 0.0 };

                    let factor_1 = 11.5
                        * r_1.powi(5)
                        * (GCGS * m_prim / a1.powi(3)).powf(2.5)
                        * (h_1 * h_1 + k_1 * k_1)
                        / GCGS;
                    let w_1 = if factor_1 > 0.0 { w_1 / factor_1 } else { 0.0 };

                    (h_0, k_0, a0, h_1, k_1, a1, w_0, w_1)
                } else {
                    let a_0 = world_states[idx_0].a_orb;
                    let e_0 = world_states[idx_0].e_orb;
                    let lamb_0 = world_states[idx_0].lambda;
                    let omeg_0 = world_states[idx_0].omega;

                    let a_1 = world_states[idx_1].a_orb;
                    let e_1 = world_states[idx_1].e_orb;
                    let lamb_1 = world_states[idx_1].lambda;
                    let omeg_1 = world_states[idx_1].omega;

                    let factor_0 =
                        11.5 * r_0.powi(5) * (GCGS * m_prim / a_0.powi(3)).powf(2.5) * e_0.powi(2)
                            / GCGS;
                    let w_0 = if factor_0 > 0.0 { w_0 / factor_0 } else { 0.0 };

                    let factor_1 =
                        11.5 * r_1.powi(5) * (GCGS * m_prim / a_1.powi(3)).powf(2.5) * e_1.powi(2)
                            / GCGS;
                    let w_1 = if factor_1 > 0.0 { w_1 / factor_1 } else { 0.0 };

                    let sigma_0 = j * lamb_1 + (1.0 - j) * lamb_0 - omeg_0;
                    let l_0 = (m_0 * GCGS * m_0 * m_prim * a_0).sqrt();
                    let sigma_val_0 = l_0 * (1.0 - (1.0 - e_0 * e_0).max(0.0).sqrt());

                    let sigma_1 = j * lamb_1 + (1.0 - j) * lamb_0 - omeg_1;
                    let l_1 = (m_1 * GCGS * m_1 * m_prim * a_1).sqrt();
                    let sigma_val_1 = l_1 * (1.0 - (1.0 - e_1 * e_1).max(0.0).sqrt());

                    let lambda_0 = l_0 - (1.0 - j) * (sigma_val_0 + sigma_val_1);
                    let lambda_1 = l_1 - j * (sigma_val_0 + sigma_val_1);

                    let h_0 = e_0 * sigma_0.cos();
                    let k_0 = e_0 * sigma_0.sin();
                    let a0 = lambda_0 * lambda_0 / (GCGS * m_0 * m_0 * m_prim);

                    let h_1 = e_1 * sigma_1.cos();
                    let k_1 = e_1 * sigma_1.sin();
                    let a1 = lambda_1 * lambda_1 / (GCGS * m_1 * m_1 * m_prim);

                    (h_0, k_0, a0, h_1, k_1, a1, w_0, w_1)
                };

                let ring_mass = self.primary_world.ring.mass;
                let aring_in = self.primary_world.ring.inner;
                let aring_out = self.primary_world.ring.outer;
                let ring_surface_density = if ring_mass > 0.0 && aring_out > aring_in {
                    ring_mass
                        / (std::f64::consts::PI * (aring_out * aring_out - aring_in * aring_in))
                } else {
                    0.0
                };

                let params = MmrAvgHamParams {
                    m: [m_0, m_1],
                    r: [r_0, r_1],
                    w: [w_0, w_1],
                    j,
                    m_prim,
                    r_prim,
                    j2_prim: self.primary_world.j2,
                    j4_prim: self.primary_world.j4,
                    k2_prim,
                    q_prim,
                    cs_ee,
                    cs_eep,
                    cr_e,
                    cr_ep,
                    cr_ee,
                    cr_eep,
                    cr_epep,
                    speedup: self.grid.speedup,
                    psgn: [psgn_0, psgn_1],
                    reslock: self.primary_world.tidal_resonant,
                    ttide: [ttide_0, ttide_1],
                    alpha_lind: if ring_surface_density <= 2.0 {
                        2.0e-5
                    } else {
                        1.0e-4
                    },
                    ring_surface_density,
                    aring_in,
                    aring_out,
                };

                let mut y = vec![h_0, k_0, a0, h_1, k_1, a1];
                let dydx = vec![0.0; 6];

                let steps = (d_time / self.grid.speedup / ORB_D_TIME) as i64;
                for _q in 0..steps {
                    if y[0] * y[0] + y[1] * y[1] < 1.0e-14 {
                        if y[0].abs() > y[1].abs() {
                            y[0] = 1.3e-4;
                        } else {
                            y[1] = 1.3e-4;
                        }
                    }
                    if y[3] * y[3] + y[4] * y[4] < 1.0e-14 {
                        if y[3].abs() > y[4].abs() {
                            y[3] = 1.3e-4;
                        } else {
                            y[4] = 1.3e-4;
                        }
                    }
                    y = mmid(&y, &dydx, &params, 0.0, ORB_D_TIME, 10, mmr_avg_ham);
                }

                let h_0 = y[0];
                let k_0 = y[1];
                let a0 = y[2];
                let h_1 = y[3];
                let k_1 = y[4];
                let a1 = y[5];

                let e0 = (h_0 * h_0 + k_0 * k_0).sqrt();
                let e1 = (h_1 * h_1 + k_1 * k_1).sqrt();

                let e0_clamped = e0.min(0.99999);
                let e1_clamped = e1.min(0.99999);

                let lambda_0 = (a0 * GCGS * m_0 * m_0 * m_prim).sqrt();
                let lambda_1 = (a1 * GCGS * m_1 * m_1 * m_prim).sqrt();

                let term_e0 = (1.0 - e0_clamped * e0_clamped).sqrt();
                let term_e1 = (1.0 - e1_clamped * e1_clamped).sqrt();

                let num = lambda_1 + j * (lambda_0 + lambda_1) * (1.0 / term_e0 - 1.0);
                let den = 1.0 - j * (1.0 - term_e1 / term_e0);
                let l_1 = if den.abs() > 1.0e-30 { num / den } else { num };
                let l_0 = (lambda_0 + lambda_1) / term_e0 - l_1 * term_e1 / term_e0;

                let a_0 = l_0 * l_0 / (GCGS * m_0 * m_0 * m_prim);
                let a_1 = l_1 * l_1 / (GCGS * m_1 * m_1 * m_prim);

                let (
                    final_a_im,
                    final_e_im,
                    final_h_im,
                    final_k_im,
                    final_a__im,
                    final_a_i,
                    final_e_i,
                    final_h_i,
                    final_k_i,
                    final_a__i,
                ) = if world_idx == idx_0 {
                    (a_0, e0, h_0, k_0, a0, a_1, e1, h_1, k_1, a1)
                } else {
                    (a_1, e1, h_1, k_1, a1, a_0, e0, h_0, k_0, a0)
                };

                let (left, right) = world_states.split_at_mut(world_idx);
                let w_im = &mut right[0];
                let w_i = &mut left[i];

                w_im.a_orb = final_a_im;
                w_im.e_orb = final_e_im;
                w_i.a_orb = final_a_i;
                w_i.e_orb = final_e_i;

                w_im.h_old = final_h_im;
                w_im.k_old = final_k_im;
                w_im.a_old = final_a__im;
                w_i.h_old = final_h_i;
                w_i.k_old = final_k_i;
                w_i.a_old = final_a__i;

                w_im.cs_ee_old = cs_ee;
                w_i.cs_ee_old = cs_ee;
                w_im.cs_eep_old = cs_eep;
                w_i.cs_eep_old = cs_eep;
                w_im.cr_e_old = cr_e;
                w_i.cr_e_old = cr_e;
                w_im.cr_ep_old = cr_ep;
                w_i.cr_ep_old = cr_ep;
                w_im.cr_ee_old = cr_ee;
                w_i.cr_ee_old = cr_ee;
                w_im.cr_eep_old = cr_eep;
                w_i.cr_eep_old = cr_eep;
                w_im.cr_epep_old = cr_epep;
                w_i.cr_epep_old = cr_epep;

                if [e0, e1].iter().any(|&e| (0..=1).contains(&(e as i32))) {
                    use std::fs::OpenOptions;
                    use std::io::Write;
                    if let Ok(mut file) = OpenOptions::new()
                        .append(true)
                        .open("Outputs/Resonances.txt")
                    {
                        let _ = writeln!(
                            file,
                            "Time {} Myr, eccentricity out of bounds. Stopping.",
                            real_time / MYR2SEC
                        );
                    }
                    eprintln!(
                        "Time {} Myr, eccentricity out of bounds. Stopping. e1={}, e2={}",
                        real_time / MYR2SEC,
                        e0,
                        e1
                    );
                    std::process::exit(1);
                }

                w_im.n_orb = if w_im.a_orb > 0.0 {
                    (GCGS * m_prim / w_im.a_orb.powi(3)).sqrt()
                } else {
                    0.0
                };
            }
        }

        if !resorbevol {
            world_states[world_idx].h_old = 0.0;
            world_states[world_idx].k_old = 0.0;
            world_states[world_idx].a_old = 0.0;
            world_states[world_idx].cs_ee_old = 0.0;
            world_states[world_idx].cs_eep_old = 0.0;
            world_states[world_idx].cr_e_old = 0.0;
            world_states[world_idx].cr_ep_old = 0.0;
            world_states[world_idx].cr_ee_old = 0.0;
            world_states[world_idx].cr_eep_old = 0.0;
            world_states[world_idx].cr_epep_old = 0.0;

            let ctl = matches!(self.world_spec.ecc_model, crate::input::EccModel::E10Ctl);

            let mass_moon = self.worlds[world_idx].mass();
            let aorb_im = world_states[world_idx].a_orb;
            let eorb_im = world_states[world_idx].e_orb;
            let w_tot =
                world_states[world_idx].w_tide_tot + world_states[world_idx].w_fluidtide_tot;

            let m_prim = self.primary_world.mass;
            let r_prim = self.primary_world.rad;
            let k2_prim = self.primary_world.k2;

            let (d_eorb_moon, d_aorb_moon, d_eorb_pl, d_aorb_pl) = if !ctl {
                let d_e_moon = if eorb_im > 0.0 {
                    -w_tot * aorb_im / (GCGS * m_prim * mass_moon * eorb_im)
                } else {
                    0.0
                };
                let d_a_moon = -2.0 * w_tot * aorb_im * aorb_im / (GCGS * m_prim * mass_moon);

                let prim_sign_im = if self.worlds[world_idx].retrograde {
                    -1.0
                } else {
                    1.0
                };
                let d_e_pl =
                    57.0 / 8.0 * k2_prim * (GCGS / m_prim).sqrt() * r_prim.powi(5) * mass_moon
                        / q_prim
                        * aorb_im.powf(-6.5)
                        * eorb_im;

                let reslock = self.primary_world.tidal_resonant;
                let t_tide_im =
                    self.worlds[world_idx].t_reslock * MYR2SEC * real_time / (4568.2 * MYR2SEC);

                let d_a_pl = if !reslock {
                    prim_sign_im
                        * 3.0
                        * k2_prim
                        * (GCGS / m_prim).sqrt()
                        * r_prim.powi(5)
                        * mass_moon
                        / q_prim
                        * aorb_im.powf(-5.5)
                } else if t_tide_im > 0.0 {
                    prim_sign_im * aorb_im / t_tide_im
                } else {
                    0.0
                };

                (d_e_moon, d_a_moon, d_e_pl, d_a_pl)
            } else {
                let rad_moon = self.worlds[world_idx].planetary_rad;
                let i_m = world_states[world_idx].moi * mass_moon * rad_moon * rad_moon;
                let i_p = self.primary_world.moi_coef * m_prim * r_prim * r_prim;

                let norb_im = world_states[world_idx].n_orb;
                let q_tide_im = world_states[world_idx].q_tide;

                let e2 = eorb_im * eorb_im;
                let term_e = (1.0 - e2).max(0.0).sqrt();

                let tau_p = if norb_im > 0.0 {
                    1.0 / (2.0 * norb_im * q_prim)
                } else {
                    0.0
                };
                let tau_m = if norb_im > 0.0 && q_tide_im > 0.0 {
                    1.0 / (2.0 * norb_im * q_tide_im)
                } else {
                    0.0
                };

                let k_p = 3.0
                    * k2_prim
                    * tau_p
                    * (GCGS * mass_moon * mass_moon / r_prim)
                    * (r_prim / aorb_im).powi(6)
                    * norb_im.powi(2);
                let k_m = 3.0
                    * world_states[world_idx].k2
                    * tau_m
                    * (GCGS * m_prim.powi(2) / rad_moon)
                    * (rad_moon / aorb_im).powi(6)
                    * norb_im.powi(2);

                let x_p = world_states[world_idx].i_orb.cos();
                let x_m = world_states[world_idx].obl.cos();

                let n_prim = (2.0 * std::f64::consts::PI
                    / (self.primary_world.spin_period * 3600.0))
                    .max(0.);

                let eta_p = if aorb_im > 0.0 && norb_im > 0.0 && term_e > 0.0 {
                    (m_prim + mass_moon) / (m_prim * mass_moon) * i_p * n_prim
                        / (aorb_im * aorb_im * norb_im * term_e)
                } else {
                    0.0
                };

                let eta_m = if aorb_im > 0.0 && norb_im > 0.0 && term_e > 0.0 {
                    (m_prim + mass_moon) / (m_prim * mass_moon) * i_m * world_states[world_idx].spin
                        / (aorb_im * aorb_im * norb_im * term_e)
                } else {
                    0.0
                };

                let e4 = e2 * e2;
                let e6 = e2 * e4;
                let e8 = e4 * e4;

                let n_val =
                    (1.0 + 7.5 * e2 + 45.0 / 8.0 * e4 + 5.0 / 16.0 * e6) / (1.0 - e2).powi(6);
                let n_a =
                    (1.0 + 16.5 * e2 + 255.0 / 8.0 * e4 + 185.0 / 16.0 * e6 + 25.0 / 64.0 * e8)
                        / (1.0 - e2).powf(7.5);
                let o_e = (1.0 + 1.5 * e2 + 1.0 / 8.0 * e4) / (1.0 - e2).powi(5);
                let n_e =
                    (1.0 + 3.25 * e2 + 15.0 / 8.0 * e4 + 5.0 / 64.0 * e6) / (1.0 - e2).powf(6.5);
                let o_e_val = (1.0 + 3.0 * e2 + 3.0 / 8.0 * e4) / (1.0 - e2).powf(4.5);

                let d_a_pl = if norb_im > 0.0 {
                    4.0 * aorb_im.powi(2) / (GCGS * m_prim * mass_moon)
                        * k_p
                        * (n_val * x_p * n_prim / norb_im - n_a)
                } else {
                    0.0
                };
                let d_a_moon = if norb_im > 0.0 {
                    4.0 * aorb_im.powi(2) / (GCGS * m_prim * mass_moon)
                        * k_m
                        * (n_val * x_m * world_states[world_idx].spin / norb_im - n_a)
                } else {
                    0.0
                };

                let d_e_pl = if norb_im > 0.0 {
                    11.0 * aorb_im * eorb_im / (GCGS * m_prim * mass_moon)
                        * k_p
                        * (o_e * x_p * n_prim / norb_im - 18.0 / 11.0 * n_e)
                } else {
                    0.0
                };
                let d_e_moon = if norb_im > 0.0 {
                    11.0 * aorb_im * eorb_im / (GCGS * m_prim * mass_moon)
                        * k_m
                        * (o_e * x_m * world_states[world_idx].spin / norb_im - 18.0 / 11.0 * n_e)
                } else {
                    0.0
                };

                let d_spin_moon = if norb_im > 0.0 {
                    -k_m / (i_m * norb_im)
                        * ((1.0 + x_m * x_m) * o_e_val * world_states[world_idx].spin / norb_im
                            - 2.0 * x_m * n_val)
                } else {
                    0.0
                };

                let d_iorb = if n_prim > 0.0 && norb_im > 0.0 {
                    world_states[world_idx].i_orb.sin() * k_p / (i_p * n_prim * norb_im)
                        * ((x_p - eta_p) * o_e_val * n_prim / norb_im - 2.0 * n_val)
                } else {
                    0.0
                };
                let d_obl = if world_states[world_idx].spin > 0.0 && norb_im > 0.0 {
                    world_states[world_idx].obl.sin() * k_m
                        / (i_m * world_states[world_idx].spin * norb_im)
                        * ((x_m - eta_m) * o_e_val * world_states[world_idx].spin / norb_im
                            - 2.0 * n_val)
                } else {
                    0.0
                };

                world_states[world_idx].spin += d_time * d_spin_moon;
                world_states[world_idx].i_orb += d_time * d_iorb;
                world_states[world_idx].obl += d_time * d_obl;

                (d_e_moon, d_a_moon, d_e_pl, d_a_pl)
            };

            let ring_mass = self.primary_world.ring.mass;
            let aring_in = self.primary_world.ring.inner;
            let aring_out = self.primary_world.ring.outer;
            let ring_surface_density = if ring_mass > 0.0 && aring_out > aring_in {
                ring_mass / (std::f64::consts::PI * (aring_out * aring_out - aring_in * aring_in))
            } else {
                0.0
            };

            let (d_aorb_ring, d_eorb_ring) = if ring_surface_density > 0.0 {
                let mut ring_torque = 0.0;
                let den_in = 1.0 - (aring_in / aorb_im).powf(1.5);
                let den_out = 1.0 - (aring_out / aorb_im).powf(1.5);
                if den_in.abs() > 1.0e-30 && den_out.abs() > 1.0e-30 {
                    let kmin = (1.0 / den_in).floor() as i32 + 1;
                    let kmax = (1.0 / den_out).floor() as i32;
                    let alpha_lind: f64 = if ring_surface_density <= 2.0 {
                        2.0e-5
                    } else {
                        1.0e-4
                    };
                    let kmax_limit = (1.0 / alpha_lind.sqrt()).floor() as i32;
                    if kmin <= kmax && kmax <= kmax_limit {
                        for k_idx in kmin..=kmax {
                            ring_torque += std::f64::consts::PI * std::f64::consts::PI / 3.0
                                * ring_surface_density
                                * GCGS
                                * mass_moon
                                * mass_moon
                                * (k_idx as f64)
                                * ((k_idx - 1) as f64)
                                * aorb_im
                                / m_prim;
                        }
                    }
                }
                (
                    2.0 * ring_torque / mass_moon * (aorb_im / (GCGS * m_prim)).sqrt(),
                    0.0,
                )
            } else {
                (0.0, 0.0)
            };

            let total_d_eorb = d_eorb_moon + d_eorb_pl + d_eorb_ring;
            if -d_time * total_d_eorb < world_states[world_idx].e_orb {
                world_states[world_idx].e_orb += d_time * total_d_eorb;
            } else {
                let d_eorb_moon_adjusted = -world_states[world_idx].e_orb / d_time;
                world_states[world_idx].w_tide_tot = (-d_eorb_moon_adjusted
                    + 171.0 / 16.0 * (GCGS / m_prim).sqrt() * r_prim.powi(5) * mass_moon / q_prim
                        * aorb_im.powf(-6.5)
                        * world_states[world_idx].e_orb)
                    * GCGS
                    * m_prim
                    * mass_moon
                    * world_states[world_idx].e_orb
                    / aorb_im;
                world_states[world_idx].w_fluidtide_tot = 0.0;
                world_states[world_idx].e_orb = MIN_ECC;
            }

            let total_d_aorb = d_aorb_moon + d_aorb_pl + d_aorb_ring;
            if -d_time * total_d_aorb < world_states[world_idx].a_orb {
                world_states[world_idx].a_orb += d_time * total_d_aorb;
            } else {
                let name = &world_states[world_idx].name;
                use std::fs::OpenOptions;
                use std::io::Write;
                let path = format!("Outputs/{}_{}_Orbit.csv", world_idx, name);
                let crash_msg = format!(
                    "Orbit: time={} Gyr, -dtime*d_aorb_moon (= {} m) + -dtime*d_aorb_pl (= {} m) + -dtime*d_aorb_ring (= {} m) > aorb = {} m, moon {} crashes into planet",
                    real_time / GYR2SEC,
                    -d_time * d_aorb_moon,
                    -d_time * d_aorb_pl,
                    -d_time * d_aorb_ring,
                    world_states[world_idx].a_orb,
                    world_idx
                );
                println!("{}", crash_msg);
                if let Ok(mut file) = OpenOptions::new().write(true).append(true).open(&path) {
                    let _ = writeln!(file, "{}", crash_msg);
                }
                world_states[world_idx].a_orb = -1.0e-10;
                world_states[world_idx].e_orb = 0.0;
            }

            world_states[world_idx].n_orb = (GCGS * m_prim / world_states[world_idx].a_orb.powi(3))
                .sqrt()
                .max(0.) // we can do this because if the sqrt() returns NaN, then 0 is returned instead.
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
        resonance[im] = vec![0.; nmoons];

        let m_prim = self.primary_world.mass;
        let r_prim = self.primary_world.rad;
        let k2_prim = self.primary_world.k2;
        let prim_sign = if self.worlds[im].retrograde {
            -1.0
        } else {
            1.0
        };
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
                            (a_old[inner] / a_old[outer]).powf(-1.5)
                        } else {
                            (world_states[inner].a_orb / world_states[outer].a_orb).powf(-1.5)
                        } * j
                            / (j + l);

                        let tol_lower = if is_already_captured { 0.985 } else { 0.99 };
                        let tol_upper = if is_already_captured { 1.015 } else { 1.01 };

                        if commensurability > tol_lower && commensurability < tol_upper {
                            resonance[inner][outer] = j;
                            resonance[outer][inner] = j;

                            let dnorb_dt = |idx: usize| {
                                -1.5 * (GCGS * m_prim).sqrt()
                                    * world_states[idx].a_orb.powf(-2.5)
                                    * prim_sign
                                    * if reslock {
                                        world_states[idx].a_orb / t_tide[idx]
                                    } else {
                                        3.0 * k2_prim
                                            * (GCGS / m_prim).sqrt()
                                            * r_prim.powi(5)
                                            * self.worlds[idx].mass()
                                            / q_prim
                                            * world_states[idx].a_orb.powf(-5.5)
                                    }
                            };

                            let dnorb_dt_inner = dnorb_dt(inner);
                            let dnorb_dt_outer = dnorb_dt(outer);

                            if j * dnorb_dt_inner <= (j + l) * dnorb_dt_outer {
                                let p_cap =
                                    self.mmr_capture_probability(world_states, inner, outer, j);
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
    pub fn resscreen(&self, resonance: &[f64], res_acct_for: &mut [f64], res_acct_for_old: &[f64]) {
        let nmoons = resonance.len();
        res_acct_for.fill(0.0);

        let mut res_min = IJMAX as f64;
        let mut nbres = 0;

        // Find the min order of resonance for each moon and the number of moons
        // involved in resonances of this order

        // Copy only the lowest-order resonances for each moon
        for i in 0..nmoons {
            if resonance[i] > 0.0 && resonance[i] <= res_min {
                res_min = resonance[i];
            }
            if resonance[i] == res_min {
                nbres += 1;
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

pub fn mmid(
    y: &[f64],
    dydx: &[f64],
    params: &MmrAvgHamParams,
    xs: f64,
    htot: f64,
    nstep: usize,
    derivs: fn(f64, &[f64], &MmrAvgHamParams) -> Vec<f64>,
) -> Vec<f64> {
    let nv = y.len();
    let mut ym = vec![0.0; nv];
    let mut yn = vec![0.0; nv];

    let h = htot / (nstep as f64);

    for i in 0..nv {
        ym[i] = y[i];
        yn[i] = y[i] + h * dydx[i];
    }

    let mut x = xs + h;
    let mut dy = derivs(x, &yn, params);
    let h2 = 2.0 * h;

    for _n in 1..nstep {
        let mut swap = vec![0.0; nv];
        for i in 0..nv {
            swap[i] = ym[i] + h2 * dy[i];
            ym[i] = yn[i];
            yn[i] = swap[i];
        }
        x += h;
        dy = derivs(x, &yn, params);
    }

    let mut yout = vec![0.0; nv];
    for i in 0..nv {
        yout[i] = 0.5 * (ym[i] + yn[i] + h * dy[i]);
    }
    yout
}

pub struct MmrAvgHamParams {
    pub m: [f64; 2],
    pub r: [f64; 2],
    pub w: [f64; 2],
    pub j: f64,
    pub m_prim: f64,
    pub r_prim: f64,
    pub j2_prim: f64,
    pub j4_prim: f64,
    pub k2_prim: f64,
    pub q_prim: f64,
    pub cs_ee: f64,
    pub cs_eep: f64,
    pub cr_e: f64,
    pub cr_ep: f64,
    pub cr_ee: f64,
    pub cr_eep: f64,
    pub cr_epep: f64,
    pub speedup: f64,
    pub psgn: [f64; 2],
    pub reslock: bool,
    pub ttide: [f64; 2],
    pub alpha_lind: f64,
    pub ring_surface_density: f64,
    pub aring_in: f64,
    pub aring_out: f64,
}

pub fn mmr_avg_ham(_x: f64, y: &[f64], params: &MmrAvgHamParams) -> Vec<f64> {
    let h = [y[0], y[3]];
    let k = [y[1], y[4]];
    let a_ = [y[2], y[5]];

    let mut e2 = [0.0; 2];
    let mut n_ = [0.0; 2];
    let mut sigbar = [0.0; 2];
    let mut lambda = [0.0; 2];
    let mut l_ang = [0.0; 2];
    let mut a = [0.0; 2];
    let mut n = [0.0; 2];

    for im in 0..2 {
        e2[im] = h[im] * h[im] + k[im] * k[im];
        n_[im] = (GCGS * params.m_prim / a_[im].powi(3)).sqrt();
        sigbar[im] = 0.5 * e2[im];
        lambda[im] = (a_[im] * GCGS * params.m[im] * params.m[im] * params.m_prim).sqrt();

        let denom = 1.0 - (1.0 - e2[im]).max(0.0).sqrt();
        let l_ratio = if e2[im] < 1.0e-12 {
            1.0 / (1.0 + 0.25 * e2[im])
        } else {
            sigbar[im] / denom
        };
        l_ang[im] = lambda[im] * l_ratio;
        a[im] = l_ang[im] * l_ang[im] / (GCGS * params.m[im] * params.m[im] * params.m_prim);
        n[im] = (GCGS * params.m_prim / a[im].powi(3)).sqrt();
    }

    let dhk = [(1.0 - params.j) * n[0] + params.j * n[1]; 2];

    let mut delta_n = [0.0; 2];
    let mut omdot = [0.0; 2];
    for im in 0..2 {
        let r_a2 = (params.r_prim / a_[im]).powi(2);
        let r_a4 = r_a2 * r_a2;
        delta_n[im] = n_[im]
            * (3.0 * params.j2_prim * r_a2
                + (45.0 / 4.0 * params.j2_prim * params.j2_prim - 15.0 / 4.0 * params.j4_prim)
                    * r_a4);
        omdot[im] = n_[im]
            * (1.5 * params.j2_prim * r_a2
                + (63.0 / 8.0 * params.j2_prim * params.j2_prim - 15.0 / 4.0 * params.j4_prim)
                    * r_a4);
    }

    let mut delta_sigdot = [0.0; 2];
    for im in 0..2 {
        delta_sigdot[im] = (1.0 - params.j) * delta_n[0] + params.j * delta_n[1] - omdot[im];
    }

    let mut c = [0.0; 2];
    let mut d = [0.0; 2];
    let mut dk_tide = [0.0; 2];
    let mut dh_tide = [0.0; 2];
    let mut da_tide = [0.0; 2];

    for im in 0..2 {
        c[im] = 3.0 * params.k2_prim / params.q_prim * params.m[im] / params.m_prim
            * (GCGS * params.m_prim).sqrt()
            * params.r_prim.powi(5);
        d[im] = params.w[im] / (params.k2_prim / params.q_prim)
            * (params.m_prim / params.m[im]).powi(2)
            * (params.r[im] / params.r_prim).powi(5);
        dk_tide[im] = -3.5 * c[im] * d[im] * a[im].powf(-6.5) * k[im] * params.speedup;
        dh_tide[im] = -3.5 * c[im] * d[im] * a[im].powf(-6.5) * h[im] * params.speedup;

        if !params.reslock {
            da_tide[im] = c[im]
                * (params.psgn[im] - 7.0 * d[im] * e2[im])
                * a[im].powf(-5.5)
                * params.speedup;
        } else {
            let t_tide_val = params.ttide[im];
            let first_term = if t_tide_val > 0.0 {
                params.psgn[im] * a[im] / t_tide_val
            } else {
                0.0
            };
            da_tide[im] =
                (first_term - 7.0 * c[im] * d[im] * e2[im] * a[im].powf(-5.5)) * params.speedup;
        }

        if params.ring_surface_density > 0.0 {
            let mut ring_torque = 0.0;
            let den_in = 1.0 - (params.aring_in / a[im]).powf(1.5);
            let den_out = 1.0 - (params.aring_out / a[im]).powf(1.5);
            if den_in.abs() > 1.0e-30 && den_out.abs() > 1.0e-30 {
                let kmin = (1.0 / den_in).floor() as i32 + 1;
                let kmax = (1.0 / den_out).floor() as i32;
                let alpha_lind: f64 = if params.ring_surface_density <= 2.0 {
                    2.0e-5
                } else {
                    1.0e-4
                };
                let kmax_limit = (1.0 / alpha_lind.sqrt()).floor() as i32;
                if kmin <= kmax && kmax <= kmax_limit {
                    for i in kmin..=kmax {
                        ring_torque += std::f64::consts::PI * std::f64::consts::PI / 3.0
                            * params.ring_surface_density
                            * GCGS
                            * params.m[im]
                            * params.m[im]
                            * (i as f64)
                            * ((i - 1) as f64)
                            * a[im]
                            / params.m_prim;
                    }
                }
            }
            da_tide[im] += 2.0 * ring_torque / params.m[im]
                * (a[im] / (GCGS * params.m_prim)).sqrt()
                * params.speedup;
        }
    }

    let mut dk = [0.0; 2];
    let mut dh = [0.0; 2];

    dk[0] = (dhk[0] + delta_sigdot[0]) * h[0]
        - GCGS * params.m[0] * params.m[1] / (a_[1] * lambda[0])
            * (params.cs_eep * h[1]
                + 2.0 * params.cs_ee * h[0]
                + params.cr_e
                + 2.0 * params.cr_ee * h[0]
                + params.cr_eep * h[1])
        + dk_tide[0];
    dh[0] = (-dhk[0] - delta_sigdot[0]) * k[0]
        - GCGS * params.m[0] * params.m[1] / (a_[1] * lambda[0])
            * (-params.cs_eep * k[1] - 2.0 * params.cs_ee * k[0]
                + 2.0 * params.cr_ee * k[0]
                + params.cr_eep * k[1])
        + dh_tide[0];
    dk[1] = (dhk[1] + delta_sigdot[1]) * h[1]
        - GCGS * params.m[0] * params.m[1] / (a_[1] * lambda[1])
            * (2.0 * params.cs_ee * h[1]
                + params.cs_eep * h[0]
                + params.cr_ep
                + 2.0 * params.cr_epep * h[1]
                + params.cr_eep * h[0])
        + dk_tide[1];
    dh[1] = (-dhk[1] - delta_sigdot[1]) * k[1]
        - GCGS * params.m[0] * params.m[1] / (a_[1] * lambda[1])
            * (-2.0 * params.cs_ee * k[1] - params.cs_eep * k[0]
                + 2.0 * params.cr_epep * k[1]
                + params.cr_eep * k[0])
        + dh_tide[1];

    let mut dl_tide = [0.0; 2];
    let mut d_sigbar_tide = [0.0; 2];
    for im in 0..2 {
        dl_tide[im] = l_ang[im] / (2.0 * a[im]) * da_tide[im];
        d_sigbar_tide[im] = h[im] * dh_tide[im] + k[im] * dk_tide[im];
    }

    let mut d_lambda_tide = [0.0; 2];
    d_lambda_tide[0] = (1.0 + params.j * sigbar[1]) * dl_tide[0]
        - (1.0 - params.j) * lambda[0] * d_sigbar_tide[0]
        - (1.0 - params.j) * lambda[1] * d_sigbar_tide[1]
        - (1.0 - params.j) * dl_tide[1] * sigbar[1];
    d_lambda_tide[1] = (1.0 + (1.0 - params.j) * sigbar[0]) * dl_tide[1]
        - params.j * lambda[0] * d_sigbar_tide[0]
        - params.j * lambda[1] * d_sigbar_tide[1]
        - params.j * dl_tide[0] * sigbar[0];

    let denom = 1.0 + (1.0 - params.j) * sigbar[0] + params.j * sigbar[1];
    d_lambda_tide.iter_mut().for_each(|n| {
        *n /= denom;
    });

    let da_ = [0, 1].map(|im| 2.0 * a_[im] / lambda[im] * d_lambda_tide[im]);

    vec![dh[0], dk[0], da_[0], dh[1], dk[1], da_[1]]
}
