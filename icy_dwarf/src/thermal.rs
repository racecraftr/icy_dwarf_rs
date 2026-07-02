use std::{f64::consts::FRAC_PI_2, fs, process::exit};

use crate::{
    FloatMat,
    consts::{
        CM, E_YOUNG_OLIV, E_YOUNG_SERP, GCGS, GRAM, KM, KM2CM, NU_POISSON_OLIV, NU_POISSON_SERP,
        PA2BA, PI_GREEK,
    },
    input::{Fracs, IcyDwarfInput},
    planet_system::{RHO_ADHS_TH, RHO_H2OL_TH, RHO_H2OS_TH, RHO_NH3L_TH, WorldState, ZoneState},
    to_faer_mat,
};
use faer::{Mat, linalg::solvers::DenseSolveCore, prelude::Solve};
use itertools::Itertools;
use num::{
    complex::{Complex64, ComplexFloat},
    pow::Pow,
    traits::Inv,
};

const K: f64 = 200.0e9 / GRAM * CM;

impl IcyDwarfInput {
    pub fn thermal(&self, world_states: &mut [WorldState], dtime: f64) {
        let tilt_scale: f64 = 1.0;

        for world_state in world_states.iter_mut() {
            // 1. Calculate Pressure
            self.calculate_pressure(world_state);

            // 2. Update Porosity & Radii (Creep & Compaction)
            self.update_porosity(world_state, dtime);

            // Further thermal logic to be implemented
        }
    }

    // TODO: finish this
    pub fn tide(&self, world_state: &mut WorldState) {
        const D_EPS: f64 = 2.22e-16;
        let base_vec = vec![0_f64; world_state.zones.len()];
        let mut shearmod = vec![Complex64::ZERO; world_state.zones.len()];
        let mut rho = base_vec.clone();

        let mu_rigid_ice = 4.0e9 / GRAM * CM;

        // accumulated mass throughout all the zones.
        let mut m_acc = base_vec.clone();
        let mut g_vec = base_vec.clone();
        drop(base_vec);
        let alpha_andr = 0.3;
        for (i, zone) in world_state.zones.iter().enumerate() {
            rho[i] = zone.mass_total / zone.volumes().0; // density is just mass over volume
            m_acc[i] = zone.mass_total + if i == 0 { 0. } else { m_acc[i - 1] };
            g_vec[i] = GCGS * m_acc[i] / (zone.radius + zone.dr).powi(2);
            let z = if i < world_state.zones.len() - 1 {
                zone
            } else {
                &world_state.zones[world_state.zones.len() - 2]
            };
            let mut mu_visc = PA2BA * z.pressure * z.creep();
            if zone.mass_ice > 0.
                && zone.mass_ammonia_solid + zone.mass_ammonia_solid >= 0.01 * zone.mass_ice
                && zone.temp > 140.
            {
                mu_visc = (mu_visc
                    * if z.temp < 176. {
                        1.0e-2
                    } else if z.temp < 250. {
                        1.0e-6
                    } else if z.temp < 271. {
                        1.0e-14
                    } else {
                        1.
                    })
                .max(1.0e3);
            }
            let mu_rigid_rock = (z.x_hydr * E_YOUNG_SERP / (2. * (1. + NU_POISSON_SERP))
                + (1. - z.x_hydr) * E_YOUNG_OLIV / (2. * (1. + NU_POISSON_OLIV)))
                / GRAM
                * CM;
            let mut mu_rigid =
                if z.mass_ice + z.mass_ammonia_solid + z.mass_water + z.mass_ammonia_liquid > 0. {
                    let phi = 1. - zone.fracs().0;
                    mu_rigid_ice
                        * if phi < 0.3 {
                            mu_rigid_rock * 0.3_f64.exp()
                        } else {
                            1.
                        }
                } else {
                    mu_rigid_rock
                };

            if z.mass_water + z.mass_ammonia_liquid > 0.9 * z.mass_total {
                mu_visc = 1.0e2 * PA2BA;
                mu_rigid = 1.03 * PA2BA;
            }
            let gamma_andr = match alpha_andr {
                0.2 => 0.918169,
                0.3 => 0.897471,
                0.4 => 0.887264,
                0.5 => 0.886227,
                _ => exit(0),
            };
            // let cond = world_state.omega.abs() < 100. * D_EPS;
            let cond_i = |n: f64| {
                if world_state.omega.abs() < 100. * D_EPS {
                    Complex64::ZERO
                } else {
                    Complex64::I * n
                }
            };
            shearmod[i] = match self.world_spec.rhelogy {
                crate::input::TidalModel::Maxwell => {
                    mu_rigid * world_state.omega.powi(2) * mu_visc.powi(2)
                        / (mu_rigid.powi(2) + (world_state.omega * mu_visc).powi(2))
                        + cond_i(
                            mu_rigid.powi(2) * world_state.omega * mu_visc
                                / (mu_rigid.powi(2) + (world_state.omega * mu_visc).powi(2)),
                        )
                }
                crate::input::TidalModel::Burgers => {
                    let mu2 = 0.02 * mu_visc;
                    let c_1 = 2. / mu_rigid + mu2 / (mu_rigid * mu_visc);
                    let c_2 = 1. / mu_visc + mu2 * (world_state.omega / mu_rigid).powi(2);
                    let d_burgers = c_1 * c_2 + c_2 * c_2 + world_state.omega.powi(2);
                    world_state.omega.powi(2) * (c_1 - mu2 * c_2 / mu_rigid) / d_burgers
                        + cond_i(
                            (c_2 + mu2 * world_state.omega.powi(2) * c_1 / mu_rigid) / d_burgers,
                        )
                }
                crate::input::TidalModel::Andr => {
                    let beta_andr = 1.0 / (mu_rigid * (mu_visc / mu_rigid).powf(alpha_andr));
                    let a_andr = mu_rigid.inv()
                        + world_state.omega.powf(-alpha_andr)
                            * beta_andr
                            * (alpha_andr * FRAC_PI_2).cos();
                    let b_andr = 1.0 / (mu_visc * world_state.omega)
                        + world_state.omega.powf(-alpha_andr)
                        + beta_andr * (alpha_andr * FRAC_PI_2).sin();
                    let d_andr = a_andr.powi(2) + b_andr.powi(2);
                    a_andr / d_andr + cond_i(b_andr / d_andr)
                }
                crate::input::TidalModel::SunCoop => {
                    let (voigt_comp_offset, voigt_visc_offset, zeta_andr) = (0.43, 0.02, 1.);
                    let comp_maxwell = mu_rigid.inv();
                    let comp_voigt = voigt_comp_offset * comp_maxwell;
                    let visc_voigt = voigt_visc_offset * mu_visc;
                    let sine_andr =
                        (alpha_andr * FRAC_PI_2).cos() + cond_i((alpha_andr * FRAC_PI_2).sin());
                    let c_comp_maxwell = comp_maxwell + cond_i(1. / (world_state.omega * mu_visc));
                    let c_comp_sub_andr = comp_maxwell
                        * (world_state.omega * comp_maxwell * mu_visc * zeta_andr).pow(-alpha_andr)
                        * sine_andr;
                    let c_comp_voigt = cond_i(comp_voigt.powi(2) * visc_voigt * world_state.omega)
                        * (comp_voigt.powi(2)
                            + visc_voigt.powi(2)
                            + world_state.omega.powi(2)
                            + 1.)
                            .inv();

                    (c_comp_maxwell + c_comp_sub_andr + c_comp_voigt).inv()
                }
            }
        }

        let y_tide = prop_mtx(
            &world_state
                .zones
                .iter()
                .map(|x| x.radius)
                .collect::<Vec<_>>(),
            &rho,
            &g_vec,
            &shearmod,
            0,
        );

        let (e2, e4, e6, e8, e10) = world_state.ecc();
        let eterm_1 = e10 * (2555911.0 / 122880.0) - e8 * (63949.0 / 2304.0) + e6 * (551.0 / 12.0)
            - e4 * (101.0 / 4.0)
            + e2 * 7.0;
        let eterm_2 = e10 * (-171083.0 / 320.0) + e8 * (339187.0 / 576.0) - e6 * (3847.0 / 12.0)
            + e4 * (605.0 / 8.0);
        let eterm_3 =
            e10 * (368520907.0 / 81920.0) - e8 * (1709915.0 / 768.0) + e6 * (2855.0 / 6.0);
        let eterm_4 = e10 * (-66268493.0 / 5760.0) + e8 * (2592379.0 / 1152.0);
        let eterm_5 = e10 * (6576742601.0 / 737280.0);
        let eterm = match self.world_spec.ecc_model {
            crate::input::EccModel::E2 => e2,
            crate::input::EccModel::E10Cpl => eterm_1 + eterm_2 + eterm_3 + eterm_4 + eterm_5,
            crate::input::EccModel::E10Ctl => {
                eterm_1 + 2. * eterm_2 + 3. * eterm_3 + 4. * eterm_4 + 5. * eterm_5
            }
        } / 7.;

        for (idx, zone) in world_state.zones.iter().enumerate().skip(1) {
            let r_out = zone.radius + zone.dr;
            let x = 2. * y_tide[idx][0] - 6. * y_tide[idx][1];
            let h_mu = 4. / 3.
                * (r_out / (K + 4. / 3. * shearmod[idx]).abs()
                    * (y_tide[idx][2] - (K - 2. / 3. * shearmod[idx]) / r_out * x)
                        .abs()
                        .powi(2)
                    - r_out * ((y_tide[idx][0].conj() - y_tide[idx][1].conj()) / r_out * x).re
                    + 1. / 3. * x.abs().powi(2)
                    + 6. * (r_out * y_tide[idx][3].abs() / shearmod[idx].abs()).powi(2)
                    + 24. * y_tide[idx][1].abs().powi(2));

            let w_tide = if self.world_spec.tidal_heating {
                zone.volumes().0
                    * 2.
                    * world_state.omega.powi(5)
                    * world_state.zones.last().unwrap().radius.powi(4)
                    * (eterm + world_state.obl.sin() / 7.)
                    / r_out.powi(2)
                    * h_mu
                    * shearmod[idx].im
            } else {
                0.
            };
            world_state.w_tide_tot += w_tide;
        }
    }

    pub fn calculate_pressure(&self, world: &mut WorldState) {
        let mut cumulative_mass = 0.0;
        let mut gravity = vec![0.0; self.grid.n_zones];

        for (i, zone) in world.zones.iter().enumerate() {
            cumulative_mass += zone.mass_total;
            gravity[i] =
                GCGS * cumulative_mass * GRAM / zone.radius.powi(2) * KM2CM.powi(2) / KM.powi(2);
        }

        if let Some(last) = world.zones.last_mut() {
            last.pressure = 0.0;
        }

        for ir in (0..self.grid.n_zones - 1).rev() {
            let next_g = gravity[ir + 1];
            let curr_g = gravity[ir];
            let next_r = world.zones[ir + 1].radius;
            let curr_r = world.zones[ir].radius;
            let next_pressure = world.zones[ir + 1].pressure;

            let next_zone = &world.zones[ir + 1];
            let Fracs(frock, fh2os, fadhs, fh2ol, fnh3l) = next_zone.fracs();

            let avg_density = frock
                * (next_zone.x_hydr * self.world_spec.rho_hydr_th()
                    + (1.0 - next_zone.x_hydr) * self.world_spec.rho_rock_th())
                + fh2os * RHO_H2OS_TH
                + fh2ol * RHO_H2OL_TH
                + fadhs * RHO_ADHS_TH
                + fnh3l * RHO_NH3L_TH;

            world.zones[ir].pressure = next_pressure
                + 0.5 * (next_g + curr_g) * (next_r - curr_r) / crate::consts::KM2CM
                    * crate::consts::KM
                    * avg_density;
        }
    }

    fn update_porosity(&self, world: &mut WorldState, dtime: f64) {
        let mut r_old = vec![0.0; self.grid.n_zones + 1];
        r_old[0] = 0.0;
        for (ir, zone) in world.zones.iter().enumerate() {
            r_old[ir + 1] = zone.radius;
        }
        let mut pore_old = vec![0.0; self.grid.n_zones];
        for (ir, zone) in world.zones.iter().enumerate() {
            pore_old[ir] = zone.porosity;
        }

        for zone in world.zones.iter_mut() {
            let v_rock = zone.mass_rock
                / (zone.x_hydr * self.world_spec.rho_hydr_th()
                    + (1.0 - zone.x_hydr) * self.world_spec.rho_rock_th());
            let d_vol = zone.volumes().0;

            let x_ice = 1.0 - v_rock / d_vol;
            let c_rate = creep(zone.temp, zone.pressure, x_ice, zone.porosity, zone.x_hydr);
            zone.porosity -= dtime * (1.0 - zone.porosity) * c_rate;
            if zone.porosity < 0. || zone.mass_rock < 0.01 && zone.mass_water > 0.01 {
                zone.porosity = 0.0;
            }
        }

        // Update radii
        let mut current_r: f64 = 0.0;
        for ir in 0..self.grid.n_zones {
            let zone = &mut world.zones[ir];
            let next_r = (current_r.powi(3)
                + (r_old[ir + 1].powi(3) - r_old[ir].powi(3)) * (1.0 - pore_old[ir])
                    / (1.0 - zone.porosity))
                .powf(1.0 / 3.0);
            zone.radius = next_r;
            current_r = next_r;
        }
    }
}

impl WorldState {
    pub fn ecc(&self) -> (f64, f64, f64, f64, f64) {
        (1..=5)
            .map(|n| self.e_orb.powi(n * 2))
            .next_tuple()
            .unwrap()
    }
}

impl ZoneState {
    pub fn creep(&self) -> f64 {
        creep(
            self.temp,
            self.pressure,
            1.0 - self.fracs().0,
            self.porosity,
            self.x_hydr,
        )
    }
}

impl ZoneState {
    pub fn viscosity(&self) -> f64 {
        let x = self.mass_ammonia_liquid / self.mass_water;
        let (a, b) = if self.temp > 240.0 {
            [
                (-10.8143, 1819.86),
                (0.711062, 250.822),
                (-22.4943, 6505.25),
                (41.8343, 14923.4),
                (18.5149, 7141.76),
            ]
        } else {
            [
                (-13.8628, 2701.73),
                (-68.7617, 14973.3),
                (230.038, -46174.5),
                (-249.897, 45967.7),
                (0., 0.),
            ]
        }
        .iter()
        .enumerate()
        .fold((0., 0.), |(a, b), (k, (a_t, b_t))| {
            let x_pow = x.powi(k as i32);
            (a + a_t * x_pow, b + b_t * x_pow)
        });
        (a + b / self.temp).exp()
    }
}

fn creep(t: f64, p: f64, x_ice: f64, porosity: f64, x_hydr: f64) -> f64 {
    use crate::consts::{D_FLOW_LAW, MPA, R_G};
    let eff_p = p / MPA / (1.0 - porosity);

    let eps_disl = 4.0e5 * eff_p.powi(4) * (-60.0e3 / (R_G * t)).exp();

    let eps_basal = if t < 255.0 {
        3.9e-3 * eff_p.powf(1.8) * D_FLOW_LAW.powf(-1.4) * (-49.0e3 / (R_G * t)).exp()
    } else {
        3.0e26 * eff_p.powf(1.8) * D_FLOW_LAW.powf(-1.4) * (-192.0e3 / (R_G * t)).exp()
    };

    let eps_gbs = 5.5e7 * eff_p.powf(2.4) * (-60.0e3 / (R_G * t)).exp();
    let eps_diff = 3.02e-14 * eff_p.powi(1) * D_FLOW_LAW.powi(-2) * (-59.4e3 / (R_G * t)).exp();

    let creep_rate_ice = eps_diff + 1.0 / (1.0 / eps_basal + 1.0 / eps_gbs) + eps_disl;

    creep_rate_ice
        * if x_ice > 0.3 {
            1.
        } else {
            let t_eff = t.max(140.);
            let creep_rate_hydr =
                416869.3834703355 * eff_p * D_FLOW_LAW.powi(-3) * (-240.0e3 / (R_G * t_eff)).exp();
            let creep_rate_dry = 177827.94100389228
                * eff_p
                * D_FLOW_LAW.powf(-2.98)
                * (-261.0e3 / (R_G * t_eff)).exp();

            (x_hydr * creep_rate_hydr + (1. - x_hydr) * creep_rate_dry)
                * (0.3 + 7. / 3. * x_ice).exp()
        }
}

pub fn prop_mtx(
    r: &[f64],
    rho: &[f64],
    g: &[f64],
    shearmod: &[Complex64],
    ir_core: usize,
) -> Vec<Vec<Complex64>> {
    let nr = rho.len();
    if nr == 0 {
        return Vec::new();
    }
    assert!(r.len() >= nr);
    assert!(g.len() >= nr);
    assert!(shearmod.len() >= nr);

    if ir_core >= nr {
        return vec![vec![Complex64::ZERO; 6]; nr];
    }

    let mut ypropmtx = vec![[[Complex64::ZERO; 6]; 6]; nr];
    let mut ypropinv = vec![[[Complex64::ZERO; 6]; 6]; nr];
    let mut bpropmtx = vec![[[Complex64::ZERO; 3]; 6]; nr];

    for ir in ir_core..nr {
        let r_val = r[ir + 1];
        let r_val_2 = r_val * r_val;
        let r_val_3 = r_val_2 * r_val;
        let r_val_4 = r_val_3 * r_val;
        let r_val_5 = r_val_4 * r_val;

        let rho_g_r = Complex64::from(rho[ir] * g[ir] * r_val);
        let sm = shearmod[ir];
        let four_pi_g_rho = Complex64::from(4.0 * PI_GREEK * GCGS * rho[ir]);
        let two_pi_g_rho = Complex64::from(2.0 * PI_GREEK * GCGS * rho[ir]);

        let rho_g_r_over_sm = rho_g_r / sm;
        let r_over_sm = Complex64::from(r_val) / sm;

        ypropmtx[ir] = [
            [
                Complex64::from(r_val_3 / 7.0),
                Complex64::from(r_val),
                Complex64::ZERO,
                Complex64::from(1.0 / (2.0 * r_val_2)),
                Complex64::from(1.0 / r_val_4),
                Complex64::ZERO,
            ],
            [
                Complex64::from(5.0 * r_val_3 / 42.0),
                Complex64::from(r_val / 2.0),
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::from(-1.0 / (3.0 * r_val_4)),
                Complex64::ZERO,
            ],
            [
                (rho_g_r - sm) * Complex64::from(r_val_2 / 7.0),
                rho_g_r + sm * 2.0,
                Complex64::from(-rho[ir] * r_val_2),
                (rho_g_r - sm * 6.0) / Complex64::from(2.0 * r_val_3),
                (rho_g_r - sm * 8.0) / Complex64::from(r_val_5),
                Complex64::from(-rho[ir] / r_val_3),
            ],
            [
                sm * Complex64::from(8.0 * r_val_2 / 21.0),
                sm,
                Complex64::ZERO,
                sm / Complex64::from(2.0 * r_val_3),
                sm * Complex64::from(8.0 / (3.0 * r_val_5)),
                Complex64::ZERO,
            ],
            [
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::from(-r_val_2),
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::from(-1.0 / r_val_3),
            ],
            [
                four_pi_g_rho * Complex64::from(r_val_3 / 7.0),
                four_pi_g_rho * Complex64::from(r_val),
                Complex64::from(-5.0 * r_val),
                two_pi_g_rho / Complex64::from(r_val_2),
                four_pi_g_rho / Complex64::from(r_val_4),
                Complex64::ZERO,
            ],
        ];

        ypropinv[ir] = [
            [
                rho_g_r_over_sm - 8.0,
                Complex64::from(16.0),
                -r_over_sm,
                r_over_sm * 2.0,
                Complex64::from(rho[ir]) * r_over_sm,
                Complex64::ZERO,
            ],
            [
                -rho_g_r_over_sm + 6.0,
                Complex64::from(-6.0),
                r_over_sm,
                Complex64::ZERO,
                -Complex64::from(rho[ir]) * r_over_sm,
                Complex64::ZERO,
            ],
            [
                four_pi_g_rho,
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::from(-1.0),
            ],
            [
                rho_g_r_over_sm + 2.0,
                Complex64::from(6.0),
                -r_over_sm,
                r_over_sm * -3.0,
                Complex64::from(rho[ir]) * r_over_sm,
                Complex64::ZERO,
            ],
            [
                -rho_g_r_over_sm + 1.0,
                Complex64::from(-16.0),
                r_over_sm,
                r_over_sm * 5.0,
                -Complex64::from(rho[ir]) * r_over_sm,
                Complex64::ZERO,
            ],
            [
                four_pi_g_rho * Complex64::from(r_val),
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::ZERO,
                Complex64::from(5.0),
                Complex64::from(-r_val),
            ],
        ];

        for j in 0..6 {
            ypropinv[ir][0][j] *= Complex64::from(3.0 / (5.0 * r_val_3));
            ypropinv[ir][1][j] *= Complex64::from(1.0 / (5.0 * r_val));
            ypropinv[ir][2][j] *= Complex64::from(1.0 / (5.0 * r_val));
            ypropinv[ir][3][j] *= Complex64::from(2.0 * r_val_2 / 5.0);
            ypropinv[ir][4][j] *= Complex64::from(3.0 * r_val_4 / 35.0);
            ypropinv[ir][5][j] *= Complex64::from(-r_val_3 / 5.0);
        }
    }

    bpropmtx[ir_core][2][0] = Complex64::from(1.0);
    bpropmtx[ir_core][3][1] = Complex64::from(1.0);
    bpropmtx[ir_core][5][2] = Complex64::from(1.0);

    for ir in (ir_core + 1)..nr {
        let mut btemp = [[Complex64::ZERO; 3]; 6];
        for i in 0..6 {
            for j in 0..3 {
                let mut sum = Complex64::ZERO;
                for k in 0..6 {
                    sum += ypropinv[ir - 1][i][k] * bpropmtx[ir - 1][k][j];
                }
                btemp[i][j] = sum;
            }
        }
        for i in 0..6 {
            for j in 0..3 {
                let mut sum = Complex64::ZERO;
                for k in 0..6 {
                    sum += ypropmtx[ir][i][k] * btemp[k][j];
                }
                bpropmtx[ir][i][j] = sum;
            }
        }
    }

    let mbc = [2, 3, 5]
        .iter()
        .map(|&idx| bpropmtx[nr - 1][idx][0..=2].to_vec())
        .collect::<Vec<_>>();
    let bsurf = vec![
        Complex64::ZERO,
        Complex64::ZERO,
        Complex64::from(-5.0 / r[nr - 1]),
    ];

    let mut ytide = vec![vec![Complex64::ZERO; 6]; nr];

    if let Some((x, _)) = gauss_jordan(&mbc, &bsurf) {
        for ir in 0..nr {
            for i in 0..6 {
                let mut sum = Complex64::ZERO;
                for j in 0..3 {
                    sum += bpropmtx[ir][i][j] * x[(j, 0)];
                }
                ytide[ir][i] = sum;
            }
        }
    } else {
        eprintln!("Thermal: Singular matrix in GaussJordan (prop_mtx)");
        exit(1);
    }

    ytide
}

/// Calculates the single value decomposition of a mxn matrix.
/// Retunrs, in order: U, Sigma, V^*.
fn svd(mat: &[Vec<f64>]) -> Option<(FloatMat, FloatMat, FloatMat)> {
    let rows = mat.len();
    let cols = mat[0].len();
    let svd = to_faer_mat(mat)?.svd().ok()?;
    let (u, s, v) = (svd.U().to_owned(), svd.S(), svd.V().to_owned());
    let mut sigma = Mat::<f64>::zeros(rows, cols);
    for (i, &n) in s.column_vector().iter().enumerate() {
        sigma[(i, i)] = n;
    }
    Some((u, sigma, v))
}

/// Gauss jordan elimination on a complex matrix.
/// Given a matrix and a vector,
/// Uses full pivoting to solve Ax = b.
/// Returns both x and the inverse of A.
fn gauss_jordan(
    mat: &[Vec<Complex64>],
    b: &[Complex64],
) -> Option<(Mat<Complex64>, Mat<Complex64>)> {
    let rows = mat.len();
    let cols = mat[0].len();
    if rows != cols || b.len() != cols {
        return None; // avoids unneccessary computation
    }
    let a = Mat::from_fn(rows, cols, |i, j| mat[i][j]);
    if a.determinant() == Complex64::from(0.0) {
        // matrix is not invertible
        return None;
    }
    let b = Mat::from_fn(b.len(), 1, |i, _| b[i]);
    let lu = a.full_piv_lu();
    let x = lu.solve(b);
    let a_inv = lu.inverse();
    Some((x, a_inv))
}

#[derive(Clone, Debug, Default)]
pub struct ThermalOut {
    pub radius_km: f64,
    pub temp_kelvin: f64,
    pub mass_rock: f64,
    pub mass_ice: f64,
    pub mass_ammonia_solid: f64,
    pub mass_water: f64,
    pub mass_ammonia_liquid: f64,
    pub nusselt_num: f64,
    pub ice_frac_amorphous: f64,
    pub thermal_cond: f64,
    pub deg_of_hydr: f64,
    pub porosity: f64,
    pub crack: bool,
    pub tidal_heating_rate: f64,
}

type ThermVol = (f64, (f64, f64, f64, f64, f64));

impl ThermalOut {
    pub fn from_line(ln: &str) -> Option<Self> {
        let parts = ln.split_whitespace().collect::<Vec<_>>();
        let radius_km = parts[0].parse::<f64>().ok()?;
        let radius_km = radius_km * KM2CM;
        Some(Self {
            radius_km,
            temp_kelvin: parts[1].parse().ok()?,
            mass_rock: parts[2].parse().ok()?,
            mass_ice: parts[3].parse().ok()?,
            mass_ammonia_solid: parts[4].parse().ok()?,
            mass_water: parts[5].parse().ok()?,
            mass_ammonia_liquid: parts[6].parse().ok()?,
            nusselt_num: parts[7].parse().ok()?,
            ice_frac_amorphous: parts[8].parse().ok()?,
            thermal_cond: parts[9].parse().ok()?,
            deg_of_hydr: parts[10].parse().ok()?,
            porosity: parts[11].parse().ok()?,
            crack: parts[12].parse::<u8>().map(|n| n == 1).ok()?,
            tidal_heating_rate: parts[13].parse().ok()?,
        })
    }

    pub fn mass_total(&self) -> f64 {
        self.mass_rock
            + self.mass_ice
            + self.mass_ammonia_solid
            + self.mass_ammonia_liquid
            + self.mass_water
    }

    pub fn vol(&self, input: &IcyDwarfInput) -> (f64, Fracs) {
        let vol_rock = self.mass_rock
            / (self.deg_of_hydr * input.world_spec.rho_hydr_th()
                + (1.0 - self.deg_of_hydr) * input.world_spec.rho_rock_th());
        let vol_ice = self.mass_ice / RHO_H2OL_TH;
        let vol_adhs = self.mass_ammonia_solid / RHO_ADHS_TH;
        let vol_water = self.mass_water / RHO_H2OL_TH;
        let vol_nh3l = self.mass_ammonia_liquid / RHO_NH3L_TH;
        (
            vol_rock + vol_adhs + vol_water + vol_nh3l,
            Fracs(vol_rock, vol_ice, vol_adhs, vol_water, vol_nh3l),
        )
    }

    pub fn fracs(&self) -> Fracs {
        let mass_total = self.mass_total();
        Fracs(
            self.mass_rock / mass_total,
            self.mass_ice / mass_total,
            self.mass_ammonia_solid / mass_total,
            self.mass_water / mass_total,
            self.mass_ammonia_liquid / mass_total,
        )
    }
}

fn read_thermal_out(path: &str, output_time_step: usize) -> Option<Vec<Vec<ThermalOut>>> {
    let Ok(lines) =
        fs::read_to_string(path).map(|s| s.lines().map(str::to_owned).collect::<Vec<_>>())
    else {
        return None;
    };
    Some(
        lines
            .chunks(output_time_step)
            // we don't need to define NT
            // as a pamaeter, as it is sized dynamically.
            .map(|chunk| {
                chunk
                    .iter()
                    .filter_map(|ln| ThermalOut::from_line(ln))
                    .collect()
            })
            .collect(),
    )
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_svd() {
        // matrix from https://en.wikipedia.org/wiki/Singular_value_decomposition.
        let mat: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0, 0.0, 2.0],
            vec![0.0, 0.0, 3.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0, 0.0, 0.0],
        ];
        let Some((u, sigma, v)) = svd(&mat) else {
            panic!("No SVD available for {:?}", mat);
        };

        let (rows, cols) = (mat.len(), mat[0].len());
        let rec = (u.as_ref()) * (sigma.as_ref()) * (v.adjoint().as_ref());
        println!("Reconstructed matrix = {:?}", &rec);
        println!("Actual matrix = {:?}", &mat);

        for i in 0..rows {
            for j in 0..cols {
                assert!((rec[(i, j)] - mat[i][j]).abs() < 1e-10);
            }
        }

        let u_ortho = u.adjoint() * u.as_ref();
        println!("U^* x U = {:?}", &u_ortho);
        let v_ortho = v.adjoint() * v.as_ref();
        println!("V^* x V = {:?}", &v_ortho);

        for i in 0..rows {
            for j in 0..rows {
                let expected = if i == j { 1. } else { 0. };
                assert!((u_ortho[(i, j)] - expected).abs() < 1e-10);
            }
        }
        for i in 0..cols {
            for j in 0..cols {
                let expected = if i == j { 1. } else { 0. };
                assert!((v_ortho[(i, j)] - expected).abs() < 1e-10);
            }
        }
    }
}
