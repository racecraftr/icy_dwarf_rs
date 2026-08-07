use std::{f64::consts::FRAC_PI_2, process::exit};

use num::{
    complex::{Complex64, ComplexFloat},
    pow::Pow,
    traits::Inv,
};

use crate::{
    consts::*,
    input::IcyDwarfInput,
    planet_system::{RHO_H2OL_TH, WorldState},
    thermal::{C0, K, ThermalWorkspace, prop_mtx},
};

impl IcyDwarfInput {
    pub fn tide(&self, world_state: &mut WorldState, workspace: &mut ThermalWorkspace) {
        const D_EPS: f64 = 2.22e-16;
        let n_zones = world_state.zones.len();
        let shearmod = &mut workspace.shearmod;
        shearmod.fill(C0);
        let rho = &mut workspace.rho;
        rho.fill(0.0);
        let m_acc = &mut workspace.m_acc;
        m_acc.fill(0.0);
        let g_vec = &mut workspace.g_vec;
        g_vec.fill(0.0);

        const MU_RIGID_ICE: f64 = 4.0e9 / GRAM * CM;

        let alpha_andr = 0.3;
        for (i, zone) in world_state.zones.iter().enumerate() {
            rho[i] = zone.mass_total / zone.volumes().0; // density is just mass over volume
            m_acc[i] = zone.mass_total + if i == 0 { 0. } else { m_acc[i - 1] };
            g_vec[i] = GCGS * m_acc[i] / (zone.radius + zone.dr).powi(2);
            let z = if i < n_zones - 1 {
                zone
            } else {
                &world_state.zones[n_zones - 2]
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
                    MU_RIGID_ICE
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
                    C0
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
                            * (alpha_andr * FRAC_PI_2).cos()
                            * gamma_andr;
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
                    let sine_andr = (alpha_andr * FRAC_PI_2).cos()
                        + cond_i((alpha_andr * FRAC_PI_2).sin()) * gamma_andr;
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

        let r_grid = &mut workspace.r_grid;
        for (i, zone) in world_state.zones.iter().enumerate() {
            r_grid[i] = zone.radius;
        }
        if let Some(last_zone) = world_state.zones.last() {
            r_grid[n_zones] = last_zone.radius + last_zone.dr;
        }

        let y_tide = prop_mtx(r_grid, rho, g_vec, shearmod, 0);

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

        let last_radius = world_state.zones.last().map(|z| z.radius).unwrap_or(0.0);
        for (idx, zone) in world_state.zones.iter_mut().enumerate().skip(1) {
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
                    * last_radius.powi(4)
                    * (eterm + world_state.obl.sin() / 7.)
                    / r_out.powi(2)
                    * h_mu
                    * shearmod[idx].im
            } else {
                0.
            };
            world_state.w_tide_tot += w_tide;
            zone.tide_heat_rate = w_tide / 1.0e7;
        }
    }
    pub fn fluid_tide(&self, world_state: &mut WorldState, qth: &mut [f64]) {
        world_state.cesq = 0.0;
        world_state.til_t = 0.0;
        world_state.w_fluidtide_tot = 0.0;

        let n_zones = world_state.zones.len();
        if n_zones == 0 {
            return;
        }

        // Find ircore: highest index zone containing rock
        let mut ircore = 0;
        for (ir, zone) in world_state.zones.iter().enumerate().rev() {
            if zone.mass_rock > 0.0 {
                ircore = ir;
                break;
            }
        }

        // Find ir_ocean: top of ocean (liquid water > 0 and mass_rock <= 0)
        let mut ir_ocean = ircore;
        for (ir, zone) in world_state.zones.iter().enumerate().rev() {
            if zone.mass_water > 0.0 && zone.mass_rock <= 0.0 {
                ir_ocean = ir;
                break;
            }
        }

        let m_prim = self.primary_world.mass;
        let eorb = world_state.e_orb;
        let obl = world_state.obl;

        if m_prim > 0.0 && (eorb > 0.0 || obl > 0.0) && ir_ocean > ircore {
            let ocean_mid = (ir_ocean + ircore) / 2;
            let mut r_ocean = world_state.zones[ocean_mid].radius;
            if r_ocean <= 0.0 {
                r_ocean = f64::EPSILON;
            }

            let h_ocean = world_state.zones[ir_ocean].radius - world_state.zones[ircore].radius;

            // Accumulated mass up to ocean_mid
            let m_acc: f64 = world_state.zones[..=ocean_mid]
                .iter()
                .map(|z| z.mass_total)
                .sum();

            let g_ocean = GCGS * m_acc / (r_ocean * r_ocean);

            let omega_tide = world_state.n_orb;
            if omega_tide <= 0.0 {
                return;
            }

            let cesq = (g_ocean * h_ocean) / (2.0 * omega_tide * r_ocean).powi(2);

            let last_radius = world_state.zones.last().map(|z| z.radius).unwrap_or(1.0);
            let til_t_scale = 1.0;
            let ocean_top_r = world_state.zones[ir_ocean].radius;
            let r_diff = last_radius - ocean_top_r;
            let til_t = if r_diff.abs() > 1.0e-12 {
                til_t_scale * last_radius / r_diff
            } else {
                0.0
            };

            world_state.cesq = cesq;
            world_state.til_t = til_t;

            let tilom = Complex64::new(0.5, 0.0);
            let mut w_fluidtide = [0.0; 5];

            if eorb > f64::EPSILON {
                let p0 = self.tropf(cesq, til_t, tilom, 0, 0.5);
                w_fluidtide[0] = p0.re;

                let p1 = self.tropf(cesq, til_t, -tilom, 2, 3.0);
                w_fluidtide[1] = p1.re;

                let p2 = self.tropf(cesq, til_t, tilom, 2, 3.0);
                w_fluidtide[2] = p2.re;
            }

            if obl > f64::EPSILON {
                let p3 = self.tropf(cesq, til_t, -tilom, 1, 1.5);
                w_fluidtide[3] = p3.re;

                let p4 = self.tropf(cesq, til_t, tilom, 1, 1.5);
                w_fluidtide[4] = p4.re;
            }

            let mut sum_w_fluidtide = (-1.5 * eorb * 0.5) * (-1.5 * eorb * 0.5) * w_fluidtide[0]
                + (-1.0 / 8.0 * eorb * 3.0) * (-1.0 / 8.0 * eorb * 3.0) * w_fluidtide[1]
                + (7.0 / 8.0 * eorb * 3.0) * (7.0 / 8.0 * eorb * 3.0) * w_fluidtide[2]
                + (0.5 * obl * 1.5) * (0.5 * obl * 1.5) * w_fluidtide[3]
                + (0.5 * obl * 1.5) * (0.5 * obl * 1.5) * w_fluidtide[4];

            let a_orb = world_state.a_orb;
            let factor = (GCGS * m_prim / a_orb / r_ocean)
                * (last_radius / a_orb)
                * (last_radius / a_orb)
                * (r_ocean / last_radius)
                * (r_ocean / last_radius);

            sum_w_fluidtide = sum_w_fluidtide * RHO_H2OL_TH * factor.powi(2) / (2.0 * omega_tide);

            let mut w_fluidtide_tot = 0.0;
            for ir in ircore..ir_ocean {
                let vol = world_state.zones[ir].volumes().0;
                let heating = sum_w_fluidtide * vol;
                if ir < qth.len() {
                    qth[ir] += heating;
                }
                w_fluidtide_tot += heating;
            }

            world_state.w_fluidtide_tot = w_fluidtide_tot;
            world_state.heat_fluidtide += w_fluidtide_tot;
        }
    }
}
