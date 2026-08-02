//! This module solves 1D thermal evolution, tidal dissipation, and heat conduction equations for planetary bodies.

use std::{
    f64::consts::{FRAC_PI_2, PI},
    fs,
    process::exit,
};

use crate::{
    Args, FloatMat,
    consts::{
        CM, E_YOUNG_OLIV, E_YOUNG_SERP, GCGS, GRAM, KM, KM2CM, NU_POISSON_OLIV, NU_POISSON_SERP,
        PA2BA, PI_GREEK,
    },
    crack::{creep, read_data},
    input::{ChondriteType, Fracs, IcyDwarfInput},
    planet_system::{RHO_ADHS_TH, RHO_H2OL_TH, RHO_H2OS_TH, RHO_NH3L_TH, WorldState, ZoneState},
    to_faer_mat,
    traits::float_traits::FloatExt,
};
use faer::{Mat, linalg::solvers::DenseSolveCore, prelude::Solve};
use itertools::Itertools;
use num::{
    complex::{Complex64, ComplexFloat},
    pow::Pow,
    traits::Inv,
};

const K: f64 = 200.0e9 / GRAM * CM;
const C0: Complex64 = Complex64::ZERO;

impl IcyDwarfInput {
    /// Execute 1D thermal evolution and heat transport calculations for all worlds.
    ///
    /// # Parameters
    /// - `world_states`: The mutable array of world states for all moons.
    /// - `dtime`: The duration of the time step in seconds.
    /// - `real_time`: The current simulation time in seconds.
    /// - `args`: The [`Args`] containing the data path.
    pub fn thermal(
        &self,
        world_states: &mut [WorldState],
        dtime: f64,
        real_time: f64,
        args: &Args,
    ) {
        for world_state in world_states.iter_mut() {
            // 1. Calculate Pressure
            self.calculate_pressure(world_state);

            // 2. Update Porosity & Radii (Creep & Compaction)
            self.update_porosity(world_state, dtime);

            // 3. Cracking
            for zone in world_state.zones.iter_mut() {
                self.crack(zone, dtime, &read_data(args).unwrap(), false);
            }

            // 4. Tidal heating (solid)
            self.tide(world_state);
            world_state.heat_tide += world_state.w_tide_tot;

            let n_zones = world_state.zones.len();

            // 5. Radioactive Decay Heating (Qth)
            let mut qth = vec![0.0; n_zones];
            let frac_k_leached = 0.0;

            let mut q_kleached_tot = 0.0;
            let mut m_liq_tot = 0.0;
            let mut irh2os = n_zones.saturating_sub(1);

            for (i, zone) in world_state.zones.iter().enumerate() {
                let (q_zone, q_k) = zone.decay(real_time, frac_k_leached, &self.world_spec);
                qth[i] += q_zone;
                q_kleached_tot += q_k;
                m_liq_tot += zone.mass_water;
                if zone.mass_ice > 0.0 && i < irh2os {
                    irh2os = i;
                }
            }

            // Distribute heat from leached radionuclides among layers containing liquid water
            if m_liq_tot > 0.0 {
                for (i, zone) in world_state.zones.iter().enumerate() {
                    qth[i] += q_kleached_tot * zone.mass_water / m_liq_tot;
                }
            } else if n_zones > 0 {
                qth[irh2os] += q_kleached_tot;
            }

            // Track total radioactive heat
            world_state.heat_radio += qth.iter().sum::<f64>();

            // Add tidal heating to Qth
            for (i, zone) in world_state.zones.iter().enumerate() {
                qth[i] += zone.tide_heat_rate * 1.0e7;
            }

            // Add fluid tidal heating to Qth (using tropf)
            self.fluid_tide(world_state, &mut qth);

            // 6. Thermal Conductivity
            for zone in world_state.zones.iter_mut() {
                zone.kappa = zone.kapcond(&self.world_spec);
            }

            // TODO: Convection routines (convect) when ported

            // 7. Conductive Fluxes (rrflux)
            let mut rrflux = vec![0.0; n_zones + 1];
            for i in 1..n_zones {
                let r_i = world_state.zones[i - 1].radius;
                let r_next = world_state.zones[i].radius;
                let r_prev = if i > 1 {
                    world_state.zones[i - 2].radius
                } else {
                    0.0
                };
                let dr_denom = r_next - r_prev;
                if dr_denom > 0.0 {
                    rrflux[i] = -(r_i * r_i)
                        * (world_state.zones[i].kappa + world_state.zones[i - 1].kappa)
                        * (world_state.zones[i].temp - world_state.zones[i - 1].temp)
                        / dr_denom;
                }
            }

            // 8. Memorize Old Temperatures
            for zone in world_state.zones.iter_mut() {
                zone.temp_old = zone.temp;
            }

            // 9. Solve Heat Equation (Energy Update)
            for i in 0..n_zones.saturating_sub(1) {
                world_state.zones[i].energy_total +=
                    dtime * qth[i] + 4.0 * PI * dtime * (rrflux[i] - rrflux[i + 1]);
            }

            // 10. Phase Equilibrium / State Update
            let x_salt = if self
                .worlds
                .iter()
                .any(|w| w.briny && w.name == world_state.name)
            {
                0.01
            } else {
                0.0
            };
            for i in 0..n_zones {
                let _ = world_state.zones[i].apply_state(x_salt);
            }

            // 11. Surface Boundary Condition
            let temp_surf = self
                .worlds
                .iter()
                .find(|w| w.name == world_state.name)
                .map(|w| w.temp_surf)
                .unwrap_or_else(|| world_state.zones.last().map(|z| z.temp).unwrap_or(70.0));

            if let Some(last_zone) = world_state.zones.last_mut() {
                let e_rock = last_zone.mass_rock * crate::planet_system::heat_rock(temp_surf);
                let e_h2os = last_zone.mass_ice * crate::consts::QH2O * temp_surf.powi(2) * 0.5;
                let e_slush =
                    last_zone.mass_ammonia_solid * crate::consts::QADH * temp_surf.powi(2) * 0.5;
                last_zone.energy_total = e_rock + e_h2os + e_slush;
                last_zone.temp = temp_surf;
            }
        }
    }

    pub fn read_thermal_out(&self, path: &str) -> Option<Vec<Vec<ThermalOut>>> {
        let output_time_step = self.grid.time_step as usize;
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

    // TODO: finish this
    pub fn tide(&self, world_state: &mut WorldState) {
        const D_EPS: f64 = 2.22e-16;
        let base_vec = vec![0_f64; world_state.zones.len()];
        let mut shearmod = vec![C0; world_state.zones.len()];
        let mut rho = base_vec.clone();

        const MU_RIGID_ICE: f64 = 4.0e9 / GRAM * CM;

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

    /// Computes per-zone radiogenic heating from radioactive decay.
    ///
    /// Ports the C `decay()` function from `Thermal.h`. Calculates heating from
    /// long-lived radionuclides (²³⁵U, ²³⁸U, ²³²Th), short-lived ²⁶Al, and ⁴⁰K,
    /// scaled for hydration (hydrated rock has more mass from -OH but no extra
    /// radionuclides).
    ///
    /// Returns `(qth, q_kleached)`:
    /// - `qth`: radiogenic heating in this zone (erg/s), including non-leached ⁴⁰K
    /// - `q_kleached`: leached ⁴⁰K heating contribution from this zone (erg/s),
    ///   to be distributed among liquid-bearing zones by the caller
    pub fn decay(
        &self,
        t: f64,
        frac_k_leached: f64,
        world_spec: &crate::input::WorldSpec,
    ) -> (f64, f64) {
        use crate::consts::{GYR2SEC, MEV2ERG};

        const LN2: f64 = 0.6931;
        // Grams^-1 / # of Si atoms: 1e6 atoms * nucleon mass in grams * avg molar mass of rock
        const SI: f64 = 1.0 / (1.0e6 * 1.67e-24 * 151.0);

        // Half-lives in Gyr
        const T_HALF_U235: f64 = 0.704;
        const T_HALF_U238: f64 = 4.47;
        const T_HALF_TH232: f64 = 14.0;
        const T_HALF_K40: f64 = 1.265;
        const T_HALF_AL26: f64 = 0.000716;

        // Decay energies in MeV (parent-daughter mass deficit minus ~1 MeV per neutrino)
        const DE_U235: f64 = 46.74 - 4.0;
        const DE_U238: f64 = 52.07 - 6.0;
        const DE_TH232: f64 = 42.96 - 4.0;
        const DE_K40: f64 = 0.6087;
        const DE_AL26: f64 = 3.117;

        // Initial abundances per 1e6 Si atoms
        // CI (default)
        const CI_U235: f64 = 0.00592;
        const CI_U238: f64 = 0.01871;
        const CI_TH232: f64 = 0.04399;
        const CI_K40: f64 = 5.244;
        // CO
        const CO_U235: f64 = 0.00619;
        const CO_U238: f64 = 0.01942;
        const CO_TH232: f64 = 0.04293;
        const CO_K40: f64 = 2.219;

        // 26Al initial abundance: (26Al/27Al) * Al per 1e6 Si = 5e-5 * 8.41e4
        const AL26_INIT: f64 = 5.0e-5 * 8.41e4;

        let exp_decay = |t_half: f64| (-t * LN2 / (t_half * GYR2SEC)).exp();

        // Long-lived radionuclides: select abundances by chondrite type
        // C code: chondr==1 → CO, chondr==2 → CV, else → CI
        // Rust: chondrite bool — true → CO, false → CI (CV not in Rust enum)
        let (n_u235, n_u238, n_th232, n_k40) = if world_spec.chondrite == ChondriteType::CO {
            (CO_U235, CO_U238, CO_TH232, CO_K40)
        } else {
            (CI_U235, CI_U238, CI_TH232, CI_K40)
        };

        let mut s = n_u235 * DE_U235 / T_HALF_U235 * exp_decay(T_HALF_U235)
            + n_u238 * DE_U238 / T_HALF_U238 * exp_decay(T_HALF_U238)
            + n_th232 * DE_TH232 / T_HALF_TH232 * exp_decay(T_HALF_TH232);

        // Short-lived 26Al
        s += AL26_INIT * DE_AL26 / T_HALF_AL26 * exp_decay(T_HALF_AL26);

        // Convert to erg/s/g
        s *= SI * MEV2ERG / GYR2SEC * LN2;
        let s_k =
            n_k40 * DE_K40 / T_HALF_K40 * exp_decay(T_HALF_K40) * SI * MEV2ERG / GYR2SEC * LN2;

        // Per-zone radiogenic heating, scaled for hydration:
        // Hydrated rock has extra mass from -OH groups but no extra radionuclides,
        // so we subtract the water mass bound in the hydrated mineral structure.
        let rho_h2ol_th = crate::planet_system::RHO_H2OL_TH;
        let rho_rock_th = world_spec.rho_rock_th();
        let rho_hydr_th = world_spec.rho_hydr_th();

        let m_rock_dry = self.mass_rock
            - rho_h2ol_th
                * (self.mass_rock
                    / (self.x_hydr * rho_hydr_th + (1.0 - self.x_hydr) * rho_rock_th)
                    - self.mass_rock / rho_rock_th);

        let qth = m_rock_dry * (s + (1.0 - frac_k_leached) * s_k);
        let q_kleached = m_rock_dry * frac_k_leached * s_k;

        (qth, q_kleached)
    }
}

impl ZoneState {
    /// Calculates the viscosity of a water-ammonia liquid depending on
    /// temperature and ammonia mass fraction (Kargel et al. 1991).
    /// The viscosity is returned in Pa s.
    #[allow(dead_code)]
    pub fn viscosity(&self) -> f64 {
        if self.temp <= 0.0 {
            return 0.0;
        }
        let x = if self.mass_water > 0.0 {
            self.mass_ammonia_liquid / self.mass_water
        } else {
            0.0
        };
        let (a, b) = if self.temp > 240.0 {
            let a_val = -10.8143 + 0.711062 * x - 22.4943 * x * x + 41.8343 * x.powi(3)
                - 18.5149 * x.powi(4);
            let b_val =
                1819.86 + 250.822 * x + 6505.25 * x * x - 14923.4 * x.powi(3) + 7141.46 * x.powi(4);
            (a_val, b_val)
        } else {
            let a_val = -13.8628 - 68.7617 * x + 230.083 * x * x - 249.897 * x.powi(3);
            let b_val = 2701.73 + 14973.3 * x - 46174.5 * x * x + 45967.6 * x.powi(3);
            (a_val, b_val)
        };
        (a + b / self.temp).exp()
    }

    pub fn kapcond(&self, world_spec: &crate::input::WorldSpec) -> f64 {
        use crate::consts::*;
        let kaph2os = 5.67e7 / self.temp;
        let vrock = self.mass_rock
            / (self.x_hydr * world_spec.rho_hydr_th()
                + (1.0 - self.x_hydr) * world_spec.rho_rock_th());
        let vh2os = self.mass_ice / crate::planet_system::RHO_H2OS_TH;
        let vadhs = self.mass_ammonia_solid / crate::planet_system::RHO_ADHS_TH;
        let vh2ol = self.mass_water / crate::planet_system::RHO_H2OL_TH;
        let vnh3l = self.mass_ammonia_liquid / crate::planet_system::RHO_NH3L_TH;
        let total_vol =
            4.0 / 3.0 * PI_GREEK * ((self.radius + self.dr).powi(3) - self.radius.powi(3));

        let frock = vrock / total_vol;
        let fh2os = vh2os / total_vol;
        let fadhs = vadhs / total_vol;
        let fh2ol = vh2ol / total_vol;
        let fnh3l = vnh3l / total_vol;

        let mut kap: f64;
        if frock >= 1.0 - 1.0e-5 {
            kap = self.x_hydr * KAPHYDR + (1.0 - self.x_hydr) * KAPROCK;
            kap *= ((-4.0 * self.porosity / 0.08).exp()
                + (-4.4 - 4.0 * self.porosity / 0.17).exp())
            .powf(0.25);
        } else {
            let mut kapice = 0.0;
            if fh2os + fadhs + fh2ol + fnh3l > 0.0 {
                kapice = fh2os * kaph2os.ln()
                    + fadhs * KAPADHS.ln()
                    + fh2ol * KAPH2OL.ln()
                    + fnh3l * KAPNH3L.ln();
                kapice /= fh2os + fadhs + fh2ol + fnh3l;
                kapice = kapice.exp();
            }
            let krock = self.x_hydr * KAPHYDR + (1.0 - self.x_hydr) * KAPROCK;
            let b1 = -krock * (3.0 * frock - 1.0) - kapice * (2.0 - 3.0 * frock);
            let c1 = -krock * kapice;
            kap = (-b1 + (b1 * b1 - 8.0 * c1).sqrt()) * 0.25;
            if self.porosity < 0.7 {
                kap *= (1.0 - self.porosity / 0.7).powf(4.1 * self.porosity + 0.22);
            } else {
                kap = 0.0;
            }
        }
        kap
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
    if ir_core >= nr {
        return vec![vec![C0; 6]; nr];
    }

    assert!(r.len() >= nr);
    assert!(g.len() >= nr);
    assert!(shearmod.len() >= nr);

    let mut ypropmtx = vec![[[C0; 6]; 6]; nr];
    let mut ypropinv = vec![[[C0; 6]; 6]; nr];
    let mut bpropmtx = vec![[[C0; 3]; 6]; nr];

    for ir in ir_core..nr {
        let r_val = r[ir + 1];
        let r_val_2 = r_val * r_val;
        let r_val_3 = r_val_2 * r_val;
        let r_val_4 = r_val_2 * r_val_2;
        let r_val_5 = r_val_4 * r_val;

        let rho_g_r = Complex64::from(rho[ir] * g[ir] * r_val);
        let sm = shearmod[ir];
        let four_pi_g_rho = Complex64::from(4.0 * PI_GREEK * GCGS * rho[ir]);
        let two_pi_g_rho = four_pi_g_rho * 0.5;

        let rho_g_r_over_sm = rho_g_r / sm;
        let r_over_sm = r_val / sm;

        ypropmtx[ir] = [
            [
                Complex64::from(r_val_3 / 7.0),
                Complex64::from(r_val),
                C0,
                Complex64::from((2.0 * r_val_2).inv()),
                Complex64::from(r_val_4.inv()),
                C0,
            ],
            [
                5.0 * r_val_3 / 42.0.as_cplx(),
                r_val * 0.5.as_cplx(),
                C0,
                C0,
                Complex64::from(-1.0 / (3.0 * r_val_4)),
                C0,
            ],
            [
                (rho_g_r - sm) * Complex64::from(r_val_2 / 7.0),
                rho_g_r + sm * 2.0,
                Complex64::from(-rho[ir] * r_val_2),
                (rho_g_r - sm * 6.0) * 0.5 * r_val_3,
                (rho_g_r - sm * 8.0) / r_val_5,
                -rho[ir] / r_val_3.as_cplx(),
            ],
            [
                sm * 8.0 * r_val_2 / 21.0,
                sm,
                C0,
                sm / (2.0 * r_val_3),
                sm * 8.0 / (3.0 * r_val_5),
                C0,
            ],
            [
                C0,
                C0,
                Complex64::from(-r_val_2),
                C0,
                C0,
                Complex64::from(-r_val_3.inv()),
            ],
            [
                four_pi_g_rho * Complex64::from(r_val_3 / 7.0),
                four_pi_g_rho * Complex64::from(r_val),
                Complex64::from(-5.0 * r_val),
                two_pi_g_rho / Complex64::from(r_val_2),
                four_pi_g_rho / Complex64::from(r_val_4),
                C0,
            ],
        ];

        ypropinv[ir] = [
            [
                rho_g_r_over_sm - 8.0,
                Complex64::from(16.0),
                -r_over_sm,
                r_over_sm * 2.0,
                Complex64::from(rho[ir]) * r_over_sm,
                C0,
            ],
            [
                -rho_g_r_over_sm + 6.0,
                Complex64::from(-6.0),
                r_over_sm,
                C0,
                -Complex64::from(rho[ir]) * r_over_sm,
                C0,
            ],
            [four_pi_g_rho, C0, C0, C0, C0, Complex64::from(-1.0)],
            [
                rho_g_r_over_sm + 2.0,
                Complex64::from(6.0),
                -r_over_sm,
                r_over_sm * -3.0,
                Complex64::from(rho[ir]) * r_over_sm,
                C0,
            ],
            [
                -rho_g_r_over_sm + 1.0,
                -16.0.as_cplx(),
                r_over_sm,
                r_over_sm * 5.0,
                -rho[ir] * r_over_sm,
                C0,
            ],
            [
                four_pi_g_rho * r_val,
                C0,
                C0,
                C0,
                5.0.as_cplx(),
                -r_val.as_cplx(),
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

    bpropmtx[ir_core][2][0] = 1_f64.as_cplx();
    bpropmtx[ir_core][3][1] = 1_f64.as_cplx();
    bpropmtx[ir_core][5][2] = 1_f64.as_cplx();

    for ir in (ir_core + 1)..nr {
        let mut btemp = [[C0; 3]; 6];
        for i in 0..6 {
            for j in 0..3 {
                let mut sum = C0;
                for k in 0..6 {
                    sum += ypropinv[ir - 1][i][k] * bpropmtx[ir - 1][k][j];
                }
                btemp[i][j] = sum;
            }
        }
        for i in 0..6 {
            for j in 0..3 {
                let mut sum = C0;
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
    let bsurf = vec![C0, C0, Complex64::from(-5.0 / r[nr - 1])];

    let mut ytide = vec![vec![C0; 6]; nr];

    if let Some((x, _)) = gauss_jordan(&mbc, &bsurf) {
        for ir in 0..nr {
            for i in 0..6 {
                ytide[ir][i] = (0_usize..3).map(|j| bpropmtx[ir][i][j] * x[(j, 0)]).sum();
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
#[allow(dead_code)]
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
    if a.determinant() == C0 {
        // matrix is not invertible
        return None;
    }
    let b = Mat::from_fn(b.len(), 1, |i, _| b[i]);
    let lu = a.full_piv_lu();
    let x = lu.solve(b);
    let a_inv = lu.inverse();
    Some((x, a_inv))
}

/// This struct stores a snapshot of radial zone state loaded from output CSV files.
#[allow(dead_code)]
#[derive(Clone, Debug, Default)]
pub struct ThermalOut {
    /// Zone radius in centimeters.
    pub radius_km: f64,
    /// Temperature in Kelvin.
    pub temp_kelvin: f64,
    /// Rock mass in grams.
    pub mass_rock: f64,
    /// Water ice mass in grams.
    pub mass_ice: f64,
    /// Solid ammonia dihydrate mass in grams.
    pub mass_ammonia_solid: f64,
    /// Liquid water mass in grams.
    pub mass_water: f64,
    /// Liquid ammonia solution mass in grams.
    pub mass_ammonia_liquid: f64,
    /// Nusselt convection number.
    pub nusselt_num: f64,
    /// Amorphous ice fraction.
    pub ice_frac_amorphous: f64,
    /// Thermal conductivity.
    pub thermal_cond: f64,
    /// Degree of rock hydration.
    pub deg_of_hydr: f64,
    /// Matrix porosity fraction.
    pub porosity: f64,
    /// Core cracking flag.
    pub crack: bool,
    /// Tidal heating power rate in erg per second.
    pub tidal_heating_rate: f64,
}

impl ThermalOut {
    /// Parse a line of output text into a [`ThermalOut`] struct.
    ///
    /// # Parameters
    /// - `ln`: The line string to parse.
    ///
    /// # Returns
    /// An optional [`ThermalOut`] struct containing the parsed fields.
    pub fn from_line(ln: &str) -> Option<Self> {
        let parts = ln.split_whitespace().collect::<Vec<_>>();
        let radius_km = parts[0].parse::<f64>().ok()? * KM2CM;
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

    /// Calculate total zone mass in grams.
    ///
    /// # Returns
    /// The total mass sum.
    pub fn mass_total(&self) -> f64 {
        self.mass_rock
            + self.mass_ice
            + self.mass_ammonia_solid
            + self.mass_ammonia_liquid
            + self.mass_water
    }

    /// Calculate total zone volume and phase volumes.
    ///
    /// # Parameters
    /// - `input`: The [`IcyDwarfInput`] configuration reference.
    ///
    /// # Returns
    /// A tuple containing total volume and a [`Fracs`] struct of phase volumes.
    #[allow(dead_code)]
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

    /// Calculate phase mass fractions for the output zone.
    ///
    /// # Returns
    /// A [`Fracs`] struct containing phase mass fractions.
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

    #[test]
    fn test_viscosity() {
        let zone = ZoneState {
            temp: 273.15,
            mass_water: 1000.0,
            mass_ammonia_liquid: 100.0,
            ..Default::default()
        };
        let visc = zone.viscosity();
        assert!(visc > 0.0 && visc.is_finite());
    }
}
