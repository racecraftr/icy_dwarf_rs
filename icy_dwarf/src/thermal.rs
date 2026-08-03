//! This module solves 1D thermal evolution, tidal dissipation, and heat conduction equations for planetary bodies.

use std::{
    f64::consts::{FRAC_PI_3, PI},
    fs,
    process::exit,
};

use crate::{
    Args, FloatMat,
    consts::*,
    crack::{creep, read_data},
    input::{ChondriteType, Fracs, IcyDwarfInput},
    planet_system::{RHO_ADHS_TH, RHO_H2OL_TH, RHO_H2OS_TH, RHO_NH3L_TH, WorldState, ZoneState},
    to_faer_mat,
    traits::float_traits::FloatExt,
};
use faer::{Mat, linalg::solvers::DenseSolveCore, prelude::Solve};
use itertools::Itertools;
use num::{complex::Complex64, traits::Inv};

pub mod thermal_out;
use thermal_out::ThermalOut;
mod tide;

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

            // Rock Hydration & Dehydration
            self.process_hydration(world_state, &mut qth, dtime);

            // Layer Differentiation (melting / Rayleigh-Taylor)
            let world_info = self.worlds.iter().find(|w| w.name == world_state.name);
            let x_p = world_info.map(|w| w.ammonia).unwrap_or(0.0);
            let t_liq = if x_p >= 1.0e-2 { 174.0 } else { 271.0 };
            let irdiff_old = world_state.irdiff;

            for i in 0..n_zones.saturating_sub(1) {
                if i > world_state.irdiff && world_state.zones[i].temp > t_liq {
                    world_state.irdiff = i;
                }
            }
            if world_state.irdiff > n_zones / 2 {
                for i in 0..n_zones.saturating_sub(1) {
                    if i > world_state.irdiff && world_state.zones[i].temp > TDIFF {
                        world_state.irdiff = i;
                    }
                }
            }
            if world_state.irdiff > 0 && world_state.irdiff != irdiff_old {
                let irdiff = world_state.irdiff;
                self.separate(world_state, irdiff);
            }

            // Gravitational Heating Release
            self.gravitational_heating(world_state, &mut qth, dtime);

            // Parameterized Convection
            self.convect(world_state);

            // 6. Thermal Conductivity
            for zone in world_state.zones.iter_mut() {
                zone.kappa = zone.kapcond(&self.world_spec) * zone.nusselt;
            }

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
                let e_h2os = last_zone.mass_ice * QH2O * temp_surf.powi(2) * 0.5;
                let e_slush = last_zone.mass_ammonia_solid * QADH * temp_surf.powi(2) * 0.5;
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

    /// Perform core separation / unmixing of differentiated layers.
    pub fn separate(&self, world_state: &mut WorldState, irdiff: usize) {
        let n_zones = world_state.zones.len();
        if n_zones == 0 || irdiff >= n_zones {
            return;
        }

        let world_info = self.worlds.iter().find(|w| w.name == world_state.name);
        let x_fines = world_info.map(|w| w.rock_frac).unwrap_or(0.0);
        let x_pores = world_info.map(|w| w.rock_h20).unwrap_or(0.0);

        let mut m_rock_new = vec![0.0; n_zones];
        let mut m_h2os_new = vec![0.0; n_zones];
        let mut m_adhs_new = vec![0.0; n_zones];
        let mut m_h2ol_new = vec![0.0; n_zones];
        let mut m_nh3l_new = vec![0.0; n_zones];

        let mut v_rock_new = vec![0.0; n_zones];
        let mut v_h2os_new = vec![0.0; n_zones];
        let mut v_adhs_new = vec![0.0; n_zones];
        let mut v_h2ol_new = vec![0.0; n_zones];
        let mut v_nh3l_new = vec![0.0; n_zones];

        let mut vol_cell = vec![0.0; n_zones];
        for jr in 0..=irdiff {
            vol_cell[jr] = world_state.zones[jr].volumes().0;
        }

        let mut jr = 0;
        let mut ircore = jr;

        for ir in 0..=irdiff {
            let v_rock_ir = world_state.zones[ir].volumes().1.0;
            let m_rock_ir = world_state.zones[ir].mass_rock;

            if v_rock_new[jr] >= vol_cell[jr] * (1.0 - x_pores) && v_rock_ir > 0.0 {
                let q = (v_rock_new[jr] - vol_cell[jr] * (1.0 - x_pores)) / v_rock_ir;
                v_rock_new[jr] = vol_cell[jr] * (1.0 - x_pores);
                m_rock_new[jr] -= q * m_rock_ir;
                vol_cell[jr] *= x_pores;
                jr = (jr + 1).min(n_zones - 1);
                ircore = jr;
                v_rock_new[jr] = q * v_rock_ir;
                m_rock_new[jr] = q * m_rock_ir;
            }

            v_rock_new[jr] += v_rock_ir * (1.0 - x_fines);
            m_rock_new[jr] += m_rock_ir * (1.0 - x_fines);

            if v_rock_new[jr] >= vol_cell[jr] * (1.0 - x_pores) && v_rock_ir > 0.0 {
                let q = (v_rock_new[jr] - vol_cell[jr] * (1.0 - x_pores)) / v_rock_ir;
                v_rock_new[jr] = vol_cell[jr] * (1.0 - x_pores);
                m_rock_new[jr] -= q * m_rock_ir;
                vol_cell[jr] *= x_pores;
                jr = (jr + 1).min(n_zones - 1);
                ircore = jr;
                v_rock_new[jr] = q * v_rock_ir;
                m_rock_new[jr] = q * m_rock_ir;
            }
        }
        vol_cell[ircore] = (vol_cell[ircore] - v_rock_new[ircore]).max(0.0);

        let mut v_ice = vol_cell[ircore];
        let mut m_fines = 0.0;
        let mut v_fines = 0.0;

        for ir in 0..=irdiff {
            if ir >= ircore + 1 {
                v_ice += world_state.zones[ir].volumes().0;
            }
            let v_rock_ir = world_state.zones[ir].volumes().1.0;
            let m_rock_ir = world_state.zones[ir].mass_rock;
            m_fines += m_rock_ir * x_fines;
            v_fines += v_rock_ir * x_fines;
        }

        if v_ice > 0.0 {
            m_rock_new[ircore] += m_fines * vol_cell[ircore] / v_ice;
            v_rock_new[ircore] += v_fines * vol_cell[ircore] / v_ice;
            vol_cell[ircore] = (vol_cell[ircore] - v_fines * vol_cell[ircore] / v_ice).max(0.0);

            for ir in (ircore + 1)..=irdiff {
                let cell_v = world_state.zones[ir].volumes().0;
                m_rock_new[ir] += m_fines * cell_v / v_ice;
                v_rock_new[ir] += v_fines * cell_v / v_ice;
                vol_cell[ir] = (vol_cell[ir] - v_fines * cell_v / v_ice).max(0.0);
            }
        }

        if x_pores > 0.0 {
            jr = 0;
        }

        for ir in 0..=irdiff {
            let zone = &world_state.zones[ir];
            let vol_slush_ir = zone.volumes().1.2 + zone.volumes().1.3 + zone.volumes().1.4;

            let volume1 = v_adhs_new[jr] + v_h2ol_new[jr] + v_nh3l_new[jr];
            if volume1 >= vol_cell[jr] && vol_slush_ir > 0.0 {
                let q = (volume1 - vol_cell[jr]) / vol_slush_ir;
                v_h2ol_new[jr] -= q * zone.volumes().1.3;
                v_nh3l_new[jr] -= q * zone.volumes().1.4;
                v_adhs_new[jr] -= q * zone.volumes().1.2;
                m_h2ol_new[jr] -= q * zone.mass_water;
                m_nh3l_new[jr] -= q * zone.mass_ammonia_liquid;
                m_adhs_new[jr] -= q * zone.mass_ammonia_solid;
                vol_cell[jr] = 0.0;
                jr = (jr + 1).min(n_zones - 1);
                v_adhs_new[jr] = q * zone.volumes().1.2;
                v_h2ol_new[jr] = q * zone.volumes().1.3;
                v_nh3l_new[jr] = q * zone.volumes().1.4;
                m_adhs_new[jr] = q * zone.mass_ammonia_solid;
                m_h2ol_new[jr] = q * zone.mass_water;
                m_nh3l_new[jr] = q * zone.mass_ammonia_liquid;
            }

            v_adhs_new[jr] += zone.volumes().1.2;
            v_h2ol_new[jr] += zone.volumes().1.3;
            v_nh3l_new[jr] += zone.volumes().1.4;
            m_adhs_new[jr] += zone.mass_ammonia_solid;
            m_h2ol_new[jr] += zone.mass_water;
            m_nh3l_new[jr] += zone.mass_ammonia_liquid;
        }

        let irice = jr;
        vol_cell[irice] =
            (vol_cell[irice] - v_adhs_new[irice] - v_h2ol_new[irice] - v_nh3l_new[irice]).max(0.0);

        for ir in 0..=irdiff {
            let zone = &world_state.zones[ir];
            let v_h2os_ir = zone.volumes().1.1;

            if v_h2os_new[jr] >= vol_cell[jr] && v_h2os_ir > 0.0 {
                let q = (v_h2os_new[jr] - vol_cell[jr]) / v_h2os_ir;
                v_h2os_new[jr] -= q * v_h2os_ir;
                m_h2os_new[jr] -= q * zone.mass_ice;
                vol_cell[jr] = 0.0;
                jr = (jr + 1).min(n_zones - 1);
                v_h2os_new[jr] = q * v_h2os_ir;
                m_h2os_new[jr] = q * zone.mass_ice;
            }

            v_h2os_new[jr] += v_h2os_ir;
            m_h2os_new[jr] += zone.mass_ice;
        }

        for ir in 0..=irdiff {
            let zone = &mut world_state.zones[ir];
            zone.mass_rock = m_rock_new[ir];
            zone.mass_ice = m_h2os_new[ir];
            zone.mass_ammonia_solid = m_adhs_new[ir];
            zone.mass_water = m_h2ol_new[ir];
            zone.mass_ammonia_liquid = m_nh3l_new[ir];
            zone.mass_total = zone.mass_rock
                + zone.mass_ice
                + zone.mass_ammonia_solid
                + zone.mass_water
                + zone.mass_ammonia_liquid;
        }
    }

    /// Process rock hydration and dehydration across grid zones and update heats of reaction.
    pub fn process_hydration(&self, world_state: &mut WorldState, qth: &mut [f64], dtime: f64) {
        let n_zones = world_state.zones.len();
        if n_zones == 0 || dtime <= 0.0 {
            return;
        }

        const F_MEM: f64 = 0.9;
        let rho_rock_th = self.world_spec.rho_rock_th();
        let rho_hydr_th = self.world_spec.rho_hydr_th();
        let rho_h2ol_th = RHO_H2OL_TH;

        let mut dont_dehydrate = vec![false; n_zones];

        let mut ircore = 0;
        for (ir, zone) in world_state.zones.iter().enumerate() {
            if zone.mass_rock > 0.0 {
                ircore = ir;
            }
        }

        let mut ircrack = n_zones;
        for ir in (0..=ircore).rev() {
            if world_state.zones[ir].crack > 0.0 {
                ircrack = ir;
            } else {
                break;
            }
        }

        if ircore > 0 && ircrack < ircore {
            for ir in (ircrack..ircore).rev() {
                let temp = world_state.zones[ir].temp;
                let x_hydr_old = world_state.zones[ir].x_hydr;

                if temp < crate::consts::TDEHYDR_MAX && x_hydr_old <= 0.99 {
                    let target_x = if temp < crate::consts::TDEHYDR_MIN {
                        1.0
                    } else {
                        1.0 - (temp - crate::consts::TDEHYDR_MIN)
                            / (crate::consts::TDEHYDR_MAX - crate::consts::TDEHYDR_MIN)
                    };
                    let new_x =
                        (F_MEM * x_hydr_old + (1.0 - F_MEM) * target_x).clamp(x_hydr_old, 1.0);

                    if new_x > x_hydr_old {
                        world_state.zones[ir].x_hydr = new_x;
                        dont_dehydrate[ir] = true;

                        let v_moved = (world_state.zones[ir].mass_rock
                            / (new_x * rho_hydr_th + (1.0 - new_x) * rho_rock_th)
                            - world_state.zones[ir].mass_rock
                                / (x_hydr_old * rho_hydr_th + (1.0 - x_hydr_old) * rho_rock_th))
                            .min(world_state.zones[ir].mass_water / rho_h2ol_th)
                            .max(0.0);

                        let m_water_trans = v_moved * rho_h2ol_th;
                        world_state.zones[ir].mass_water =
                            (world_state.zones[ir].mass_water - m_water_trans).max(0.0);
                        world_state.zones[ir].mass_rock += m_water_trans;
                    }
                }
            }
        }

        for ir in 0..=ircore {
            if ir >= n_zones {
                break;
            }
            if dont_dehydrate[ir] {
                continue;
            }

            let temp = world_state.zones[ir].temp;
            let x_hydr_old = world_state.zones[ir].x_hydr;

            if temp > TDEHYDR_MIN && x_hydr_old >= 0.01 {
                let target_x = if temp >= crate::consts::TDEHYDR_MAX {
                    0.0
                } else {
                    1.0 - (temp - crate::consts::TDEHYDR_MIN)
                        / (crate::consts::TDEHYDR_MAX - crate::consts::TDEHYDR_MIN)
                };
                let new_x = (F_MEM * x_hydr_old + (1.0 - F_MEM) * target_x).clamp(0.0, x_hydr_old);

                if new_x < x_hydr_old {
                    world_state.zones[ir].x_hydr = new_x;

                    let d_vol = world_state.zones[ir].volumes().0;
                    let v_rock_new = world_state.zones[ir].mass_rock
                        / (new_x * rho_hydr_th + (1.0 - new_x) * rho_rock_th);
                    let v_h2ol_new = (d_vol - v_rock_new).max(0.0);

                    let m_h2ol_new = v_h2ol_new * rho_h2ol_th;
                    world_state.zones[ir].mass_water += m_h2ol_new;
                    world_state.zones[ir].mass_rock = (world_state.zones[ir].mass_total
                        - world_state.zones[ir].mass_water)
                        .max(0.0);
                }
            }
        }

        for (ir, zone) in world_state.zones.iter_mut().enumerate() {
            let delta_x = zone.x_hydr - zone.x_hydr_old;
            if delta_x.abs() > 1.0e-10 {
                let q_rxn = delta_x * zone.mass_rock * crate::consts::HHYDR / dtime;
                if ir < qth.len() {
                    qth[ir] += q_rxn;
                }
                if delta_x > 0.0 {
                    world_state.heat_serp += q_rxn;
                } else {
                    world_state.heat_dehydr += -q_rxn;
                }
            }
            zone.x_hydr_old = zone.x_hydr;
        }
    }

    /// Calculate gravitational potential energy release during differentiation.
    pub fn gravitational_heating(&self, world_state: &mut WorldState, qth: &mut [f64], dtime: f64) {
        let n_zones = world_state.zones.len();
        if n_zones == 0 || dtime <= 0.0 {
            return;
        }

        let phi_old = world_state.phi;
        let mut phi =
            0.6 * GCGS * world_state.zones[0].mass_total.powi(2) / world_state.zones[0].radius;

        let mut m_acc = world_state.zones[0].mass_total;
        for ir in 1..n_zones {
            let r_prev = world_state.zones[ir - 1].radius;
            let r_curr = world_state.zones[ir].radius;
            let r_avg = (r_prev + r_curr) * 0.5;
            phi += GCGS * m_acc * world_state.zones[ir].mass_total / r_avg;
            m_acc += world_state.zones[ir].mass_total;
        }

        if (phi - phi_old).abs() < 1.0e-5 * phi_old {
            phi = phi_old;
        }
        world_state.phi = phi;

        let irdiff = world_state.irdiff;
        if irdiff > 0 && irdiff < n_zones && phi > phi_old {
            let total_diff_vol: f64 = world_state.zones[..=irdiff]
                .iter()
                .map(|z| z.volumes().0)
                .sum();

            if total_diff_vol > 0.0 {
                let delta_phi = (phi - phi_old) / dtime;
                for ir in 0..=irdiff {
                    let zone_vol = world_state.zones[ir].volumes().0;
                    let heat_rate = delta_phi * (zone_vol / total_diff_vol);
                    if ir < qth.len() {
                        qth[ir] += heat_rate;
                    }
                    world_state.heat_grav += heat_rate;
                }
            }
        }
    }

    /// Calculate convective heat transport (hydrothermal, ocean mud, subsolidus ice shell).
    pub fn convect(&self, world_state: &mut WorldState) {
        let n_zones = world_state.zones.len();
        if n_zones == 0 {
            return;
        }

        for zone in world_state.zones.iter_mut() {
            zone.nusselt = 1.0;
        }

        let mut ircore = 0;
        for (ir, zone) in world_state.zones.iter().enumerate() {
            if zone.mass_rock > 0.0 {
                ircore = ir;
            }
        }
        let irdiff = world_state.irdiff.min(n_zones - 1);

        if irdiff > ircore + 2 {
            let r_in = world_state.zones[ircore].radius;
            let r_out = world_state.zones[irdiff].radius;
            let d_shell = (r_out - r_in).max(1.0e-5);

            let t_in = world_state.zones[ircore].temp;
            let t_out = world_state.zones[irdiff].temp;
            let delta_t = (t_in - t_out).max(0.0);

            if delta_t > 0.0 {
                let g = GCGS * world_state.zones[irdiff].mass_total / (r_out * r_out);
                let alpha_ice = 1.0e-4;
                let kappa_ice = 1.0e-2;
                let eta_ice = world_state.zones[ircore].creep().max(1.0e13);

                let ra = (RHO_H2OS_TH * g * alpha_ice * delta_t * d_shell.powi(3))
                    / (kappa_ice * eta_ice);

                let ra_crit = 1000.0;
                if ra > ra_crit {
                    let nu = (ra / ra_crit).powf(0.3).clamp(1.0, 100.0);
                    for ir in ircore..=irdiff {
                        world_state.zones[ir].nusselt = nu;
                    }
                }
            }
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
        let total_vol = 4.0 * FRAC_PI_3 * ((self.radius + self.dr).powi(3) - self.radius.powi(3));

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
        let four_pi_g_rho = Complex64::from(4.0 * PI * GCGS * rho[ir]);
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
