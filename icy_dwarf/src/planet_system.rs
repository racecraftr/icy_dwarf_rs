//! This module defines interior zone states, planetary world structures, and main simulation progression functions.

use std::f64::consts::{FRAC_PI_3, PI};

use itertools::Itertools;

use crate::{
    GLOBAL_ARGS, append_output,
    consts::*,
    crack, create_output,
    input::{Fracs, IcyDwarfInput, IcyWorld, QMode, TidalQ, WorldSpec},
    traits::float_traits::FloatExt,
};

const TODAY: f64 = 4568.2 * MYR2SEC;

/// The density of water ice in grams per cubic centimeter converted to mass units.
pub const RHO_H2OS_TH: f64 = RHO_H2OS * GRAM;
/// The density of liquid water in grams per cubic centimeter converted to mass units.
pub const RHO_H2OL_TH: f64 = RHO_H2OS_TH;
/// The density of solid ammonia dihydrate converted to mass units.
pub const RHO_ADHS_TH: f64 = RHO_ADHS * GRAM;
/// The density of liquid ammonia solution converted to mass units.
pub const RHO_NH3L_TH: f64 = XC / (1.0 / RHO_H2OL_TH) + (1.0 / RHO_ADHS_TH - 1.0 / RHO_H2OS_TH);

/// Represents the state of a radial layer of an [`IcyWorld`].
#[derive(Clone, Debug, Default)]
pub struct ZoneState {
    /// Inner radius of shell in centimeters.
    pub radius: f64,

    /// Radial distance from inner radius to outer radius in centimeters.
    pub dr: f64,
    /// Temperature of the zone in Kelvin.
    pub temp: f64,
    /// Temperature of the zone at previous step in Kelvin.
    pub temp_old: f64,
    /// Total mass of the zone in grams.
    pub mass_total: f64,
    /// Mass of rock in the zone in grams.
    pub mass_rock: f64,
    /// Initial mass of rock in the zone in grams.
    pub mass_rock_init: f64,
    /// Mass of water ice in the zone in grams.
    pub mass_ice: f64,
    /// Mass of solid ammonia dihydrate in grams.
    pub mass_ammonia_solid: f64,
    /// Mass of liquid water in grams.
    pub mass_water: f64,
    /// Mass of liquid ammonia solution in grams.
    pub mass_ammonia_liquid: f64,
    /// Total thermal energy in erg units.
    pub energy_total: f64,
    /// Matrix porosity fraction.
    pub porosity: f64,
    /// Hydrostatic pressure in pascals.
    pub pressure: f64,
    /// Rock hydration fraction.
    pub x_hydr: f64,
    /// Previous rock hydration fraction.
    pub x_hydr_old: f64,
    /// Thermal conductivity.
    pub kappa: f64,
    /// Nusselt convection number.
    pub nusselt: f64,
    /// Crack state identifier code.
    pub crack: f64,
    /// Crack radius in meters.
    pub crack_size: f64,
    /// Pore fluid pressure in pascals.
    pub p_pore: f64,
    /// Hydration stress in pascals.
    pub p_hydr: f64,
    /// Chemical activity quotients array.
    pub act: [f64; 3],
    /// Stress components tensor array.
    pub stress: [f64; 12],
    /// Open fracture volume fraction.
    pub frac_open: f64,
    /// Tidal heating power rate in erg per second.
    pub tide_heat_rate: f64,
}

impl ZoneState {
    /// Calculate the total radial volume and phase volumes of the zone.
    ///
    /// # Returns
    /// A tuple containing total volume and a [`Fracs`] struct of phase volumes.
    pub fn volumes(&self) -> (f64, Fracs) {
        // Volume = outer vol - inner vol
        //        = pi (r + dr)^2 - pi r^2
        //        = pi ((r + dr)^2 - r^2)
        //        = pi (r^2 + 2dr + dr^2 - r^2) = pi(2dr + dr^2)
        let total_vol = 4.0 * FRAC_PI_3 * ((self.radius + self.dr).powi(3) - self.radius.powi(3));
        let Fracs(f_rock, f_ice, f_as, f_water, f_al) = self.fracs();
        (
            total_vol,
            Fracs(
                total_vol * f_rock,
                total_vol * f_ice,
                total_vol * f_as,
                total_vol * f_water,
                total_vol * f_al,
            ),
        )
    }

    /// Calculate mass fractions for all material phases in the zone.
    ///
    /// # Returns
    /// A [`Fracs`] struct containing phase mass fractions.
    pub fn fracs(&self) -> Fracs {
        Fracs(
            self.mass_rock / self.mass_total,
            self.mass_ice / self.mass_total,
            self.mass_ammonia_solid / self.mass_total,
            self.mass_water / self.mass_total,
            self.mass_ammonia_liquid / self.mass_total,
        )
    }

    /// Calculate total internal thermal energy of the zone.
    ///
    /// # Parameters
    /// - `t_init`: The reference initial temperature in Kelvin.
    ///
    /// # Returns
    /// A tuple containing total energy and component energies for rock, ice, and slush.
    pub fn _internal_energy(&self, t_init: f64) -> (f64, (f64, f64, f64)) {
        let n = t_init.powi(2) * 0.5;
        let energy_rock = self.mass_rock * heat_rock(t_init);
        let energy_ice = self.mass_ice * QH2O * n;
        let energy_slush = self.mass_ammonia_solid * n;
        let sum = energy_rock + energy_ice + energy_slush;
        (sum, (energy_rock, energy_ice, energy_slush))
    }

    /// Calculate brittle strength of the rock matrix in the zone.
    ///
    /// # Returns
    /// The brittle strength value in pascals.
    pub fn _brittle_strength(&self) -> f64 {
        crack::strain(self.pressure, self.x_hydr, self.temp, self.porosity).0
    }

    /// Recalculate zone temperature and phase state from total internal energy.
    ///
    /// # Parameters
    /// - `x_salt`: The salinity mass fraction parameter.
    pub fn apply_state(&mut self, x_salt: f64) -> Result<(), String> {
        let specific_energy = self.energy_total / self.mass_total;
        let Fracs(frock, mut fh2os, mut fadhs, mut fh2ol, mut fnh3l) = self.fracs();

        let [mut gh2os, mut gadhs, mut gh2ol, mut gnh3l] =
            [fh2os, fadhs, fh2ol, fnh3l].map(|n| (n / (1.0 - frock)).unwrap_or_nan(0.));
        let x = gnh3l + XC * gadhs;

        let mut t_lo = 20.0;
        let mut t_hi = 5000.0;
        let mut t_md = (t_lo + t_hi) / 2.0;

        for _ in 0..30 {
            // Calculate Elo
            let tp = t_lo;
            let erock_lo = heat_rock(tp);
            let [eice_lo, ..] = heat_ice(tp, x, x_salt);
            let elo = frock * erock_lo + (1.0 - frock) * eice_lo;

            // Calculate Emd
            let tp = t_md;
            let erock_md = heat_rock(tp);
            let [eice_md, gh2os_md, gadhs_md, gh2ol_md, gnh3l_md] = heat_ice(tp, x, x_salt);
            let emd = frock * erock_md + (1.0 - frock) * eice_md;

            // Calculate Ehi
            let tp = t_hi;
            let erock_hi = heat_rock(tp);
            let [eice_hi, ..] = heat_ice(tp, x, x_salt);
            let ehi = frock * erock_hi + (1.0 - frock) * eice_hi;

            if specific_energy >= elo
                && specific_energy <= ehi
                && elo > 0.0
                && ehi > 0.0
                && emd > 0.0
            {
                if specific_energy <= emd {
                    t_hi = t_md;
                } else {
                    t_lo = t_md;
                }
                t_md = (t_lo + t_hi) / 2.0;

                // Keep updating fractions with the current midpoint
                gh2os = gh2os_md;
                gadhs = gadhs_md;
                gh2ol = gh2ol_md;
                gnh3l = gnh3l_md;
            } else {
                return Err(format!(
                    "Could not compute temperature. Tlo={}, Thi={}, Tmd={}, Elo={}, Ehi={}, Emd={}, E={}",
                    t_lo, t_hi, t_md, elo, ehi, emd, specific_energy
                ));
            }
        }

        self.temp = t_md;

        fh2os = (1.0 - frock) * gh2os;
        fadhs = (1.0 - frock) * gadhs;
        fh2ol = (1.0 - frock) * gh2ol;
        fnh3l = (1.0 - frock) * gnh3l;

        self.mass_ice = self.mass_total * fh2os;
        self.mass_ammonia_solid = self.mass_total * fadhs;
        self.mass_water = self.mass_total * fh2ol;
        self.mass_ammonia_liquid = self.mass_total * fnh3l;

        Ok(())
    }
}

/// Represents the state of an [`IcyWorld`].
#[derive(Clone, Debug)]
pub struct WorldState {
    /// The name of the satellite world.
    pub name: String,
    /// Semi-major axis of the orbit in centimeters.
    pub a_orb: f64,
    /// Orbital eccentricity.
    pub e_orb: f64,
    /// Orbital inclination in degrees.
    pub i_orb: f64,
    /// Obliquity in degrees.
    pub obl: f64,
    /// Mean motion orbital frequency in radians per second.
    pub n_orb: f64,
    /// Mean longitude in radians.
    pub lambda: f64,
    /// Longitude of pericenter in radians.
    pub omega: f64,
    /// Rotation rate in radians per second.
    pub spin: f64,
    /// Second Love number of the body.
    pub k2: f64,
    /// Tidal dissipation factor Q.
    pub q_tide: f64,
    /// Dimensionless moment of inertia coefficient.
    pub moi: f64,
    /// Vector of radial zone states in the body.
    pub zones: Vec<ZoneState>,
    /// Active accounted orbital resonances vector.
    pub res_acct_for: Vec<f64>,
    /// Previous step first Poincaré variable.
    pub h_old: f64,
    /// Previous step second Poincaré variable.
    pub k_old: f64,
    /// Previous step semi-major axis in centimeters.
    pub a_old: f64,
    /// Previous step secular eccentricity coefficient.
    pub cs_ee_old: f64,
    /// Previous step mutual eccentricity coefficient.
    pub cs_eep_old: f64,
    /// Previous step resonant inner eccentricity coefficient.
    pub cr_e_old: f64,
    /// Previous step resonant outer eccentricity coefficient.
    pub cr_ep_old: f64,
    /// Previous step resonant squared inner eccentricity coefficient.
    pub cr_ee_old: f64,
    /// Previous step resonant product eccentricity coefficient.
    pub cr_eep_old: f64,
    /// Previous step resonant squared outer eccentricity coefficient.
    pub cr_epep_old: f64,
    /// Cumulative solid body tidal work in erg units.
    pub w_tide_tot: f64,
    /// Cumulative ocean fluid tidal work in erg units.
    pub w_fluidtide_tot: f64,
    /// Radiogenic heating rate in erg per second.
    pub heat_radio: f64,
    /// Gravitational differentiation heating rate in erg per second.
    pub heat_grav: f64,
    /// Serpentinization reaction heating rate in erg per second.
    pub heat_serp: f64,
    /// Dehydration endothermic heat absorption rate in erg per second.
    pub heat_dehydr: f64,
    /// Solid body tidal heating rate in erg per second.
    pub heat_tide: f64,
    /// Fluid ocean tidal heating rate in erg per second.
    pub heat_fluidtide: f64,
    /// Equivalent gravity wave speed squared.
    pub cesq: f64,
    /// Fluid tidal dissipation time factor.
    pub til_t: f64,
    /// Total cracked rock mass in grams.
    pub m_cracked_rock: f64,
    /// Total liquid water mass in grams.
    pub m_liq: f64,
    /// Gravitational binding potential energy in erg units.
    pub phi: f64,
    /// Outermost differentiated grid zone index.
    pub irdiff: usize,
}

impl IcyDwarfInput {
    /// Execute the main simulation loop across all satellites over time.
    ///
    /// # Parameters
    /// - `output_path`: The file system directory path string for output files.
    pub fn planet_system(&self, output_path: &String) {
        let dtime = self.grid.time_step * 1.0e-6 * MYR2SEC;
        let n_time = (self.grid.time_total / self.grid.time_step) as usize;
        let _n_steps = (self.grid.output_every / self.grid.time_step) as usize;

        let t_form_min = self
            .worlds
            .iter()
            .fold(f64::MAX, |acc, w| acc.min(w.t_form));

        //-------------------------------------------------------------------
        //                Normal init (from input file only)
        //-------------------------------------------------------------------
        let mut world_states: Vec<WorldState> = self
            .worlds
            .iter()
            .enumerate()
            .map(|(idx, world)| {
                let name = &world.name;
                for s in ["Thermal", "Heats", "Crack_depth_WR", "Crack_stresses"] {
                    let _ = create_output(output_path, format!("{}_{}_{}.txt", idx, name, s));
                }
                if self.primary_world.mass > 0.0 {
                    let _ = create_output(output_path, format!("{}_{}_Orbit.txt", idx, name));
                }

                let rho_ice = 1.0
                    / ((world.ammonia / XC) / RHO_ADHS_TH
                        + (1.0 - world.ammonia / XC) / RHO_H2OS_TH);
                let rho_rock_hydr = world.hydr_init * self.world_spec.rho_hydr_th()
                    + (1.0 - world.hydr_init) * self.world_spec.rho_rock_th();
                let frock_pm =
                    (1.0 - rho_ice / world.planetary_dens) / (1.0 - rho_ice / rho_rock_hydr);

                let dr_grid = world.planetary_rad / (self.grid.n_zones as f64);
                let zones = (0..self.grid.n_zones)
                    .map(|ir| {
                        let current_r = dr_grid * ir as f64;
                        let next_r = current_r + dr_grid;
                        let d_vol = 4.0 * FRAC_PI_3 * (next_r.powi(3) - current_r.powi(3));
                        let d_m = d_vol * world.planetary_dens;
                        let m_rock = d_m * frock_pm;
                        let m_h2os = d_m * (1.0 - frock_pm) * (1.0 - world.ammonia / XC);
                        let m_adhs = d_m * (1.0 - frock_pm) * (world.ammonia / XC);

                        let e_rock = m_rock * heat_rock(world.temp_init);
                        let e_h2os = m_h2os * QH2O * world.temp_init.powi(2) * 0.5;
                        let e_slush = m_adhs * QADH * world.temp_init.powi(2) * 0.5;
                        let d_e = e_rock + e_h2os + e_slush;

                        ZoneState {
                            radius: next_r,
                            dr: dr_grid,
                            temp: world.temp_init,
                            temp_old: world.temp_init,
                            mass_total: d_m,
                            mass_rock: m_rock,
                            mass_rock_init: m_rock,
                            mass_ice: m_h2os,
                            mass_ammonia_solid: m_adhs,
                            mass_water: 0.0,
                            mass_ammonia_liquid: 0.0,
                            energy_total: d_e,
                            porosity: world.por_init,
                            pressure: 0.0,
                            x_hydr: world.hydr_init,
                            x_hydr_old: world.hydr_init,
                            kappa: 0.0,
                            nusselt: 1.0,
                            crack: 0.0,
                            crack_size: 0.0,
                            p_pore: 0.0,
                            p_hydr: 0.0,
                            act: [0.0; 3],
                            stress: [0.0; 12],
                            frac_open: 0.0,
                            tide_heat_rate: 0.0,
                        }
                    })
                    .collect_vec();

                let n_orb = (GCGS * self.primary_world.mass / world.orb_a_init.powi(3))
                    .sqrt()
                    .unwrap_or_nan(0.);

                WorldState {
                    name: world.name.clone(),
                    a_orb: world.orb_a_init,
                    e_orb: world.orb_e_init,
                    i_orb: world.orb_i_init,
                    obl: world.orb_o_init,
                    n_orb,
                    lambda: 0.32 + (idx as f64) * 0.77,
                    omega: 0.58 + (idx as f64) * 0.27,
                    spin: n_orb,
                    k2: 0.001,
                    q_tide: 1.0e10,
                    moi: 0.4,
                    zones,
                    res_acct_for: vec![0.; self.n_moons()],
                    h_old: 0.0,
                    k_old: 0.0,
                    a_old: 0.0,
                    cs_ee_old: 0.0,
                    cs_eep_old: 0.0,
                    cr_e_old: 0.0,
                    cr_ep_old: 0.0,
                    cr_ee_old: 0.0,
                    cr_eep_old: 0.0,
                    cr_epep_old: 0.0,
                    w_tide_tot: 0.0,
                    w_fluidtide_tot: 0.0,
                    heat_radio: 0.0,
                    heat_grav: 0.0,
                    heat_serp: 0.0,
                    heat_dehydr: 0.0,
                    heat_tide: 0.0,
                    heat_fluidtide: 0.0,
                    cesq: 0.0,
                    til_t: 0.0,
                    m_cracked_rock: 0.0,
                    m_liq: 0.0,
                    phi: 0.0,
                    irdiff: 0,
                }
            })
            .collect();

        // Handle initial differentiation if requested
        for (idx, world) in self.worlds.iter().enumerate() {
            if world.start_diff {
                let n_zones = world_states[idx].zones.len();
                if n_zones > 0 {
                    let max_diff = n_zones - 1;
                    world_states[idx].irdiff = max_diff;
                    self.separate(&mut world_states[idx], max_diff);
                }
            }
        }

        if self.primary_world.mass > 0.0 {
            for s in [
                "Primary",
                "Resonances",
                "ResAcctFor",
                "PCapture",
                "icydwarf_outputs_1",
            ] {
                let _ = create_output(output_path, format!("{}.txt", s));
            }
        }

        if self.housekeeping.recover {
            let Some((rec_states, _)) = self.recover(output_path) else {
                eprintln!("Could not recover input");
                return;
            };
            world_states = rec_states
        }
        //-------------------------------------------------------------------
        //                       Initialize time loop
        //-------------------------------------------------------------------
        let mut real_time = if !self.housekeeping.recover {
            t_form_min
        } else {
            0.0 // trecover would be here
        } - dtime;

        let mut isteps = 0;
        let n_steps = (self.grid.output_every / self.grid.time_step) as usize;

        for _ in 0..=n_time {
            real_time += dtime;
            let q_prim = self.primary_world.tidal_q.q_prim(&self.worlds, real_time);

            if self.primary_world.mass > 0.0 {
                let nmoons = self.n_moons();

                world_states
                    .iter_mut()
                    .for_each(|w| w.e_orb.max_assign(MIN_ECC));

                let mut resonance = vec![vec![0.0; nmoons]; nmoons];
                let mut p_capture = vec![vec![0.0; nmoons]; nmoons];
                let a_old: Vec<f64> = world_states.iter().map(|w| w.a_old).collect();
                let t_tide: Vec<f64> = self
                    .worlds
                    .iter()
                    .map(|w| w.t_reslock * MYR2SEC * real_time / (4568.2 * MYR2SEC))
                    .collect();

                for im in 0..nmoons {
                    let tzero_im = self.worlds[im].t_form * MYR2SEC;
                    if real_time >= tzero_im {
                        self.res_check(
                            im,
                            &world_states,
                            &mut resonance,
                            &mut p_capture,
                            real_time,
                            &a_old,
                            &t_tide,
                            q_prim,
                        );
                    }
                }

                for im in 0..nmoons {
                    let mut res_acct_for = world_states[im].res_acct_for.clone();
                    let res_acct_for_old = world_states[im].res_acct_for.clone();
                    self.resscreen(&resonance[im], &mut res_acct_for, &res_acct_for_old);
                    world_states[im].res_acct_for = res_acct_for;

                    for i in 0..nmoons {
                        if i != im && world_states[im].res_acct_for[i] == 0.0 {
                            world_states[i].res_acct_for[im] = 0.0;
                        }
                    }

                    let tzero_im = self.worlds[im].t_form * MYR2SEC;
                    if real_time >= tzero_im {
                        self.orbit(im, &mut world_states, dtime, real_time, q_prim);
                    }
                }

                world_states
                    .iter_mut()
                    .for_each(|w| w.e_orb.max_assign(MIN_ECC));
            }

            // Call Thermal logic
            self.thermal(&mut world_states, dtime, real_time, &GLOBAL_ARGS);

            isteps += 1;

            if self.primary_world.mass > 0.0 {
                if isteps % (n_steps / 10).max(1) == 0 {
                    for (im, world) in world_states.iter().enumerate() {
                        let orbit_output = [
                            real_time / GYR2SEC,
                            world.a_orb / KM2CM,
                            world.a_old / KM2CM,
                            world.e_orb,
                            world.i_orb * 180.0 / PI,
                            world.obl * 180.0 / PI,
                            world.spin,
                            world.h_old,
                            world.k_old,
                            if world.h_old != 0.0 {
                                let mut angle = (world.k_old / world.h_old).atan();
                                if world.h_old < 0.0 {
                                    angle += PI * world.k_old.signum();
                                }
                                angle * 180.0 / PI
                            } else {
                                0.0
                            },
                            (world.w_tide_tot + world.w_fluidtide_tot) / 1.0e7,
                            (world.w_tide_tot + world.w_fluidtide_tot)
                                / (11.5
                                    * self.worlds[im].planetary_rad.powi(5)
                                    * ((GCGS * self.primary_world.mass / world.a_orb.powi(3))
                                        .sqrt())
                                    .powi(5)
                                    * world.e_orb.powi(2)
                                    / GCGS),
                        ];
                        let _ = append_output(
                            output_path,
                            &format!("{}_{}_Orbit.csv", im, world.name),
                            &orbit_output,
                        );
                    }
                    let rebound_output = world_states
                        .iter()
                        .flat_map(|w| {
                            [
                                w.a_orb / KM2CM,
                                w.e_orb,
                                w.i_orb,
                                w.zones.last().unwrap().radius / KM2CM,
                                w.k2,
                                w.moi,
                                w.q_tide,
                            ]
                        })
                        .collect_vec();
                    let _ = append_output(output_path, "icydwarf_outputs_1.txt", &rebound_output);
                }

                if isteps == n_steps {
                    let mut mring = self.primary_world.ring.mass;
                    for w in &self.worlds {
                        if w.from_ring && real_time >= w.t_form * MYR2SEC {
                            mring -= w.mass();
                        }
                    }
                    let primary_output = [real_time / GYR2SEC, q_prim, mring * GRAM];
                    let _ = append_output(output_path, "Primary.txt", &primary_output);
                }
            }

            if isteps == n_steps {
                isteps = 0;

                for (im, world) in world_states.iter().enumerate() {
                    for zone in &world.zones {
                        let thermal_output = [
                            zone.radius / KM2CM,
                            zone.temp,
                            zone.mass_rock,
                            zone.mass_ice,
                            zone.mass_ammonia_solid,
                            zone.mass_water,
                            zone.mass_ammonia_liquid,
                            zone.nusselt,
                            0.0,
                            zone.kappa / 1.0e5,
                            zone.x_hydr,
                            zone.porosity,
                            zone.crack,
                            zone.tide_heat_rate,
                        ];
                        let _ = append_output(
                            output_path,
                            &format!("{}_{}_Thermal.csv", im, world.name),
                            &thermal_output,
                        );
                    }

                    let heat_output = [
                        real_time / GYR2SEC,
                        world.heat_radio * dtime,
                        world.heat_grav * dtime,
                        world.heat_serp * dtime,
                        world.heat_dehydr * dtime,
                        world.heat_tide * dtime,
                        world.heat_fluidtide * dtime,
                        world.cesq,
                        world.til_t,
                        world.w_fluidtide_tot,
                    ];
                    let _ = append_output(
                        output_path,
                        &format!("{}_{}_Heats.csv", im, world.name),
                        &heat_output,
                    );

                    let crack_depth = world
                        .zones
                        .iter()
                        .rev()
                        .find(|z| z.crack > 0.0)
                        .map_or(0.0, |z| z.radius / KM2CM);
                    let w_r = if world.m_cracked_rock < 0.000001 {
                        0.0
                    } else {
                        world.m_liq / world.m_cracked_rock
                    };
                    let crack_depth_wr_output = [real_time / GYR2SEC, crack_depth, w_r];
                    let _ = append_output(
                        output_path,
                        &format!("{}_{}_Crack_depth_WR.csv", im, world.name),
                        &crack_depth_wr_output,
                    );

                    for zone in &world.zones {
                        let mut stress_output = zone.stress;
                        stress_output[0] = zone.radius / KM2CM;
                        stress_output[10] = zone.frac_open;
                        stress_output[11] = zone.crack;
                        let _ = append_output(
                            output_path,
                            &format!("{}_{}_Crack_stresses.csv", im, world.name),
                            &stress_output,
                        );
                    }
                }
            }
        }
    }
}

#[allow(dead_code)]
impl IcyWorld {
    pub fn mass(&self) -> f64 {
        self.planetary_dens * 4.0 * FRAC_PI_3 * self.planetary_rad.powi(3)
    }

    pub fn rho_ice(&self) -> f64 {
        1.0 / (self.ammonia / XC) / (RHO_ADHS * GRAM)
    }

    pub fn frac_rock_mass(&self, world_spec: &WorldSpec) -> f64 {
        let rho_ice = self.rho_ice();
        (1. - rho_ice / self.planetary_dens)
            / (1. - rho_ice / self.hydr_init * world_spec.rho_rock_th()
                + (1. - self.hydr_init) * world_spec.rho_hydr_th())
    }

    pub fn frac_rock_volume(&self, world_spec: &WorldSpec) -> f64 {
        let frac_rock_mass = self.frac_rock_mass(world_spec);
        frac_rock_mass * self.planetary_dens
            / (self.hydr_init * world_spec.rho_hydr_th()
                + (1. - self.hydr_init * world_spec.rho_rock_th()))
    }
}

impl WorldSpec {
    pub fn rho_hydr_th(&self) -> f64 {
        self.rho_rock_hydr * GRAM
    }

    pub fn rho_rock_th(&self) -> f64 {
        self.rho_rock_dry * GRAM
    }
}

impl TidalQ {
    pub fn q_prim(&self, planets: &[IcyWorld], real_time: f64) -> f64 {
        let t_form_min = planets.iter().fold(f64::MAX, |acc, e| acc.min(e.t_form));
        let scaling = (real_time - t_form_min) / (TODAY - t_form_min);
        match self.mode {
            QMode::Lin => self.init + (self.today - self.init) * scaling,
            QMode::ExpDecay => {
                // self.init * ((self.today / self.init).ln() * scaling).exp()
                //
                // init * e^(ln(today / init) * scaling)
                // = init * (e^ln(today/init))^scaling
                // = init * (today/init) ^ scaling
                self.init * (self.today / self.init).powf(scaling)
            }
            QMode::ExpChange => {
                self.init + 1. - self.init * (self.today - self.init + 1.).powf(scaling)
            }
        }
    }
}

/// Calculate specific enthalpy of rock in erg per gram at a given temperature.
///
/// # Parameters
/// - `t`: The temperature in Kelvin.
///
/// # Returns
/// The specific enthalpy value of rock.
pub fn heat_rock(t: f64) -> f64 {
    if t > 1000.0 {
        EROCK_A * 275.0 * 275.0
            + (1000.0 - 275.0) * (EROCK_C + EROCK_D * 1000.0)
            + (EROCK_F) * (t - 1000.0)
    } else if t > 275.0 {
        EROCK_A * 275.0 * 275.0 + (t - 275.0) * (EROCK_C + EROCK_C * t)
    } else {
        EROCK_A * t * t
    }
}

/// Calculate specific enthalpy and phase mass fractions for ice mixtures.
///
/// # Parameters
/// - `t`: The temperature in Kelvin.
/// - `x`: The ammonia mass fraction.
/// - `x_salt`: The salinity mass fraction.
///
/// # Returns
/// An array containing specific enthalpy and phase mass fractions for water ice, ammonia ice, liquid water, and liquid ammonia.
pub fn heat_ice(t: f64, x: f64, x_salt: f64) -> [f64; 5] {
    let xb = XC * (2.0 / 95.0_f64).sqrt();

    // Determine transition temperatures based on concentration and salt content
    let (t_low, t_mid, t_high) = match (x <= xb, x_salt > 0.0) {
        (true, false) => (271.0, 273.0, 275.0),
        (true, true) => (248.0, 250.0, 252.0),
        (false, _) => {
            let tl = 273.0 - 95.0 * (x / XC).powi(2);
            (tl, 273.0, tl) // For x > xb, t_high = t_low = t_liq
        }
    };

    // Precalculate boundary energies to simplify the piecewise energy function
    let e_174 = 0.5 * QH2O * 174.0 * 174.0 + (x / XC) * 0.5 * (QADH - QH2O) * 174.0 * 174.0;
    let e_178 = {
        let t2 = 178.0;
        e_174
            + (1.0 - x / XC) * 0.5 * QH2O * (t2 * t2 - 174.0 * 174.0)
            + (x / XC) * (t2 - 174.0) / 4.0
                * (LADH
                    + (182.0 - t2) / 2.0 * QADH * 174.0
                    + XC * (t2 - 174.0) / 2.0 * CNH3L
                    + (1.0 - XC) * (t2 - 174.0) / 2.0 * CH2OL)
    };
    let e_low = {
        let r = ((t_mid - t_low) / 95.0).sqrt();
        e_178
            + 0.5 * QH2O * (t_low * t_low - 178.0 * 178.0)
            + x * (CNH3L - CH2OL) * (t_low - 178.0)
            + (x / XC)
                * (1.0 - r)
                * (LH2O / r + 2.0 * 95.0 * CH2OL - 2.0 * 95.0 * QH2O * t_mid
                    + 2.0 * QH2O * 95.0 * 95.0 / 3.0 * (1.0 + r + r * r))
    };
    let e_high = if t_low < t_high {
        e_low
            + (t_high - t_low) * x * (CNH3L + CH2OL * (1.0 / xb - 1.0))
            + (1.0 - x / xb) * (t_high - t_low) / 4.0 * (LH2O + 0.5 * CH2OL * (t_high - t_low))
    } else {
        e_low
    };

    // Helper function to calculate the specific internal energy e at a certain temperature t
    let e = |tp: f64| -> f64 {
        match tp {
            tp if tp <= 174.0 => 0.5 * QH2O * tp * tp + (x / XC) * 0.5 * (QADH - QH2O) * tp * tp,
            tp if tp <= 178.0 => {
                e_174
                    + (1.0 - x / XC) * 0.5 * QH2O * (tp * tp - 174.0 * 174.0)
                    + (x / XC) * (tp - 174.0) / 4.0
                        * (LADH
                            + (182.0 - tp) / 2.0 * QADH * 174.0
                            + XC * (tp - 174.0) / 2.0 * CNH3L
                            + (1.0 - XC) * (tp - 174.0) / 2.0 * CH2OL)
            }
            tp if tp <= t_low => {
                let r = ((t_mid - tp) / 95.0).sqrt();
                e_178
                    + 0.5 * QH2O * (tp * tp - 178.0 * 178.0)
                    + x * (CNH3L - CH2OL) * (tp - 178.0)
                    + (x / XC)
                        * (1.0 - r)
                        * (LH2O / r + 2.0 * 95.0 * CH2OL - 2.0 * 95.0 * QH2O * t_mid
                            + 2.0 * QH2O * 95.0 * 95.0 / 3.0 * (1.0 + r + r * r))
            }
            tp if tp <= t_high && t_low < t_high => {
                e_low
                    + (tp - t_low) * x * (CNH3L + CH2OL * (1.0 / xb - 1.0))
                    + (1.0 - x / xb) * (tp - t_low) / 4.0
                        * (LH2O + 0.5 * CH2OL * (tp - t_low) + 0.5 * QH2O * t_low * (t_high - tp))
            }
            _ => e_high + (x * CNH3L + (1.0 - x) * CH2OL) * (tp - t_high),
        }
    };

    let energy = e(t);

    // Calculate phase fractions based on current temperature using match
    let (gh2os, gadhs, gh2ol, gnh3l) = match t {
        t if t <= 174.0 => (1.0 - x / XC, x / XC, 0.0, 0.0),
        t if t <= 178.0 => {
            let f = (t - 174.0) / 4.0;
            (
                1.0 - x / XC,
                (x / XC) * (1.0 - f),
                (x / XC) * f * (1.0 - XC),
                (x / XC) * f * XC,
            )
        }
        t if t <= t_low => {
            let x_liq = XC * ((t_mid - t) / 95.0).sqrt();
            (1.0 - x / x_liq, 0.0, x / x_liq - x, x)
        }
        t if t <= t_high && t_low < t_high => {
            let f = (t - t_low) / 4.0;
            (
                (1.0 - x / xb) * (1.0 - f),
                0.0,
                (1.0 - x / xb) * f + (x / xb - x),
                x,
            )
        }
        _ => (0.0, 0.0, 1.0 - x, x),
    };

    [energy, gh2os, gadhs, gh2ol, gnh3l]
}
