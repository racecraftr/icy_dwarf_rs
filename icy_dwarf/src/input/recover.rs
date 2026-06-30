use crate::{
    consts::*,
    input::{Fracs, IcyDwarfInput},
    planet_system::{WorldState, ZoneState, heat_ice, heat_rock},
    thermal::ThermalOut,
    traits::float_traits::FloatExt,
};
use std::fs;
use std::path::Path;

#[derive(Default, Debug)]
pub struct OrbitOut {
    pub time: f64,
    pub a_orb: f64,
    pub a_old: f64,
    pub e_orb: f64,
    pub h_old: f64,
    pub k_old: f64,
    pub w_tide_tot: f64,
}

impl OrbitOut {
    pub fn from_line(line: &str) -> Option<Self> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 8 {
            return None;
        }
        Some(Self {
            time: parts[0].parse().ok()?,
            a_orb: parts[1].parse().ok()?,
            a_old: parts[2].parse().ok()?,
            e_orb: parts[3].parse().ok()?,
            h_old: parts[4].parse().ok()?,
            k_old: parts[5].parse().ok()?,
            w_tide_tot: parts[7].parse().ok()?,
        })
    }
}

#[derive(Default, Debug)]
pub struct CrackStressesOut {
    pub pressure: f64,
    pub p_pore: f64,
    pub p_hydr: f64,
    pub crack_size: f64,
}

impl CrackStressesOut {
    pub fn from_line(line: &str) -> Option<Self> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 10 {
            return None;
        }
        Some(Self {
            pressure: parts[1].parse().ok()?,
            p_pore: parts[5].parse().ok()?,
            p_hydr: parts[6].parse().ok()?,
            crack_size: parts[9].parse().ok()?,
        })
    }
}

impl IcyDwarfInput {
    pub fn recover(&self, output_folder: String) -> Option<(Vec<WorldState>, f64)> {
        let nr = self.grid.n_zones;
        let mut world_states = Vec::with_capacity(self.n_moons());
        let mut trecover = 0.0;
        let path = Path::new(&output_folder);

        for (idx, world) in self.worlds.iter().enumerate() {
            let name = &world.name;

            // Thermal
            let thermal_path = path.join(format!("{}_{}_Thermal.txt", idx, name));
            let thermal_str = fs::read_to_string(&thermal_path).ok()?;
            let thermal_lines: Vec<_> = thermal_str.lines().filter(|l| !l.is_empty()).collect();
            let mut thermal_out = Vec::with_capacity(nr);
            for line in thermal_lines.iter().rev().take(nr).rev() {
                thermal_out.push(ThermalOut::from_line(line)?);
            }

            // Orbit
            let mut orbit_out = OrbitOut::default();
            if self.primary_world.mass > 0.0 {
                let orbit_path = path.join(format!("{}_{}_Orbit.txt", idx, name));
                if let Ok(orbit_str) = fs::read_to_string(&orbit_path) {
                    let orbit_lines: Vec<_> = orbit_str.lines().filter(|l| !l.is_empty()).collect();
                    if let Some(line) = orbit_lines.last() {
                        orbit_out = OrbitOut::from_line(line)?;
                    }
                }
            }
            if idx == 0 {
                trecover = orbit_out.time * GYR2SEC;
            }

            // Crack_stresses
            let crack_path = path.join(format!("{}_{}_Crack_stresses.txt", idx, name));
            let mut crack_out = Vec::with_capacity(nr);
            if let Ok(crack_str) = fs::read_to_string(&crack_path) {
                let crack_lines: Vec<_> = crack_str.lines().filter(|l| !l.is_empty()).collect();
                for line in crack_lines.iter().rev().take(nr).rev() {
                    crack_out.push(CrackStressesOut::from_line(line)?);
                }
            } else {
                for _ in 0..nr {
                    crack_out.push(CrackStressesOut::default());
                }
            }

            let mut zones = Vec::with_capacity(nr);
            let mut current_r = 0.0;

            for ir in 0..nr {
                let t_out = &thermal_out[ir];
                let c_out = &crack_out[ir];

                let dr = t_out.radius_km - current_r;
                let mass_total = t_out.mass_total();

                let mut x = 0.0;
                let x_salt = if world.briny { 0.01 } else { 0.0 };

                if mass_total > 0.0 {
                    let Fracs(frock, _, fadhs, _, fnh3l) = t_out.fracs();

                    if frock < 1.0 {
                        let gadhs = fadhs / (1.0 - frock);
                        let gnh3l = fnh3l / (1.0 - frock);
                        x = gnh3l + XC * gadhs;
                    }
                }

                let e_rock = t_out.mass_rock * heat_rock(t_out.temp_kelvin);
                let [e_ice, ..] = heat_ice(t_out.temp_kelvin, x, x_salt);
                let e_h2os = t_out.mass_ice * e_ice;
                let e_slush =
                    (t_out.mass_ammonia_solid + t_out.mass_water + t_out.mass_ammonia_liquid)
                        * e_ice;

                zones.push(ZoneState {
                    radius: t_out.radius_km,
                    dr,
                    temp: t_out.temp_kelvin,
                    temp_old: t_out.temp_kelvin,
                    mass_total,
                    mass_rock: t_out.mass_rock,
                    mass_rock_init: t_out.mass_rock,
                    mass_ice: t_out.mass_ice,
                    mass_ammonia_solid: t_out.mass_ammonia_solid,
                    mass_water: t_out.mass_water,
                    mass_ammonia_liquid: t_out.mass_ammonia_liquid,
                    energy_total: e_rock + e_h2os + e_slush,
                    porosity: t_out.porosity,
                    pressure: c_out.pressure * MPA,
                    x_hydr: t_out.deg_of_hydr,
                    x_hydr_old: t_out.deg_of_hydr,
                    kappa: t_out.thermal_cond * 1.0e5,
                    nusselt: t_out.nusselt_num,
                    crack: if t_out.crack { 1.0 } else { 0.0 },
                    crack_size: c_out.crack_size,
                    p_pore: c_out.p_pore * MPA,
                    p_hydr: c_out.p_hydr * MPA,
                    act: [0.0; 3],
                });

                current_r = t_out.radius_km;
            }

            let aorb = orbit_out.a_orb * KM2CM;
            let norb = (GCGS * self.primary_world.mass / aorb.powi(3))
                .sqrt()
                .unwrap_or_nan(0.);

            world_states.push(WorldState {
                name: name.clone(),
                a_orb: aorb,
                e_orb: orbit_out.e_orb,
                i_orb: world.orb_i_init,
                obl: world.orb_o_init,
                n_orb: norb,
                lambda: 0.32 + (idx as f64) * 0.77,
                omega: 0.58 + (idx as f64) * 0.27,
                spin: norb,
                k2: 0.001,
                q_tide: 1.0e10,
                moi: 0.4,
                zones,
                res_acct_for: vec![0.0; self.n_moons()],
                h_old: orbit_out.h_old,
                k_old: orbit_out.k_old,
                a_old: orbit_out.a_old * KM2CM,
                cs_ee_old: 0.0,
                cs_eep_old: 0.0,
                cr_e_old: 0.0,
                cr_ep_old: 0.0,
                cr_ee_old: 0.0,
                cr_eep_old: 0.0,
                cr_epep_old: 0.0,
                w_tide_tot: orbit_out.w_tide_tot * 1.0e7,
                w_fluidtide_tot: 0.0,
            });
        }

        Some((world_states, trecover))
    }
}
