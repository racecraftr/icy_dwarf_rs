use std::fs;

use crate::{
    FloatMat,
    consts::KM2CM,
    input::{Fracs, IcyDwarfInput},
    planet_system::{RHO_ADHS_TH, RHO_H2OL_TH, RHO_H2OS_TH, RHO_NH3L_TH, WorldState},
    to_faer_mat,
};
use faer::{
    Mat,
    linalg::solvers::{DenseSolveCore, Svd},
    prelude::Solve,
};
use num::complex::Complex64;

pub const N_THERM: usize = 14;

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

    fn calculate_pressure(&self, world: &mut WorldState) {
        let mut cumulative_mass = 0.0;
        let mut gravity = vec![0.0; self.grid.n_zones];

        for (i, zone) in world.zones.iter().enumerate() {
            cumulative_mass += zone.mass_total;
            gravity[i] = crate::consts::GCGS * cumulative_mass * crate::consts::GRAM
                / zone.radius.powi(2)
                * crate::consts::KM2CM.powi(2)
                / crate::consts::KM.powi(2);
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
            let (frock, fh2os, fadhs, fh2ol, fnh3l) = next_zone.fracs();

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

        for (ir, zone) in world.zones.iter_mut().enumerate() {
            let v_rock = zone.mass_rock
                / (zone.x_hydr * self.world_spec.rho_hydr_th()
                    + (1.0 - zone.x_hydr) * self.world_spec.rho_rock_th());
            let d_vol =
                4.0 / 3.0 * std::f64::consts::PI * (r_old[ir + 1].powi(3) - r_old[ir].powi(3));

            let x_ice = 1.0 - v_rock / d_vol;
            let c_rate = creep(zone.temp, zone.pressure, x_ice, zone.porosity, zone.x_hydr);
            zone.porosity -= dtime * (1.0 - zone.porosity) * c_rate;
            if zone.porosity < 0.0 {
                zone.porosity = 0.0;
            }
            if zone.mass_rock < 0.01 && zone.mass_water > 0.01 {
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

    fn tide(&self) {
        let nr = self.grid.n_zones;
        let mut mean_layer_densities = vec![0.0; nr];
        let mut layer_masses = vec![0.0; nr];
        let mut y_tide = vec![vec![Complex64::from(0.0); 6]; nr];
    }

    fn hydrate(&self) {}

    fn prop_mtx(&self) {}
}

pub fn creep(t: f64, p: f64, x_ice: f64, porosity: f64, x_hydr: f64) -> f64 {
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

    if x_ice > 0.3 {
        creep_rate_ice
    } else {
        let t_eff = if t > 140.0 { t } else { 140.0 };
        let creep_rate_hydr = 10.0_f64.powf(5.62)
            * eff_p.powi(1)
            * D_FLOW_LAW.powi(-3)
            * (-240.0e3 / (R_G * t_eff)).exp();
        let creep_rate_dry = 10.0_f64.powf(5.25)
            * eff_p.powi(1)
            * D_FLOW_LAW.powf(-2.98)
            * (-261.0e3 / (R_G * t_eff)).exp();

        (((0.3 - x_ice) * (x_hydr * creep_rate_hydr + (1.0 - x_hydr) * creep_rate_dry).ln()
            + x_ice * creep_rate_ice.ln())
            / 0.3)
            .exp()
    }
}

/// Calculates the single value decomposition of a mxn matrix.
/// Retunrs, in order: U, Sigma, V^*.
fn svd(mat: &Vec<Vec<f64>>) -> Option<(FloatMat, FloatMat, FloatMat)> {
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

#[cfg(test)]
mod gauss_jordan_tests {
    use super::*;

    #[test]
    fn test_1() {}
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
    fn from_line(ln: &str) -> Option<Self> {
        let parts = ln.trim().split_whitespace().collect::<Vec<_>>();
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
        self.mass_rock + self.mass_ice + self.mass_ammonia_solid + self.mass_ammonia_liquid
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
