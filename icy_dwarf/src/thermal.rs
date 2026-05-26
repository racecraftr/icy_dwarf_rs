use std::{cell::OnceCell, fs, sync::Once};

use crate::{
    consts::{GRAM, KM, KM2CM, RHO_ADHS, RHO_H2OL, RHO_NH3L},
    input::IcyDwarfInput,
    planet_system::{RHO_ADHS_TH, RHO_H2OL_TH, RHO_NH3L_TH, WorldState},
};
use faer::{
    Mat,
    linalg::solvers::{DenseSolveCore, Svd},
    prelude::Solve,
};
use num::complex::Complex64;

pub const N_THERM: usize = 14;

impl IcyDwarfInput {
    pub fn thermal(&self, world_state: &mut WorldState) {
        let tilt_scale: f64 = 1.0;
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

fn svd(mat: &Vec<Vec<f64>>) -> Option<Svd<f64>> {
    let rows = mat.len();
    let cols = mat[0].len();
    Mat::from_fn(rows, cols, |i, j| mat[i][j]).svd().ok()
}

/// Gauss jordan elimination on a complex matrix.
/// Given a matrix and a vector,
/// Uses full pivoting to solve Ax = b.
/// Returns both x and the inverse of A.
fn gauss_jordan(
    mat: &Vec<Vec<Complex64>>,
    b: &Vec<Complex64>,
) -> Option<(Mat<Complex64>, Mat<Complex64>)> {
    let rows = mat.len();
    let cols = mat[0].len();
    if rows != cols || b.len() != cols {
        return None; // avoids unneccessary computation 
    }
    let A = Mat::from_fn(rows, cols, |i, j| mat[i][j]);
    if A.determinant() == Complex64::from(0.0) {
        // matrix is not invertible
        return None;
    }
    let b = Mat::from_fn(b.len(), 1, |i, _| b[i]);
    let lu = A.full_piv_lu();
    let x = lu.solve(b);
    let a_inv = lu.inverse();
    Some((x, a_inv))
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_svd() {
        let mat: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0, 0.0, 2.0],
            vec![0.0, 0.0, 3.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0, 0.0, 0.0],
        ];
        let Some(svd) = svd(&mat) else {
            assert!(false);
            return;
        };
        println!("S = {:?}", svd.S()); // S represents Sigma (the diagonal entries)
        println!("U = {:?}", svd.U()); // U represents U
        println!("V = {:?}", svd.V().conjugate()); // V represents (V^*)^*
    }
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

    pub fn vol(&self, input: &IcyDwarfInput) -> ThermVol {
        let vol_rock = self.mass_rock
            / (self.deg_of_hydr * input.world_spec.rho_hydr_th()
                + (1.0 - self.deg_of_hydr) * input.world_spec.rho_rock_th());
        let vol_ice = self.mass_ice / RHO_H2OL_TH;
        let vol_adhs = self.mass_ammonia_solid / RHO_ADHS_TH;
        let vol_water = self.mass_water / RHO_H2OL_TH;
        let vol_nh3l = self.mass_ammonia_liquid / RHO_NH3L_TH;
        (
            vol_rock + vol_adhs + vol_water + vol_nh3l,
            (vol_rock, vol_ice, vol_adhs, vol_water, vol_nh3l),
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
