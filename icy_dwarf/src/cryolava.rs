use crate::{consts::GRAM, input::IcyDwarfInput, thermal::ThermalOut};
use extendr_api::prelude::*;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

pub const N_SPECIES: i32 = 10;
pub const N_ITER_MAX: i32 = 1000;
pub const NEWT_RAPH_THRESH: f64 = 1.5e10;
pub const X_SUP_BOUND: f64 = 1.0e3;
pub const K_IC_ICE: f64 = 0.15e6;
pub const K_IC_CRUST: f64 = 0.5e6;

impl IcyDwarfInput {
    pub fn calculate_mass_liquid(&self, thermal_out: &[Vec<ThermalOut>], t: usize) -> f64 {
        thermal_out.iter().map(|v| v[t].mass_water * GRAM).sum()
    }

    pub fn calculate_seafloor(&self, thermal_out: &[Vec<ThermalOut>], t: usize) -> i32 {
        if t >= self.grid.output_time_step() {
            println!("[[CALCULATE SEAFLOOR]] :: t >= NT");
        }
        let mut r = 0;
        while r < self.grid.n_zones {
            if thermal_out[r][t].mass_rock <= 0.0 {
                return r as i32 - 1;
            }
            r += 1;
        }
        self.grid.n_zones as i32 - 1
    }

    pub fn calculate_pressure_cryolava(
        &self,
        thermal_out: &[Vec<ThermalOut>],
        t: usize,
    ) -> Vec<f64> {
        let nr = self.grid.n_zones;
        let mut pressure = vec![0.0; nr];
        let mut m = vec![0.0; nr];
        let mut frock = vec![0.0; nr];
        let mut fh2os = vec![0.0; nr];
        let mut fh2ol = vec![0.0; nr];
        let mut fadhs = vec![0.0; nr];
        let mut fnh3l = vec![0.0; nr];
        let mut g = vec![0.0; nr];

        let mut r_bound = vec![0.0; nr + 1];
        r_bound[0] = 0.0;
        for ir in 0..nr {
            r_bound[ir + 1] = thermal_out[ir][t].radius_km;
        }

        for ir in 0..nr {
            let zone = &thermal_out[ir][t];
            let d_m = zone.mass_total();
            frock[ir] = zone.mass_rock / d_m;
            fh2os[ir] = zone.mass_ice / d_m;
            fh2ol[ir] = zone.mass_water / d_m;
            fadhs[ir] = zone.mass_ammonia_solid / d_m;
            fnh3l[ir] = zone.mass_ammonia_liquid / d_m;

            m[ir] = if ir > 0 { m[ir - 1] + d_m } else { d_m };
        }

        for ir in 0..nr {
            g[ir] = crate::consts::G * m[ir] * crate::consts::GRAM
                / (r_bound[ir + 1] * r_bound[ir + 1])
                * crate::consts::KM2CM
                * crate::consts::KM2CM
                / crate::consts::KM
                / crate::consts::KM;
        }

        pressure[nr - 1] = 0.0;
        let rho_hydr = self.world_spec.rho_rock_hydr;
        let rho_dry = self.world_spec.rho_rock_dry;

        for ir in (0..nr - 1).rev() {
            let zone_curr = &thermal_out[ir][t];
            pressure[ir] = pressure[ir + 1]
                + 0.5 * (g[ir + 1] + g[ir]) * (r_bound[ir + 2] - r_bound[ir + 1])
                    / crate::consts::KM2CM
                    * crate::consts::KM
                    * (frock[ir + 1]
                        * (zone_curr.deg_of_hydr * rho_hydr
                            + (1.0 - zone_curr.deg_of_hydr) * rho_dry)
                        + fh2os[ir + 1] * crate::consts::RHO_H2OS
                        + fh2ol[ir + 1] * crate::consts::RHO_H2OL
                        + fadhs[ir + 1] * crate::consts::RHO_ADHS
                        + fnh3l[ir + 1] * crate::consts::RHO_NH3L);
        }

        pressure
    }

    pub fn cryolava(&self, thermal_out: &[Vec<ThermalOut>], path: &str) -> Result<(), String> {
        let t = self.t_cryo() as usize;

        let m_liq = self.calculate_mass_liquid(thermal_out, t);
        if m_liq <= 0.0 {
            println!("Cryolava: No liquid at t_cryolava={}", t);
            return Err(format!("Cryolava: No liquid at t_cryolava={}", t));
        }

        let nr = self.grid.n_zones;
        let pressure = self.calculate_pressure_cryolava(thermal_out, t);

        let r_seafloor = self.calculate_seafloor(thermal_out, t);
        let r_seafloor_idx = r_seafloor.max(0) as usize;

        let temp = thermal_out[r_seafloor_idx][t].temp_kelvin;

        // Find the base of the crust r_diff
        let mut r_diff = nr - 2;
        let mut r = nr - 2;
        while r > r_seafloor_idx {
            if thermal_out[r][t].mass_rock <= 0.0 {
                r_diff = r;
                break;
            }
            r -= 1;
        }

        let r_hydrostatic = (r_seafloor_idx as f64
            + crate::consts::RHO_H2OS / crate::consts::RHO_H2OL
                * (r_diff as f64 - r_seafloor_idx as f64))
            .floor() as i32;

        let species = [
            "H2", "CH4", "CH3OH", "CO", "CO2", "NH3", "N2", "H2S", "SO2", "Ar",
        ];
        let wrt_h2o = [
            1.0e-5, 0.01, 0.03, 0.2, 0.1, 0.01, 0.01, 0.005, 2.0e-5, 0.001,
        ];

        let mut abundances = vec![0.0; N_SPECIES as usize];
        for i in 0..N_SPECIES as usize {
            abundances[i] = wrt_h2o[i] * m_liq / 0.018;
            println!("{} = {:e} mol/kg", species[i], abundances[i] / m_liq);
        }

        extendr_engine::start_r();
        let _ = R!(r#"
            if (!requireNamespace("CHNOSZ", quietly = TRUE)) {
                stop("CHNOSZ package not installed in R")
            }
            library(CHNOSZ, quietly = TRUE)
            data(thermo)
            add.OBIGT("SUPCRT92")
            get_logK_safely <- function(species, state, T_c, P) {
                val <- tryCatch({
                    res <- subcrt(species, state, T = T_c, P = P)
                    res$out[[1]]$logK[1]
                }, error = function(e) {
                    NA
                })
                if (is.na(val)) 0.0 else val
            }
        "#)
        .map_err(|e| format!("Failed to initialize CHNOSZ in R: {:?}", e))?;

        let nr_active = nr - r_seafloor_idx;
        let mut partial_p = vec![vec![0.0; N_SPECIES as usize]; nr_active];
        let mut molalities = vec![vec![0.0; N_SPECIES as usize]; nr_active];
        let mut x_vap_table = vec![vec![0.0; 6]; nr_active];
        let mut k_rxn = vec![0.0; N_SPECIES as usize];

        println!("Cryolava: Calculating species molalities...");
        let r_p = self.worlds[0].planetary_rad;
        let min_temp_chnosz = self.subroutines.cryo.min_temp_chnosz;
        let warnings = self.housekeeping.warnings;

        for r in 0..nr_active {
            let abs_r = r + r_seafloor_idx;

            let mut d_int_prec = 0.0;
            let mut p_integral = 0.0;
            for i in r_seafloor_idx..abs_r {
                let mut m_inf = 0.0;
                for u in 0..i {
                    m_inf += thermal_out[u][t].mass_total();
                }
                let r_m = thermal_out[i][t].radius_km / 100.0;
                let d_int = crate::consts::RHO_H2OL * crate::consts::G / (r_m * r_m);

                let dr_m = if i > 0 {
                    (thermal_out[i][t].radius_km - thermal_out[i - 1][t].radius_km) / 100.0
                } else {
                    thermal_out[i][t].radius_km / 100.0
                };

                p_integral += (d_int + d_int_prec) / 2.0 * m_inf * crate::consts::GRAM * dr_m;
                d_int_prec = d_int;
            }

            let p_gas = pressure[r_seafloor_idx] - p_integral;

            // Use CHNOSZ to get reaction constants at given T and P
            let p_local = pressure[abs_r];
            let p_bar = p_local / crate::consts::BAR;

            let mut chnosz_t = temp;
            if chnosz_t < min_temp_chnosz {
                if warnings && r == 0 {
                    println!(
                        "Cryolava: T={} K below minimum temp for CHNOSZ. Using T={} K instead",
                        chnosz_t, min_temp_chnosz
                    );
                }
                chnosz_t = min_temp_chnosz;
            }
            let temp_c = chnosz_t - crate::consts::KELVIN;

            for i in 0..N_SPECIES as usize {
                let species_name = species[i];
                let log_reactant =
                    R!(r#"get_logK_safely({{species_name}}, "g", {{temp_c}}, {{p_bar}})"#)
                        .map_err(|e| e.to_string())?
                        .as_real()
                        .ok_or_else(|| {
                            format!(
                                "Failed to get logK for {} g at T_c={}, P={}",
                                species_name, temp_c, p_bar
                            )
                        })?;
                let log_product =
                    R!(r#"get_logK_safely({{species_name}}, "aq", {{temp_c}}, {{p_bar}})"#)
                        .map_err(|e| e.to_string())?
                        .as_real()
                        .ok_or_else(|| {
                            format!(
                                "Failed to get logK for {} aq at T_c={}, P={}",
                                species_name, temp_c, p_bar
                            )
                        })?;

                k_rxn[i] = 10f64.powf(log_product - log_reactant);
                if !k_rxn[i].is_finite() || k_rxn[i] < 0.0 {
                    println!(
                        "Cryolava: Error calculating K_rxn[{}]={} at t={}, r={}",
                        i, k_rxn[i], t, r
                    );
                }
            }

            let mut x_inf = 0.;
            let mut x_sup = X_SUP_BOUND;

            let f_inf = f(
                p_gas / crate::consts::BAR,
                m_liq,
                &abundances,
                &k_rxn,
                x_inf,
            );
            let f_sup = f(
                p_gas / crate::consts::BAR,
                m_liq,
                &abundances,
                &k_rxn,
                x_sup,
            );

            let mut x_vap_val = 0.0;
            if f_inf * f_sup > 0.0 {
                println!(
                    "Cryolava: No physical solution at depth {} km: P_gas={} bar either negative or too high",
                    (nr - abs_r) as f64 * r_p / nr as f64,
                    p_gas / crate::consts::BAR
                );
                x_vap_table[r][0] = (nr - abs_r) as f64 * r_p / nr as f64;
            } else {
                if f_inf > 0.0 {
                    std::mem::swap(&mut x_inf, &mut x_sup);
                }
                x_vap_val = 0.5 * (x_inf + x_sup);
                let mut dx_old = (x_inf - x_sup).abs();
                let mut dx = dx_old;

                let mut f_x = f(
                    p_gas / crate::consts::BAR,
                    m_liq,
                    &abundances,
                    &k_rxn,
                    x_vap_val,
                );
                let mut f_prime_x = f_prime(
                    p_gas / crate::consts::BAR,
                    m_liq,
                    &abundances,
                    &k_rxn,
                    x_vap_val,
                )
                .unwrap_or(0.0);

                let mut n_iter = 0;
                while f_x.abs() > NEWT_RAPH_THRESH {
                    if (((x_vap_val - x_sup) * f_prime_x - f_x)
                        * ((x_vap_val - x_inf) * f_prime_x - f_x)
                        > 0.0)
                        || (2.0 * f_x).abs() > (dx_old * f_prime_x).abs()
                    {
                        dx_old = dx;
                        dx = 0.5 * (x_sup - x_inf);
                        x_vap_val = x_inf + dx;
                    } else {
                        dx_old = dx;
                        dx = f_x / f_prime_x;
                        x_vap_val -= dx;
                    }

                    f_x = f(
                        p_gas / crate::consts::BAR,
                        m_liq,
                        &abundances,
                        &k_rxn,
                        x_vap_val,
                    );
                    f_prime_x = f_prime(
                        p_gas / crate::consts::BAR,
                        m_liq,
                        &abundances,
                        &k_rxn,
                        x_vap_val,
                    )
                    .unwrap_or(0.0);

                    if f_x < 0.0 {
                        x_inf = x_vap_val;
                    } else {
                        x_sup = x_vap_val;
                    }

                    n_iter += 1;
                    if n_iter >= N_ITER_MAX {
                        if warnings {
                            println!(
                                "Cryolava: could not converge towards a solution of chemical abundances after {} iterations",
                                N_ITER_MAX
                            );
                        }
                        break;
                    }
                }

                x_vap_table[r][0] = (nr - abs_r) as f64 * r_p / nr as f64;
                x_vap_table[r][1] = p_gas / crate::consts::BAR;
                x_vap_table[r][2] = x_vap_val * crate::consts::RHO_H2OL * crate::consts::R_G * temp
                    / crate::consts::BAR;
                x_vap_table[r][3] = crate::consts::RHO_H2OL / (1.0 + x_vap_table[r][2]);

                let m_inf_stress = if r > 0 {
                    let limit = r_seafloor_idx + r - 1;
                    thermal_out
                        .iter()
                        .take(limit)
                        .fold(0., |n, v| n + v[t].mass_total())
                } else {
                    0.
                };

                let r_m_stress = thermal_out[abs_r][t].radius_km / 100.0;
                let depth_m = r as f64 * r_p / nr as f64 * crate::consts::KM;

                x_vap_table[r][4] = -(x_vap_table[r][3] - crate::consts::RHO_H2OS)
                    * 2.0
                    * crate::consts::G
                    * m_inf_stress
                    * crate::consts::GRAM
                    / (r_m_stress * r_m_stress)
                    * depth_m.powf(1.5)
                    / crate::consts::PI_GREEK.sqrt();

                let is_in_ice = abs_r <= r_diff;
                let threshold = if is_in_ice { K_IC_ICE } else { K_IC_CRUST };

                x_vap_table[r][5] = if x_vap_table[r][4] > threshold {
                    1.0
                } else {
                    0.0
                };

                println!(
                    "X_VAP = {} and x_vap = V_gas/V_liq = {} found after {} iterations",
                    x_vap_val, x_vap_table[r][2], n_iter
                );
            }

            for i in 0..N_SPECIES as usize {
                molalities[r][i] = abundances[i] / (m_liq * (1.0 + x_vap_val / k_rxn[i]));
                partial_p[r][i] = molalities[r][i] / k_rxn[i];
            }
        }

        write_output(&molalities, path, "Outputs/Cryolava_molalities.txt")?;
        write_output(&partial_p, path, "Outputs/Cryolava_partialP.txt")?;
        write_output(&x_vap_table, path, "Outputs/Cryolava_xvap.txt")?;

        let nr_f = nr as f64;
        println!(
            "\n Seafloor @ radius {} km",
            r_seafloor_idx as f64 * r_p / nr_f
        );
        println!(
            " Hydrostatic level in ice @ radius {} km",
            r_hydrostatic as f64 * r_p / nr_f
        );
        if r_diff < nr - 2 {
            println!(
                " Crust starts at R_diff = {} km",
                r_diff as f64 * r_p / nr_f
            );
        } else {
            println!(" No crust");
        }

        let mut out_path = PathBuf::from(path);
        out_path.push("Outputs");
        println!(
            "\nOutputs successfully generated in {:?} directory:",
            out_path
        );
        println!(
            "1. Molalities vs. depth at t={} in mol kg-1: Cryolava_molalities.txt",
            t
        );
        println!(
            "2. Partial pressures vs. depth at t={} in bar: Cryolava_partialP.txt",
            t
        );
        println!(
            "3. Volumic vapor fraction x_vap vs. P_gas at t={}: Cryolava_xvap.txt",
            t
        );

        Ok(())
    }
}

fn write_output(data: &[Vec<f64>], base_path: &str, relative_file: &str) -> Result<(), String> {
    let mut path = PathBuf::from(base_path);
    path.push(relative_file);

    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Failed to create directories for {:?}: {:?}", path, e))?;
    }

    let mut file =
        File::create(&path).map_err(|e| format!("Failed to create file {:?}: {:?}", path, e))?;

    for row in data {
        let mut line = String::new();
        for val in row {
            line.push_str(&format!("{} \t", val));
        }
        writeln!(file, "{}", line)
            .map_err(|e| format!("Failed to write to file {:?}: {:?}", path, e))?;
    }

    Ok(())
}

fn f(p: f64, mliq: f64, abundances: &[f64], k_rxn: &[f64], x: f64) -> f64 {
    let mut f_x = p * mliq * k_rxn.iter().map(|&n| n + x).product::<f64>();
    for (i, &abundance) in abundances.iter().enumerate() {
        let lhs = k_rxn
            .iter()
            .enumerate()
            .filter_map(|(j, &n)| if j != i { Some(n + x) } else { None })
            .product::<f64>();
        f_x -= abundance * lhs;
    }
    f_x
}

fn f_prime(p: f64, mliq: f64, abundances: &[f64], k_rxn: &[f64], x: f64) -> Option<f64> {
    if abundances.len() != k_rxn.len() {
        return None;
    }
    let mut f_prime_x = p * mliq;
    let mut rhs = 0.0;
    for i in 0..k_rxn.len() {
        rhs += k_rxn
            .iter()
            .enumerate()
            .filter_map(|(j, &n)| if j != i { Some(n + x) } else { None })
            .product::<f64>();
    }

    f_prime_x *= rhs;
    let mut lhs = 0.0;
    for (i, abundance) in abundances.iter().enumerate() {
        let mut lhs2 = 0.0;
        for j in (0..k_rxn.len()).filter(|&j| j != i) {
            lhs2 += k_rxn
                .iter()
                .enumerate()
                .filter_map(|(u, &n)| if u != j && u != i { Some(n + x) } else { None })
                .product::<f64>();
        }
        lhs += lhs2 * abundance;
    }
    Some(f_prime_x - lhs)
}
