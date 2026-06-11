use crate::consts::*;
use extendr_api::prelude::*;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

/// Helper function to write output matrix to a text file in a format readable by the main code.
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

/// Calculate the integral part of equation (4) of Vance et al. (2007)
/// as a function of flaw size `a_var`, and determine the optimal `a_var`
/// that maximizes stress intensity `K_I` for various `deltaT` and `P`.
/// Outputs `Crack_integral.txt` and `Crack_aTP.txt`.
pub fn a_tp(path: &str, warnings: bool) -> Result<(), String> {
    let mut integral = vec![vec![0.0; 2]; INT_SIZE as usize];

    for j in 0..INT_SIZE as usize {
        let a_var = A_VAR_MAX * (j + 1) as f64 / INT_SIZE as f64;
        integral[j][0] = a_var;
        integral[j][1] = 0.0;

        let mut d_int_prec = 0.0;
        for i in 0..(INT_STEPS - 1) as usize {
            let x = a_var / INT_STEPS as f64 * (i + 1) as f64;

            // Calculate normal stress on grain boundary, eq (3) of Vance et al. (2007)
            let term1 =
                4.0 * L_SIZE * L_SIZE / (4.0 * L_SIZE * L_SIZE + (2.0 * L_SIZE - x).powi(2));
            let term2 = 4.0 * L_SIZE * L_SIZE / (4.0 * L_SIZE * L_SIZE + x.powi(2));
            let term3 = ((2.0 * L_SIZE - x) / x).ln();
            let term4 = 0.5
                * ((4.0 * L_SIZE * L_SIZE + (2.0 * L_SIZE - x).powi(2))
                    / (4.0 * L_SIZE * L_SIZE + x.powi(2)))
                .ln();

            let sigma_yy = term1 - term2 + term3 - term4;
            let d_int = sigma_yy * x.sqrt() / (a_var - x).sqrt();

            integral[j][1] += (d_int + d_int_prec) / 2.0 * 1.0 / INT_STEPS as f64 * a_var;
            d_int_prec = d_int;
        }
    }

    let mut delta_t = 0.0;
    let mut p_pa = 0.0;
    let mut a_tp_data = vec![vec![0.0; SIZEA_TP]; SIZEA_TP];
    let mut k_i = vec![vec![0.0; 2]; INT_SIZE as usize];

    for t in 0..SIZEA_TP {
        for p in 0..SIZEA_TP {
            for j in 0..(INT_SIZE - 1) as usize {
                let a_var = A_VAR_MAX / INT_SIZE as f64 * (j + 1) as f64;
                k_i[j][0] = a_var;

                let term_intensity = (2.0 / (PI_GREEK * a_var)).sqrt();
                let e_young_avg = 0.5 * (E_YOUNG_OLIV + E_YOUNG_SERP);
                let nu_poisson_avg = 0.5 * (NU_POISSON_OLIV + NU_POISSON_SERP);
                let den = 2.0 * PI_GREEK * (1.0 - nu_poisson_avg * nu_poisson_avg);

                k_i[j][1] = term_intensity * integral[j][1] * e_young_avg * DELTA_ALPHA / den
                    * delta_t
                    - p_pa * (PI_GREEK * a_var).sqrt();
            }

            let mut k_i_max = k_i[0][1];
            let mut k_i_max_a = k_i[0][0];
            for j in 0..(INT_SIZE - 1) as usize {
                if k_i[j][1] > k_i_max {
                    k_i_max = k_i[j][1];
                    k_i_max_a = k_i[j][0];
                }
            }

            if k_i_max_a < A_MIN {
                k_i_max_a = A_MIN;
            }
            a_tp_data[t][p] = k_i_max_a;
            p_pa += P_STEP;
        }
        p_pa = 0.0;
        delta_t += DELTA_T_STEP;
    }

    write_output(&integral, path, "Data/Crack_integral.txt")?;
    write_output(&a_tp_data, path, "Data/Crack_aTP.txt")?;

    if warnings {
        println!(
            "\n Outputs successfully generated in {}/Data/ directory:",
            path
        );
        println!("1. Crack_integral.txt");
        println!("2. Crack_aTP.txt");
    }

    Ok(())
}

/// Calculate the thermal expansivity (`alpha`) and compressibility (`beta`) of water
/// over a range of T and P using CHNOSZ.
/// Outputs `Crack_alpha.txt` and `Crack_beta.txt`.
pub fn crack_water_chnosz(path: &str, warnings: bool) -> Result<(), String> {
    let mut tempk = TEMPK_MIN;
    let mut p_bar = P_BAR_MIN;

    let mut alpha = vec![vec![0.0; SIZEA_TP]; SIZEA_TP];
    let mut beta = vec![vec![0.0; SIZEA_TP]; SIZEA_TP];

    extendr_engine::start_r();
    let _ = R!(r#"
        if (!requireNamespace("CHNOSZ", quietly = TRUE)) {
            stop("CHNOSZ package not installed in R")
        }
        library(CHNOSZ, quietly = TRUE)
        data(thermo)
        get_alpha <- function(T, P) {
            val <- as.numeric(water.SUPCRT92("alpha", T = T, P = P))
            if (is.na(val)) 0.0 else val
        }
        get_beta <- function(T, P) {
            val <- as.numeric(water.SUPCRT92("beta", T = T, P = P))
            if (is.na(val)) 0.0 else val
        }
    "#)
    .map_err(|e| format!("Failed to initialize CHNOSZ in R: {:?}", e))?;

    for t in 0..SIZEA_TP {
        for p in 0..SIZEA_TP {
            let alpha_val = R!(r#"get_alpha({{tempk}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| format!("Failed to get alpha at T={}, P={}", tempk, p_bar))?;

            let beta_val = R!(r#"get_beta({{tempk}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| format!("Failed to get beta at T={}, P={}", tempk, p_bar))?;

            alpha[t][p] = alpha_val;
            beta[t][p] = beta_val;

            p_bar += DELTA_P_BAR;
        }
        p_bar = P_BAR_MIN;
        tempk += DELTA_TEMPK;
    }

    write_output(&alpha, path, "Data/Crack_alpha.txt")?;
    write_output(&beta, path, "Data/Crack_beta.txt")?;

    if warnings {
        println!(
            "\n Outputs successfully generated in {}/Data/ directory:",
            path
        );
        println!("1. Crack_alpha.txt");
        println!("2. Crack_beta.txt");
    }

    Ok(())
}

/// Calculate log K for the dissolution of amorphous silica, chrysotile, and magnesite
/// over a range of T and P using CHNOSZ.
/// Outputs `Crack_silica.txt`, `Crack_chrysotile.txt`, and `Crack_magnesite.txt`.
pub fn crack_species_chnosz(path: &str, warnings: bool) -> Result<(), String> {
    let mut tempk = TEMPK_MIN_SPECIES;
    let mut p_bar = P_BAR_MIN;

    let mut silica = vec![vec![0.0; SIZEA_TP]; SIZEA_TP];
    let mut chrysotile = vec![vec![0.0; SIZEA_TP]; SIZEA_TP];
    let mut magnesite = vec![vec![0.0; SIZEA_TP]; SIZEA_TP];

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

    for t in 0..SIZEA_TP {
        if warnings {
            println!(
                "Crack_species_CHNOSZ: {:.1}% done",
                t as f64 / SIZEA_TP as f64 * 100.0
            );
        }
        let tempk_c = tempk - KELVIN;

        for p in 0..SIZEA_TP {
            // Silica dissolution reaction
            let silica_val = if tempk < 622.0 {
                let log_am_silica =
                    R!(r#"get_logK_safely("amorphous silica", "cr", {{tempk_c}}, {{p_bar}})"#)
                        .map_err(|e| e.to_string())?
                        .as_real()
                        .ok_or_else(|| {
                            format!(
                                "Failed to get logK for amorphous silica at T_c={}, P={}",
                                tempk_c, p_bar
                            )
                        })?;
                let log_sio2_aq = R!(r#"get_logK_safely("SiO2", "aq", {{tempk_c}}, {{p_bar}})"#)
                    .map_err(|e| e.to_string())?
                    .as_real()
                    .ok_or_else(|| {
                        format!(
                            "Failed to get logK for SiO2 aq at T_c={}, P={}",
                            tempk_c, p_bar
                        )
                    })?;
                -log_am_silica + log_sio2_aq
            } else {
                let log_quartz = R!(r#"get_logK_safely("quartz", "cr", {{tempk_c}}, {{p_bar}})"#)
                    .map_err(|e| e.to_string())?
                    .as_real()
                    .ok_or_else(|| {
                        format!(
                            "Failed to get logK for quartz at T_c={}, P={}",
                            tempk_c, p_bar
                        )
                    })?;
                let log_sio2_aq = R!(r#"get_logK_safely("SiO2", "aq", {{tempk_c}}, {{p_bar}})"#)
                    .map_err(|e| e.to_string())?
                    .as_real()
                    .ok_or_else(|| {
                        format!(
                            "Failed to get logK for SiO2 aq at T_c={}, P={}",
                            tempk_c, p_bar
                        )
                    })?;
                -log_quartz + log_sio2_aq
            };

            // Chrysotile dissolution reaction
            let log_chrysotile =
                R!(r#"get_logK_safely("chrysotile", "cr", {{tempk_c}}, {{p_bar}})"#)
                    .map_err(|e| e.to_string())?
                    .as_real()
                    .ok_or_else(|| {
                        format!(
                            "Failed to get logK for chrysotile at T_c={}, P={}",
                            tempk_c, p_bar
                        )
                    })?;
            let log_sio2_aq = R!(r#"get_logK_safely("SiO2", "aq", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!(
                        "Failed to get logK for SiO2 aq at T_c={}, P={}",
                        tempk_c, p_bar
                    )
                })?;
            let log_mg = R!(r#"get_logK_safely("Mg+2", "aq", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!(
                        "Failed to get logK for Mg+2 at T_c={}, P={}",
                        tempk_c, p_bar
                    )
                })?;
            let log_oh = R!(r#"get_logK_safely("OH-", "aq", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!("Failed to get logK for OH- at T_c={}, P={}", tempk_c, p_bar)
                })?;
            let log_h2o = R!(r#"get_logK_safely("H2O", "aq", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!(
                        "Failed to get logK for H2O aq at T_c={}, P={}",
                        tempk_c, p_bar
                    )
                })?;

            let chrysotile_val =
                -log_chrysotile + 2.0 * log_sio2_aq + 3.0 * log_mg + 6.0 * log_oh - 1.0 * log_h2o;

            // Magnesite dissolution reaction
            let log_magnesite = R!(r#"get_logK_safely("magnesite", "cr", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!(
                        "Failed to get logK for magnesite at T_c={}, P={}",
                        tempk_c, p_bar
                    )
                })?;
            let log_mg = R!(r#"get_logK_safely("Mg+2", "aq", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!(
                        "Failed to get logK for Mg+2 at T_c={}, P={}",
                        tempk_c, p_bar
                    )
                })?;
            let log_co3 = R!(r#"get_logK_safely("CO3-2", "aq", {{tempk_c}}, {{p_bar}})"#)
                .map_err(|e| e.to_string())?
                .as_real()
                .ok_or_else(|| {
                    format!(
                        "Failed to get logK for CO3-2 at T_c={}, P={}",
                        tempk_c, p_bar
                    )
                })?;

            let magnesite_val = -log_magnesite + log_mg + log_co3;

            silica[t][p] = silica_val;
            chrysotile[t][p] = chrysotile_val;
            magnesite[t][p] = magnesite_val;

            p_bar += DELTA_P_BAR;
        }
        p_bar = P_BAR_MIN;
        tempk += DELTA_TEMPK_SPECIES;
    }

    write_output(&silica, path, "Data/Crack_silica.txt")?;
    write_output(&chrysotile, path, "Data/Crack_chrysotile.txt")?;
    write_output(&magnesite, path, "Data/Crack_magnesite.txt")?;

    if warnings {
        println!(
            "\n Outputs successfully generated in {}/Data/ directory:",
            path
        );
        println!("1. Crack_silica.txt");
        println!("2. Crack_chrysotile.txt");
        println!("3. Crack_magnesite.txt\n");
    }

    Ok(())
}

#[cfg(test)]
/// These tests may take up a lot of CPU.
/// One of these tests failed. I'm not sure which one
mod crack_table_tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_a_tp() {
        // Passes
        let test_dir = "./target/test_crack_data_a_tp";
        let _ = fs::remove_dir_all(test_dir);

        let res = a_tp(test_dir, false);
        assert!(res.is_ok());

        assert!(Path::new(test_dir).join("Data/Crack_integral.txt").exists());
        assert!(Path::new(test_dir).join("Data/Crack_aTP.txt").exists());

        // Read and verify a few lines
        let integral_content =
            fs::read_to_string(Path::new(test_dir).join("Data/Crack_integral.txt")).unwrap();
        let lines: Vec<&str> = integral_content.lines().collect();
        assert_eq!(lines.len(), INT_SIZE as usize);

        let _ = fs::remove_dir_all(test_dir);
    }

    #[test]
    fn test_chnosz_tables() {
        // Since R/CHNOSZ is not thread-safe and must be initialized on a single thread,
        // we run all R-based table generation tests sequentially within a single test.
        let test_dir = "./target/test_crack_data_chnosz";
        let _ = fs::remove_dir_all(test_dir);

        // 1. Test water table generation
        // The R code successfully runs!
        let res_water = crack_water_chnosz(test_dir, false);
        assert!(res_water.is_ok());
        assert!(Path::new(test_dir).join("Data/Crack_alpha.txt").exists());
        assert!(Path::new(test_dir).join("Data/Crack_beta.txt").exists());

        // 2. Test species table generation
        let res_species = crack_species_chnosz(test_dir, false);
        assert!(res_species.is_ok());
        assert!(Path::new(test_dir).join("Data/Crack_silica.txt").exists());
        assert!(
            Path::new(test_dir)
                .join("Data/Crack_chrysotile.txt")
                .exists()
        );
        assert!(
            Path::new(test_dir)
                .join("Data/Crack_magnesite.txt")
                .exists()
        );

        let _ = fs::remove_dir_all(test_dir);
    }
}
