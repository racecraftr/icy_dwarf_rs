use crate::consts::{KELVIN, NAQ, NELTS, NGASES, NMINGAS, NVAR};
use extendr_api::prelude::*;
use extendr_engine::*;
use std::fs;

/// Simulates water-rock interactions using PHREEQC.
/// Returns the fraction of K (potassium) leached.
pub fn water_rock(path: &str, t: f64, p: f64, mut wr: f64, chondrite: i32) -> Result<f64, String> {
    let mut ph = 7.0;

    let dbase = format!("{}/PHREEQC-3.1.2/core9.dat", path);
    let infile = format!("{}/io/PHREEQCinput", path);
    let tempinput = format!("{}/io/PHREEQCinput_temp", path);

    let molmass = load_mol_mass(path)?;

    // To call R code from Rust, the recommended approach is to use the `extendr-api` crate.
    // Example (requires `extendr-api` in Cargo.toml):
    // extendr_api::extendr_engine::start_r();
    // let logfO2 = R!(r#"
    //     -3.0 * CHNOSZ::subcrt(c("quartz"), "cr", T=t_c, P=p_bar)$out$logK + ...
    // "#)
    // .unwrap()
    // .as_real()
    // .unwrap();
    extendr_engine::start_r();
    let log_o2 = R!(r#"
            -3.0 * CHNOSZ::subcrt(c("quartz"), "cr", T={{t}}, P={{p}}, )
        "#);

    // Placeholder for R/CHNOSZ calculations:
    let logf_oxygen = 0.0; // TODO: Replace with extendr-api call to CHNOSZ
    let log_ko2_h2o = 0.0; // TODO: Replace with extendr-api call to CHNOSZ

    let mut fmq = -ph + 0.25 * (logf_oxygen + log_ko2_h2o);

    if wr < 0.5 {
        println!("WR is {} < 0.5, assuming WR=0.5", wr);
        wr = 0.5;
    }

    write_phreeqc_input(&infile, t, p, ph, 0.0, fmq, wr, &tempinput)?;

    // Calling IPhreeqc requires either a `phreeqc-sys` crate or manually mapping the C FFI.
    // Example of how the FFI calls would look in Rust:
    /*
    unsafe {
        let id = CreateIPhreeqc();
        let db_c = std::ffi::CString::new(dbase).unwrap();
        LoadDatabase(id, db_c.as_ptr());
        SetSelectedOutputFileOn(id, 1);

        let input_c = std::ffi::CString::new(tempinput).unwrap();
        RunFile(id, input_c.as_ptr());

        // Extract results...

        DestroyIPhreeqc(id);
    }
    */

    let mass_water = 0.0; // TODO: Extract from IPhreeqc
    let total_k = 0.0; // TODO: Calculate based on chondrite type and extracted molmass
    let dissolved_k = 0.0; // TODO: Extract from IPhreeqc

    let frac_k_leached = if total_k > 0.0 {
        dissolved_k * mass_water / total_k
    } else {
        0.0
    };

    Ok(frac_k_leached)
}

/// Loads molar masses from Data/Molar_masses.txt
pub fn load_mol_mass(path: &str) -> Result<Vec<Vec<f64>>, String> {
    let mut molmass = vec![vec![0.0; NELTS as usize]; NVAR as usize];
    let file_path = format!("{}/Data/Molar_masses.txt", path);

    let Ok(content) = fs::read_to_string(&file_path) else {
        return Err(format!("Could not read {}", file_path));
    };

    let mut read_data = vec![];
    for line in content.lines() {
        let nums: Vec<f64> = line
            .split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();
        if nums.len() == NELTS as usize {
            read_data.push(nums);
        }
    }

    if read_data.is_empty() {
        return Err("Molar_masses.txt was empty or invalid".to_string());
    }

    // Shift to positions corresponding to simdata
    // Gas species
    let ngases = NGASES as usize;
    let nmingas = NMINGAS as usize;
    let naq = NAQ as usize;

    for i in 0..ngases {
        for j in 0..(NELTS as usize) {
            molmass[naq + 2 * (nmingas - ngases) + 5 - 1 + i][j] =
                read_data[nmingas - ngases + i][j];
        }
    }

    // Solid species
    let mut k = naq - 1;
    for datum in read_data.iter().take(nmingas - ngases) {
        for j in 0..(NELTS as usize) {
            molmass[k][j] = datum[j];
            molmass[k + 1][j] = molmass[k][j];
        }
        k += 2
    }

    // First line
    for j in 0..(NELTS as usize) {
        molmass[0][j] = read_data[0][j];
    }

    Ok(molmass)
}

fn write_phreeqc_input(
    template_file: &str,
    temp: f64,
    pressure: f64,
    ph: f64,
    rel_pe: f64,
    pe: f64,
    wr: f64,
    output_file: &str,
) -> Result<(), String> {
    let Ok(content) = fs::read_to_string(template_file) else {
        return Err(format!("Could not read template file {}", template_file));
    };

    let mut output = String::new();
    for (i, line) in content.lines().enumerate() {
        let line_no = i + 1;
        if line_no == 5 {
            output.push_str(&format!("\t pH \t \t{}\t charge\n", ph));
        } else if line_no == 6 {
            output.push_str(&format!("\t temp \t \t{}\n", temp));
        } else if line_no == 7 {
            output.push_str(&format!("\t pressure \t{}\n", pressure));
        } else if line_no == 8 {
            output.push_str(&format!("\t pe \t \t{}\n", pe));
        } else if line_no == 9 {
            output.push_str(&format!("\t -water \t{}\n", wr));
        } else if line.starts_with("-pres") {
            output.push_str(&format!("\t -pressure \t{}\n", pressure));
        } else if line.starts_with("-temp") {
            output.push_str(&format!("\t -temperature \t{}\n", temp));
        } else {
            output.push_str(line);
            output.push('\n');
        }
    }

    fs::write(output_file, output).map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod r_tests {
    use super::*;
    #[test]
    fn test_r_works() {
        test! {
            let n = R!(r#"2 + 2"#).unwrap().as_real().unwrap();
            assert!(n == 4.);
        }
    }
}
