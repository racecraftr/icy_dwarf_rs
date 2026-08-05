//! This module simulates water-rock interactions and geochemical leaching using PHREEQC and CHNOSZ.

use crate::consts::{KELVIN, NAQ, NELTS, NGASES, NMINGAS, NVAR};
use crate::input::SubroutinesGeo;
use extendr_api::prelude::*;
use std::ffi::CString;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Explore the geochemistry of water-rock interactions using PHREEQC and CHNOSZ
/// across a specified parameter grid (temperature, pressure, pe, and water-rock ratio).
///
/// # Parameters
/// - `base_dir`: Base directory path for PHREEQC files (`PathBuf`).
/// - `output_dir`: Output directory path (`PathBuf`).
/// - `geo`: Geochemical parameter ranges.
pub fn param_exploration(
    base_dir: PathBuf,
    output_dir: PathBuf,
    geo: &SubroutinesGeo,
) -> Result<(), String> {
    let pe_min = 0.25 * geo.pe.min;
    let pe_max = 0.25 * geo.pe.max;
    let pe_step = 0.25 * geo.pe.step;

    let n_temp_iter = if geo.temp.step > 0.0 {
        ((geo.temp.max - geo.temp.min) / geo.temp.step).floor() as usize
    } else {
        0
    };

    let n_pressure_iter = if geo.pressure.step > 0.0 {
        ((geo.pressure.max - geo.pressure.min) / geo.pressure.step).floor() as usize
    } else {
        0
    };

    let n_pe_iter = if pe_step > 0.0 {
        ((pe_max - pe_min) / pe_step).floor() as usize
    } else {
        0
    };

    let n_wr_iter = if geo.wr_ratio.step > 1.0 {
        ((geo.wr_ratio.max.ln() - geo.wr_ratio.min.ln()) / geo.wr_ratio.step.ln()).ceil() as usize
    } else {
        0
    };

    let phreeqc_dir = if base_dir.join("Phreeqc").exists() {
        base_dir.join("Phreeqc")
    } else if base_dir.join("PHREEQC-3.1.2").exists() {
        base_dir.join("PHREEQC-3.1.2")
    } else if PathBuf::from("Phreeqc").exists() {
        PathBuf::from("Phreeqc")
    } else {
        PathBuf::from("PHREEQC-3.1.2")
    };

    let dbase = phreeqc_dir.join("core10.dat");
    let infile = phreeqc_dir.join("io").join("PHREEQCinput");
    let solfile = phreeqc_dir.join("io").join("Sol");

    let outfile = output_dir.join("ParamExploration.txt");
    create_output(&outfile)?;

    init_chnosz()?;

    let mut simdata = vec![vec![0.0; NVAR as usize]; n_pe_iter + 1];

    for i_pressure in 0..=n_pressure_iter {
        let mut p = geo.pressure.min + geo.pressure.step * (i_pressure as f64);

        for i_temp in 0..=n_temp_iter {
            let mut t = geo.temp.min + geo.temp.step * (i_temp as f64);

            if t == 0.0 {
                t = 0.01; // PHREEQC crashes at 0 celsius
            }
            if p == 0.0 {
                p = 1.0; // 1 bar minimum
            }

            let log_quartz = get_chnosz_log_k("quartz", "cr", t, p)?;
            let log_magnetite = get_chnosz_log_k("magnetite", "cr", t, p)?;
            let log_fayalite = get_chnosz_log_k("fayalite", "cr", t, p)?;
            let log_o2 = get_chnosz_log_k("O2", "g", t, p)?;

            let log_h_plus = get_chnosz_log_k("H+", "aq", t, p)?;
            let log_e_minus = get_chnosz_log_k("e-", "aq", t, p)?;
            let log_h2o = get_chnosz_log_k("H2O", "liq", t, p)?;

            let logf_o2 = -3.0 * log_quartz - 2.0 * log_magnetite + 3.0 * log_fayalite + log_o2;
            let log_ko2_h2o = -4.0 * log_h_plus - 4.0 * log_e_minus - log_o2 + 2.0 * log_h2o;

            for i_wr in 0..=n_wr_iter {
                let wr = geo.wr_ratio.max / geo.wr_ratio.step.powi(i_wr as i32);

                println!(
                    "P={} ({} of {}) T={} ({} of {}) W:R={} ({} of {}), parallel calculations over {} values of pe",
                    p,
                    i_pressure + 1,
                    n_pressure_iter + 1,
                    t,
                    i_temp + 1,
                    n_temp_iter + 1,
                    wr,
                    i_wr + 1,
                    n_wr_iter + 1,
                    n_pe_iter + 1
                );

                for row in simdata.iter_mut() {
                    row.fill(0.0);
                }

                let tempinput = phreeqc_dir.join("io").join("PHREEQCinput_temp");

                for ipe in 0..=n_pe_iter {
                    let pe = pe_min + pe_step * (ipe as f64);

                    let mut ph = 7.0; // Neutral PH
                    let mut fmq = -ph + 0.25 * (logf_o2 + log_ko2_h2o);
                    println!("FMQ pe is {} at T={} C, P={} bar, and pH {}", fmq, t, p, ph);

                    write_phreeqc_input(&solfile, t, p, ph, 4.0 * pe, fmq + pe, wr, &tempinput)?;

                    let instance = unsafe { CreateIPhreeqc() };
                    if instance < 0 {
                        return Err("Failed to create IPhreeqc instance".to_string());
                    }

                    let c_dbase = CString::new(dbase.to_string_lossy().as_bytes())
                        .map_err(|e| e.to_string())?;
                    let c_tempinput = CString::new(tempinput.to_string_lossy().as_bytes())
                        .map_err(|e| e.to_string())?;

                    unsafe {
                        if LoadDatabase(instance, c_dbase.as_ptr()) != 0 {
                            OutputErrorString(instance);
                            DestroyIPhreeqc(instance);
                            return Err("Failed to load PHREEQC database".to_string());
                        }
                        SetSelectedOutputFileOn(instance, 1);
                        if RunFile(instance, c_tempinput.as_ptr()) != 0 {
                            OutputErrorString(instance);
                            DestroyIPhreeqc(instance);
                            return Err("Failed to run PHREEQC solution file".to_string());
                        }
                    }

                    extract_ph(instance, &mut ph);
                    extract_write_sol(instance, &mut simdata[ipe]);
                    unsafe { DestroyIPhreeqc(instance) };

                    fmq = -ph + 0.25 * (logf_o2 + log_ko2_h2o);
                    println!("FMQ pe is {} at T={} C, P={} bar, and pH {}", fmq, t, p, ph);

                    write_phreeqc_input(&infile, t, p, ph, 4.0 * pe, fmq + pe, wr, &tempinput)?;

                    let instance2 = unsafe { CreateIPhreeqc() };
                    if instance2 < 0 {
                        return Err("Failed to create IPhreeqc instance".to_string());
                    }

                    unsafe {
                        if LoadDatabase(instance2, c_dbase.as_ptr()) != 0 {
                            OutputErrorString(instance2);
                            DestroyIPhreeqc(instance2);
                            return Err("Failed to load PHREEQC database".to_string());
                        }
                        SetSelectedOutputFileOn(instance2, 1);
                        if RunFile(instance2, c_tempinput.as_ptr()) != 0 {
                            OutputErrorString(instance2);
                            DestroyIPhreeqc(instance2);
                            return Err("Failed to run PHREEQC input file".to_string());
                        }
                    }

                    simdata[ipe][1] = p;
                    simdata[ipe][3] = pe * 4.0;
                    simdata[ipe][4] = fmq;
                    extract_write(instance2, &mut simdata[ipe]);

                    unsafe { DestroyIPhreeqc(instance2) };
                }

                for row in &simdata {
                    append_output(&outfile, row)?;
                }
            }
        }
    }

    Ok(())
}

/// Simulate single water-rock interaction step using PHREEQC and CHNOSZ.
#[allow(dead_code)]
pub fn water_rock(base_dir: &Path, t: f64, p: f64, pe: f64, mut wr: f64) -> Result<f64, String> {
    let ph = 7.0;

    let phreeqc_dir = if base_dir.join("Phreeqc").exists() {
        base_dir.join("Phreeqc")
    } else if base_dir.join("PHREEQC-3.1.2").exists() {
        base_dir.join("PHREEQC-3.1.2")
    } else if PathBuf::from("Phreeqc").exists() {
        PathBuf::from("Phreeqc")
    } else {
        PathBuf::from("PHREEQC-3.1.2")
    };

    let dbase = phreeqc_dir.join("core10.dat");
    let infile = phreeqc_dir.join("io").join("PHREEQCinput");
    let tempinput = phreeqc_dir.join("io").join("PHREEQCinput_temp");

    let molmass = load_mol_mass(base_dir)?;

    let t = t - KELVIN;

    init_chnosz()?;

    let log_quartz = get_chnosz_log_k("quartz", "cr", t, p)?;
    let log_magnetite = get_chnosz_log_k("magnetite", "cr", t, p)?;
    let log_fayalite = get_chnosz_log_k("fayalite", "cr", t, p)?;
    let log_o2 = get_chnosz_log_k("O2", "g", t, p)?;

    let log_h_plus = get_chnosz_log_k("H+", "aq", t, p)?;
    let log_e_minus = get_chnosz_log_k("e-", "aq", t, p)?;
    let log_h2o = get_chnosz_log_k("H2O", "liq", t, p)?;

    let logf_o2 = -3.0 * log_quartz - 2.0 * log_magnetite + 3.0 * log_fayalite + log_o2;
    let log_ko2_h2o = -4.0 * log_h_plus - 4.0 * log_e_minus - log_o2 + 2.0 * log_h2o;

    let fmq = -ph + 0.25 * (logf_o2 + log_ko2_h2o) + pe;

    if wr < 0.5 {
        println!("WR is {} < 0.5, assuming WR=0.5", wr);
        wr = 0.5;
    }

    write_phreeqc_input(&infile, t, p, ph, fmq, wr, wr, &tempinput)?;

    let instance = unsafe { CreateIPhreeqc() };
    if instance < 0 {
        return Err("Failed to create IPhreeqc instance".to_string());
    }

    let c_dbase = CString::new(dbase.to_string_lossy().as_bytes()).map_err(|e| e.to_string())?;
    let c_tempinput =
        CString::new(tempinput.to_string_lossy().as_bytes()).map_err(|e| e.to_string())?;

    if unsafe { LoadDatabase(instance, c_dbase.as_ptr()) } != 0 {
        unsafe { OutputErrorString(instance) };
        unsafe { DestroyIPhreeqc(instance) };
        return Err("Failed to load PHREEQC database".to_string());
    }

    unsafe { SetSelectedOutputFileOn(instance, 1) };

    if unsafe { RunFile(instance, c_tempinput.as_ptr()) } != 0 {
        unsafe { OutputErrorString(instance) };
        unsafe { DestroyIPhreeqc(instance) };
        return Err("Failed to run PHREEQC input file".to_string());
    }

    let mut simdata = vec![0.0; NVAR as usize];
    extract_write_sol(instance, &mut simdata);
    extract_write(instance, &mut simdata);

    unsafe { DestroyIPhreeqc(instance) };

    let mass_water = simdata[11];
    let dissolved_k = simdata[23];

    let get_molmass = |idx: usize, col: usize| -> f64 {
        if idx < molmass.len() && col < molmass[idx].len() {
            molmass[idx][col]
        } else {
            0.0
        }
    };

    let total_k = (simdata[1488] - simdata[1489]) * get_molmass(1488, 11)
        + (simdata[1356] - simdata[1357]) * get_molmass(1356, 11)
        + (simdata[1226] - simdata[1227]) * get_molmass(1226, 11)
        + (simdata[748] - simdata[749]) * get_molmass(748, 11);

    let frac_k_leached = if total_k > 0.0 {
        dissolved_k * mass_water / total_k
    } else {
        0.0
    };

    Ok(frac_k_leached)
}

fn init_chnosz() -> Result<(), String> {
    extendr_engine::start_r();
    let _ = R!(r#"
        if (!requireNamespace("CHNOSZ", quietly = TRUE)) {
            stop("CHNOSZ package not installed in R")
        }
        library(CHNOSZ, quietly = TRUE)
        data(thermo)
        get_logK <- function(species, state, T, P) {
            res <- subcrt(species, state, T = T, P = P)
            res$out[[1]]$logK[1]
        }
    "#)
    .map_err(|e| format!("Failed to initialize CHNOSZ in R: {:?}", e))?;
    Ok(())
}

fn get_chnosz_log_k(species: &str, state: &str, temp: f64, press: f64) -> Result<f64, String> {
    let val = R!(r#"get_logK({{species}}, {{state}}, {{temp}}, {{press}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| format!("Failed to get logK for species {}", species))?;
    Ok(val)
}

/// Variant types for IPhreeqc C FFI data exchange.
#[repr(C)]
#[allow(dead_code)]
pub enum VarType {
    /// Empty variable type.
    Empty = 0,
    /// Error variable type.
    Error = 1,
    /// Long integer variable type.
    Long = 2,
    /// Double precision floating-point variable type.
    Double = 3,
    /// String pointer variable type.
    String = 4,
}

/// Raw variable values for IPhreeqc C FFI data exchange.
#[repr(C)]
pub union VarValue {
    /// Long integer value.
    pub l_val: std::ffi::c_long,
    /// Double precision floating-point value.
    pub d_val: f64,
    /// Character pointer value.
    pub s_val: *mut std::ffi::c_char,
    /// Result code integer value.
    pub vresult: i32,
}

/// Variable wrapper for IPhreeqc C FFI data exchange.
#[repr(C)]
pub struct Var {
    /// The type of the variable.
    pub vtype: VarType,
    /// The union value of the variable.
    pub val: VarValue,
}

#[link(name = "iphreeqc")]
unsafe extern "C" {
    fn CreateIPhreeqc() -> std::ffi::c_int;
    fn DestroyIPhreeqc(id: std::ffi::c_int) -> std::ffi::c_int;
    fn LoadDatabase(id: std::ffi::c_int, filename: *const std::ffi::c_char) -> std::ffi::c_int;
    fn SetSelectedOutputFileOn(id: std::ffi::c_int, status: std::ffi::c_int) -> std::ffi::c_int;
    fn RunFile(id: std::ffi::c_int, filename: *const std::ffi::c_char) -> std::ffi::c_int;
    fn VarInit(pvar: *mut Var);
    fn GetSelectedOutputValue(
        id: std::ffi::c_int,
        row: std::ffi::c_int,
        col: std::ffi::c_int,
        pvar: *mut Var,
    ) -> std::ffi::c_int;
    fn OutputErrorString(id: std::ffi::c_int);
}

fn extract_ph(instance: i32, ph: &mut f64) {
    let mut v = Var {
        vtype: VarType::Empty,
        val: VarValue { l_val: 0 },
    };
    unsafe {
        VarInit(&mut v);
        if GetSelectedOutputValue(instance, 1, 1, &mut v) == 0 && matches!(v.vtype, VarType::Double)
        {
            *ph = v.val.d_val;
        }
    }
}

fn extract_write_sol(instance: i32, data: &mut [f64]) {
    for i in 2..30 {
        let mut v = Var {
            vtype: VarType::Empty,
            val: VarValue { l_val: 0 },
        };
        unsafe {
            VarInit(&mut v);
            if GetSelectedOutputValue(instance, 1, i, &mut v) == 0 {
                let val = match v.vtype {
                    VarType::Double => v.val.d_val,
                    VarType::Long => v.val.l_val as f64,
                    _ => 0.0,
                };
                let idx = (i + 5) as usize;
                if idx < data.len() {
                    data[idx] = if val.abs() < 1.0e-50 { 0.0 } else { val };
                }
            }
        }
    }
}

fn extract_write(instance: i32, data: &mut [f64]) {
    let get_val = |row: i32, col: i32| -> f64 {
        let mut v = Var {
            vtype: VarType::Empty,
            val: VarValue { l_val: 0 },
        };
        unsafe {
            VarInit(&mut v);
            if GetSelectedOutputValue(instance, row, col, &mut v) == 0 {
                match v.vtype {
                    VarType::Double => v.val.d_val,
                    VarType::Long => v.val.l_val as f64,
                    _ => 0.0,
                }
            } else {
                0.0
            }
        }
    };

    if data.len() > 36 {
        data[0] = get_val(1, 3);
        data[2] = get_val(1, 5);
        data[5] = get_val(1, 1);
        data[6] = get_val(1, 2);
        data[35] = get_val(2, 1);
        data[36] = get_val(2, 2);
    }

    let limit = NVAR - 6 - 27;
    for i in 4..limit {
        let idx = (i + 6 + 27) as usize;
        if idx >= data.len() {
            break;
        }
        let val = get_val(2, i);
        data[idx] = if val.abs() < 1.0e-50 { 0.0 } else { val };
    }
}

/// Load molar mass vectors from the data file.
///
/// # Parameters
/// - `data_dir`: The directory path (`&Path`) pointing to the data folder.
///
/// # Returns
/// A nested vector containing molar masses of chemical species.
#[allow(dead_code)]
pub fn load_mol_mass(data_dir: &Path) -> Result<Vec<Vec<f64>>, String> {
    let mut molmass = vec![vec![0.0; NELTS as usize]; NVAR as usize];
    let file_path = data_dir.join("Molar_masses.txt");

    let Ok(content) = fs::read_to_string(&file_path) else {
        return Err(format!("Could not read {}", file_path.display()));
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

    let ngases = NGASES as usize;
    let nmingas = NMINGAS as usize;
    let naq = NAQ as usize;

    for i in 0..ngases {
        for j in 0..(NELTS as usize) {
            molmass[naq + 2 * (nmingas - ngases) + 5 - 1 + i][j] =
                read_data[nmingas - ngases + i][j];
        }
    }

    let mut k = naq - 1;
    for datum in read_data.iter().take(nmingas - ngases) {
        for j in 0..(NELTS as usize) {
            molmass[k][j] = datum[j];
            molmass[k + 1][j] = molmass[k][j];
        }
        k += 2;
    }

    for j in 0..(NELTS as usize) {
        molmass[0][j] = read_data[0][j];
    }

    Ok(molmass)
}

fn write_phreeqc_input(
    template_file: &Path,
    temp: f64,
    pressure: f64,
    ph: f64,
    _rel_pe: f64,
    pe: f64,
    wr: f64,
    output_file: &Path,
) -> Result<(), String> {
    let content = fs::read_to_string(template_file).map_err(|e| {
        format!(
            "Could not read template file {}: {}",
            template_file.display(),
            e
        )
    })?;

    let mut output = String::with_capacity(content.len() + 128);
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
        } else if line.trim_start().starts_with("-pres") {
            output.push_str(&format!("\t -pressure \t{}\n", pressure));
        } else if line.trim_start().starts_with("-temp") {
            output.push_str(&format!("\t -temperature \t{}\n", temp));
        } else {
            output.push_str(line);
            output.push('\n');
        }
    }

    fs::write(output_file, output).map_err(|e| {
        format!(
            "Could not write PHREEQC input file {}: {}",
            output_file.display(),
            e
        )
    })?;
    Ok(())
}

/// Create an output file in the output directory if it does not exist.
pub fn create_output(file_path: &Path) -> Result<(), String> {
    if let Some(parent) = file_path.parent()
        && let Err(e) = fs::create_dir_all(parent)
    {
        return Err(format!(
            "Unable to create folder {}: {}",
            parent.display(),
            e
        ));
    }
    if let Err(e) = File::create(file_path) {
        return Err(format!(
            "Unable to create file {}: {}",
            file_path.display(),
            e
        ));
    }
    Ok(())
}

/// Append a row of tab-delimited numerical data to an output file.
pub fn append_output(file_path: &Path, data: &[f64]) -> Result<(), String> {
    let mut file = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_path)
        .map_err(|e| format!("Unable to open file {}: {}", file_path.display(), e))?;

    let mut line = String::with_capacity(data.len() * 12);
    for (i, val) in data.iter().enumerate() {
        if i > 0 {
            line.push('\t');
        }
        line.push_str(&format!("{}", val));
    }
    line.push('\n');

    file.write_all(line.as_bytes())
        .map_err(|e| format!("Unable to write to file {}: {}", file_path.display(), e))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_r_works() {
        test! {
            let n = R!(r#"2 + 2"#).unwrap().as_real().unwrap();
            assert!(n == 4.);
        }
    }

    #[test]
    fn test_iphreeqc_create_and_destroy() {
        let instance = unsafe { CreateIPhreeqc() };
        assert!(
            instance >= 0,
            "CreateIPhreeqc should return a non-negative instance ID"
        );
        let res = unsafe { DestroyIPhreeqc(instance) };
        assert_eq!(res, 0, "DestroyIPhreeqc should return 0 (IPQ_OK)");
    }

    #[test]
    fn test_iphreeqc_var_init() {
        let mut v = Var {
            vtype: VarType::Error,
            val: VarValue { l_val: 42 },
        };
        unsafe { VarInit(&mut v) };
        assert!(matches!(v.vtype, VarType::Empty));
    }
}
