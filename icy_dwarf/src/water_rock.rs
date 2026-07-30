use crate::consts::{KELVIN, NAQ, NELTS, NGASES, NMINGAS, NVAR};
use extendr_api::prelude::*;
use std::fs;

/// Simulates water-rock interactions using PHREEQC.
/// Returns the fraction of K (potassium) leached.
pub fn water_rock(path: &str, t: f64, p: f64, pe: f64, mut wr: f64) -> Result<f64, String> {
    let ph = 7.0;

    let dbase = format!("{}/PHREEQC-3.1.2/core9.dat", path);
    let infile = format!("{}/io/PHREEQCinput", path);
    let tempinput = format!("{}/io/PHREEQCinput_temp", path);

    let molmass = load_mol_mass(path)?;

    let t = t - KELVIN;

    extendr_engine::start_r();
    let _ = R!(r#"
        if (!requireNamespace("CHNOSZ", quietly = TRUE)) {
            stop("CHNOSZ package not installed in R")
        }
        library(CHNOSZ, quietly = TRUE)
        data(thermo)
        get_log <- function(species, state, T, P) {
            res <- subcrt(species, state, T = T, P = P)
            res$out[[1]]$logK[1]
        }
    "#)
    .map_err(|e| format!("Failed to initialize CHNOSZ in R: {:?}", e))?;

    let log_quartz = R!(r#"get_logK("quartz", "cr", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for quartz".to_string())?;
    let log_magnetite = R!(r#"get_logK("magnetite", "cr", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for magnetite".to_string())?;
    let log_fayalite = R!(r#"get_logK("fayalite", "cr", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for fayalite".to_string())?;
    let log_o2 = R!(r#"get_logK("O2", "g", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for O2".to_string())?;

    let log_h_plus = R!(r#"get_logK("H+", "aq", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for H+".to_string())?;
    let log_e_minus = R!(r#"get_logK("e-", "aq", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for e-".to_string())?;
    let log_h2o = R!(r#"get_logK("H2O", "liq", {{t}}, {{p}})"#)
        .map_err(|e| e.to_string())?
        .as_real()
        .ok_or_else(|| "Failed to get logK for H2O".to_string())?;

    let logf_o2 = -3.0 * log_quartz - 2.0 * log_magnetite + 3.0 * log_fayalite + 1.0 * log_o2;
    let log_ko2_h2o = -4.0 * log_h_plus - 4.0 * log_e_minus - 1.0 * log_o2 + 2.0 * log_h2o;

    let fmq = -ph + 0.25 * (logf_o2 + log_ko2_h2o) + pe;

    if wr < 0.5 {
        println!("WR is {} < 0.5, assuming WR=0.5", wr);
        wr = 0.5;
    }

    write_phreeqc_input(&infile, t, p, ph, fmq, wr, &tempinput)?;

    let instance = unsafe { CreateIPhreeqc() };
    if instance < 0 {
        return Err("Failed to create IPhreeqc instance".to_string());
    }

    let c_dbase = std::ffi::CString::new(dbase).map_err(|e| e.to_string())?;
    let c_tempinput = std::ffi::CString::new(tempinput).map_err(|e| e.to_string())?;

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

    let mut simdata = vec![0.0; 2000];
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

#[repr(C)]
pub enum VarType {
    Empty = 0,
    Error = 1,
    Long = 2,
    Double = 3,
    String = 4,
}

#[repr(C)]
pub union VarValue {
    pub l_val: std::ffi::c_long,
    pub d_val: f64,
    pub s_val: *mut std::ffi::c_char,
    pub vresult: i32,
}

#[repr(C)]
pub struct Var {
    pub vtype: VarType,
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

    let limit = data.len();
    for i in 4.. {
        let idx = (i + 6 + 27) as usize;
        if idx >= limit {
            break;
        }
        let val = get_val(2, i);
        data[idx] = if val.abs() < 1.0e-50 { 0.0 } else { val };
    }
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
    let nmingas = NMINGAS as usize; // PINGAS's little-known cousin
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
    // rel_pe: f64,
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
        assert!(instance >= 0, "CreateIPhreeqc should return a non-negative instance ID");
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
