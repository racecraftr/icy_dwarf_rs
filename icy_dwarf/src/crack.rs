use std::{fs, path::Path};

use crate::consts::SIZEA_TP;

pub struct Data {
    pub atp: Vec<Vec<f64>>,
    pub alpha: Vec<Vec<f64>>,
    pub beta: Vec<Vec<f64>>,
    pub chrysotile: Vec<Vec<f64>>,
    pub integral: Vec<Vec<f64>>,
    pub magnesite: Vec<Vec<f64>>,
    pub silica: Vec<Vec<f64>>,
}

pub struct CompPlanmat {
    pub material: String,
    pub conditions: String,
    pub reference: String,
    pub index: usize,
    pub eos: f64,
    pub rho0: f64,
    pub c: f64,
    pub n: f64,
    pub k0: f64,
    pub dk_dp: f64,
    pub v0: f64,
    pub tref: f64,
    pub a0: f64,
    pub a1: f64,
    pub b1: f64,
    pub b2: f64,
}

pub fn read_data(dat_folder: &String) -> Option<Data> {
    Some(Data {
        atp: read_data_file(&format!("{}/Crack_aTP.txt", dat_folder))?,
        alpha: read_data_file(&format!("{}/Crack_alpha.txt", dat_folder))?,
        beta: read_data_file(&format!("{}/Crack_bet.txt", dat_folder))?,
        chrysotile: read_data_file(&format!("{}/Crack_chrysotile.txt", dat_folder))?,
        integral: read_data_file(&format!("{}/Crack_integral.txt", dat_folder))?,
        magnesite: read_data_file(&format!("{}/Crack_magnesite.txt", dat_folder))?,
        silica: read_data_file(&format!("{}/Crack_magnesite.txt", dat_folder))?,
    })
}

/// Reads a file in the specified Data/ folder.
fn read_data_file(dat_file_path: &String) -> Option<Vec<Vec<f64>>> {
    let Ok(s) = fs::read_to_string(Path::new(dat_file_path)) else {
        return None;
    };
    let lines = s.lines().collect::<Vec<_>>();
    if lines.len() == SIZEA_TP {
        let res = lines
            .iter()
            .map(|s| {
                s.split_whitespace()
                    .filter_map(|t| t.parse::<f64>().ok())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        if res.iter().all(|v| v.len() == SIZEA_TP) {
            Some(res)
        } else {
            None
        }
    } else {
        None
    }
}
