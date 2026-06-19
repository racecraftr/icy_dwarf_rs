use std::{fs, path::PathBuf};

use itertools::Itertools;

use crate::{input::IcyDwarfInput, thermal::ThermalOut};

impl IcyDwarfInput {
    pub fn compression(&self, thermal_out: &[Vec<ThermalOut>]) {}
}

struct PlanMatEntry {
    db_index: i32,
    eos: i32,
    rho_0: i32,
    c: i32,
    nn: f64,
    ks_0: f64,
    ks_p: f64,
    v_0: f64,
    t_ref: f64,
    a: [f64; 2],
    b: [f64; 3],
}

impl PlanMatEntry {
    pub fn from(entry: &[&str]) -> Option<Self> {
        if entry.len() != 14 {
            None
        } else {
            Some(Self {
                db_index: entry[0].parse().ok()?,
                eos: entry[1].parse().ok()?,
                rho_0: entry[2].parse().ok()?,
                c: entry[3].parse().ok()?,
                nn: entry[4].parse().ok()?,
                ks_0: entry[5].parse().ok()?,
                ks_p: entry[6].parse().ok()?,
                v_0: entry[7].parse().ok()?,
                t_ref: entry[8].parse().ok()?,
                a: [entry[9].parse().ok()?, entry[10].parse().ok()?],
                b: [
                    entry[11].parse().ok()?,
                    entry[12].parse().ok()?,
                    entry[13].parse().ok()?,
                ],
            })
        }
    }
}

pub fn planmat(data_folder: &str, n_comp: usize) -> Option<Vec<PlanMatEntry>> {
    let planmat_db = PathBuf::from(data_folder).join("Compression_planmat.txt");
    let planmat_db = fs::read_to_string(planmat_db).ok()?;
    Some(
        planmat_db
            .split_whitespace()
            .chunks(14)
            .into_iter()
            .take(n_comp)
            .filter_map(|chunk| PlanMatEntry::from(&chunk.collect::<Vec<_>>()))
            .collect(),
    )
}
