use std::{f64::consts::PI, fs, io::Write, path::PathBuf};

use itertools::Itertools;
use num::traits::Inv;

use crate::{
    consts::{G, GRAM, KM, KM2CM, M_EARTH, MPA, PI_GREEK, R_EARTH},
    input::Fracs,
    planet_system::WorldState,
    thermal::ThermalOut,
};

const FRAC_3_4PI: f64 = 0.75 / PI;

impl WorldState {
    pub fn compression(
        &self,
        data_folder: &str,
        indices: (usize, usize, usize),
        thermal_outputs: &[Vec<ThermalOut>],
        time_step: usize,
        specify: bool,
        rho_hydr: f64,
        rho_dry: f64,
    ) -> Option<()> {
        let planmat_db = planmat(data_folder, 50)?;

        let iincore = planmat_index(&planmat_db, indices.0)?;
        let ioutcore = planmat_index(&planmat_db, indices.1)?;
        let imantle = planmat_index(&planmat_db, indices.2)?;

        let nr = self.zones.len();

        let mut mp = 0.0;
        let mut m_incore = 0.0;
        let mut m_outcore = 0.0;

        for v in thermal_outputs {
            let out = &v[time_step];
            mp += out.mass_total() * GRAM;
            if out.deg_of_hydr < 0.1 {
                m_incore += out.mass_rock * GRAM;
            } else {
                m_outcore += out.mass_rock * GRAM;
            }
        }
        let m_mantle = mp - m_incore - m_outcore;

        let vol_incore = m_incore / planmat_db[iincore].rho_0;
        let r_incore = (FRAC_3_4PI * vol_incore).powf(1.0 / 3.0);

        let vol_outcore = m_outcore / planmat_db[ioutcore].rho_0;
        let r_outcore = (FRAC_3_4PI * vol_outcore + r_incore.powi(3)).powf(1.0 / 3.0);

        let vol_mantle = m_mantle / planmat_db[imantle].rho_0;
        let rp = (FRAC_3_4PI * vol_mantle + r_outcore.powi(3)).powf(1.0 / 3.0);

        let nic = (r_incore / rp * nr as f64).floor() as usize;
        let noc = (r_outcore / rp * nr as f64).floor() as usize;

        let mut icomp = vec![0; nr + 1];
        let mut r = vec![0.0; nr + 1];
        let mut rho = vec![0.0; nr + 1];

        if nic > 0 {
            for ir in 0..=nic {
                icomp[ir] = iincore;
                r[ir] = r_incore * (ir as f64) / (nic as f64);
                rho[ir] = planmat_db[icomp[ir]].rho_0;
            }
        } else {
            icomp[0] = iincore;
            r[0] = 0.0;
            rho[0] = planmat_db[icomp[0]].rho_0;
        }

        if noc > nic {
            for ir in (nic + 1)..=noc {
                icomp[ir] = ioutcore;
                r[ir] =
                    r_incore + (r_outcore - r_incore) * ((ir - nic) as f64) / ((noc - nic) as f64);
                rho[ir] = planmat_db[icomp[ir]].rho_0;
            }
        } else if noc == 0 {
            icomp[0] = ioutcore;
            r[0] = 0.0;
            rho[0] = planmat_db[icomp[0]].rho_0;
        }

        for ir in (noc + 1)..=nr {
            icomp[ir] = imantle;
            r[ir] = r_outcore + (rp - r_outcore) * ((ir - noc) as f64) / ((nr - noc) as f64);
            rho[ir] = planmat_db[icomp[ir]].rho_0;
        }

        let mut d_m = vec![0.0; nr + 1];
        let mut m = vec![0.0; nr + 1];
        for ir in 1..=nr {
            let vol = (PI / 0.75) * (r[ir].powi(3) - r[ir - 1].powi(3));
            d_m[ir] = rho[ir] * vol;
            m[ir] = m[ir - 1] + d_m[ir];
        }

        let mut iter = 0;
        let itermax = 100;
        let jmax = 20;
        let mix = 1.0;
        let mut delta = 0.0;

        let mut rhonew = vec![0.0; nr + 1];
        let mut g = vec![0.0; nr + 1];
        let mut p = vec![0.0; nr + 1];

        let mut rho_rock_comp = vec![0.0; nr + 1];
        let mut rho_hydr_comp = vec![0.0; nr + 1];
        let mut rho_h2os_comp = vec![0.0; nr + 1];
        let mut rho_adhs_comp = vec![0.0; nr + 1];
        let mut rho_h2ol_comp = vec![0.0; nr + 1];

        while iter < itermax {
            iter += 1;

            for ir in 1..=nr {
                g[ir] = G * m[ir] / r[ir].powi(2);
            }

            p[nr] = 0.0;
            for ir in (0..nr).rev() {
                p[ir] = p[ir + 1] + 0.5 * (g[ir + 1] + g[ir]) * (r[ir + 1] - r[ir]) * rho[ir + 1];
            }

            for ir in 1..=nr {
                let p_avg = 0.5 * (p[ir] + p[ir - 1]);
                let get_dy_dx = |entry: &PlanMatEntry, x_val: f64, z_val: f64| {
                    let (x2, _, x4, x5) = (2..=5).map(|n| x_val.powi(n)).next_tuple().unwrap();
                    let dy = 1.5
                        * entry.ks_0
                        * x5
                        * (z_val * x4 + (1.0 - 2.0 * z_val) * x2 + (z_val - 1.0))
                        - p_avg;
                    let dydx = 1.5
                        * entry.ks_0
                        * x4
                        * (9.0 * z_val * x4 + 7.0 * (1.0 - 2.0 * z_val) * x2 + 5.0 * (z_val - 1.0));
                    (dy, dydx)
                };
                if specify {
                    let entry = &planmat_db[icomp[ir]];
                    rhonew[ir] = match entry.eos {
                        1 => {
                            let z_val = 0.75 * (entry.ks_p - 4.0);
                            let mut x_val = 1_f64;
                            for _ in 0..jmax {
                                let (dy, dydx) = get_dy_dx(entry, x_val, z_val);
                                x_val -= dy / dydx;
                            }
                            entry.rho_0 * x_val.powi(3)
                        }
                        2 => entry.rho_0 + entry.c * p_avg.powf(entry.nn),
                        3 => {
                            let out = &thermal_outputs[ir - 1][time_step];
                            let temp_k = out.temp_kelvin;
                            let val = entry.v_0
                                * (1.0 + entry.a[0] * (entry.a[1] * (temp_k - entry.t_ref)).tanh())
                                * (entry.b[0] + entry.b[1] * (1.0 - (entry.b[2] * p[ir]).tanh()));
                            val.inv()
                        }
                        _ => {
                            println!("Compression: Error: specify EOS type in database");
                            return None;
                        }
                    }
                } else {
                    let rock_idx = planmat_index(&planmat_db, 206)?;
                    let rock_entry = &planmat_db[rock_idx];
                    let z_val = 0.75 * (rock_entry.ks_p - 4.0);
                    let mut x_val = 1_f64;
                    for _ in 0..jmax {
                        let (dy, dydx) = get_dy_dx(rock_entry, x_val, z_val);
                        x_val -= dy / dydx;
                    }
                    rho_rock_comp[ir] = rock_entry.rho_0 * x_val.powi(3);

                    let hydr_idx = planmat_index(&planmat_db, 305)?;
                    let hydr_entry = &planmat_db[hydr_idx];
                    let z_val = 0.75 * (hydr_entry.ks_p - 4.0);
                    let mut x_val = 1_f64;
                    for _ in 0..jmax {
                        let (dy, dydx) = get_dy_dx(hydr_entry, x_val, z_val);
                        x_val -= dy / dydx;
                    }
                    rho_hydr_comp[ir] = hydr_entry.rho_0 * x_val.powi(3);

                    let out = &thermal_outputs[ir - 1][time_step];
                    let temp_k = out.temp_kelvin;

                    let val = |mat: usize| {
                        let idx = planmat_index(&planmat_db, mat)?;
                        let entry = &planmat_db[idx];
                        Some(
                            (entry.v_0
                                * (1.0
                                    + entry.a[0] * ((temp_k - entry.t_ref) * entry.a[1]).tanh())
                                * (entry.b[0] + entry.b[1] * (1.0 - (entry.b[2] * p[ir]).tanh())))
                            .inv(),
                        )
                    };
                    rho_h2os_comp[ir] = val(403)?;
                    rho_adhs_comp[ir] = val(412)?;
                    rho_h2ol_comp[ir] = val(402)?;

                    let sum_mass = out.mass_total();
                    let val_mix = out.mass_rock
                        * (out.deg_of_hydr / rho_hydr_comp[ir]
                            + (1.0 - out.deg_of_hydr) / rho_rock_comp[ir])
                        + out.mass_ice / rho_h2os_comp[ir]
                        + out.mass_ammonia_solid / rho_adhs_comp[ir]
                        + (out.mass_water + out.mass_ammonia_liquid) / rho_h2ol_comp[ir];
                    rhonew[ir] = (1.0 / val_mix) * sum_mass;
                }
            }

            delta = (rhonew[1] / rho[1] - 1.).abs();

            for ir in 1..=nr {
                rho[ir] = rhonew[ir] * mix + rho[ir] * (1.0 - mix);
            }

            r[0] = 0.0;
            for ir in 1..=nr {
                r[ir] = (0.75 * d_m[ir] / rho[ir] / PI + r[ir - 1].powi(3)).powf(1.0 / 3.0);
            }

            if delta < 1.0e-9 {
                break;
            }
        }

        let rp_final = r[nr];
        let rhoavg = 0.75 * mp / PI / rp_final.powi(3);

        let output_dir = PathBuf::from("Outputs");
        fs::create_dir_all(&output_dir).ok()?;
        let output_file = output_dir.join("Compression.txt");
        let mut fout = fs::File::create(&output_file).ok()?;

        writeln!(
            fout,
            "Density (kg m-3) and pressure (MPa) profiles accounting for material compression"
        )
        .ok()?;
        writeln!(fout, "Radius (km) \t Density (kg/m3) \t Pressure (MPa)\n").ok()?;
        for ir in 1..=nr {
            writeln!(fout, "{} \t {} \t {}", r[ir] / KM, rho[ir], p[ir] / MPA).ok()?;
        }

        writeln!(fout, "\nAfter {} iterations ", iter).ok()?;
        writeln!(fout, "Convergence criterion = {}", delta).ok()?;
        writeln!(fout, "Planet mass = {} kg = {} MEarth", mp, mp / M_EARTH).ok()?;
        writeln!(
            fout,
            "Inner core mass = {} kg = {} MEarth",
            m_incore,
            m_incore / M_EARTH
        )
        .ok()?;
        writeln!(
            fout,
            "Outer core mass = {} kg = {} MEarth",
            m_outcore,
            m_outcore / M_EARTH
        )
        .ok()?;
        writeln!(
            fout,
            "Mantle mass = {} kg = {} MEarth",
            m_mantle,
            m_mantle / M_EARTH
        )
        .ok()?;
        writeln!(
            fout,
            "Inner core radius = {} km = {} REarth",
            r[nic] / KM,
            r[nic] / R_EARTH
        )
        .ok()?;
        writeln!(
            fout,
            "Outer core radius = {} km = {} REarth",
            r[noc] / KM,
            r[noc] / R_EARTH
        )
        .ok()?;
        writeln!(
            fout,
            "Planet radius = {} km = {} REarth",
            rp_final / KM,
            rp_final / R_EARTH
        )
        .ok()?;
        writeln!(fout, "Surface gravity = {} m/s2", g[nr]).ok()?;
        writeln!(fout, "Average density = {} kg/m3", rhoavg).ok()?;
        writeln!(fout, "Central density = {} kg/m3", rho[1]).ok()?;
        writeln!(fout, "Density near surface = {} kg/m3", rho[nr - 1]).ok()?;
        writeln!(fout, "Central pressure = {} MPa", p[0] / MPA).ok()?;
        writeln!(
            fout,
            "Pressure at core-mantle boundary = {} MPa",
            p[noc] / MPA
        )
        .ok()?;

        let mut d_m_thermal = vec![0.0; nr];
        let mut r_thermal = vec![0.0; nr + 1];
        let mut x_hydr = vec![0.0; nr];

        for ir in 0..nr {
            let out = &thermal_outputs[ir][time_step];
            d_m_thermal[ir] = out.mass_total();
            r_thermal[ir + 1] = out.radius_km;
            x_hydr[ir] = out.deg_of_hydr;
        }

        let mut m_cum = vec![0.0; nr];
        m_cum[0] = d_m_thermal[0];
        for ir in 1..nr {
            m_cum[ir] = m_cum[ir - 1] + d_m_thermal[ir];
        }
        let mut g_thermal = vec![0.0; nr];
        for ir in 0..nr {
            g_thermal[ir] = G * m_cum[ir] * GRAM / (r_thermal[ir + 1] / KM2CM * KM).powi(2);
        }

        let mut p_thermal = vec![0.0; nr];
        p_thermal[nr - 1] = 0.0;
        for ir in (0..=nr - 2).rev() {
            let Fracs(frock, fh2os, fadhs, fh2ol, fnh3l) = &thermal_outputs[ir][time_step].fracs();

            let avg_dens = frock * (x_hydr[ir] * rho_hydr + (1.0 - x_hydr[ir]) * rho_dry)
                + fh2os * 935.0
                + fh2ol * 1000.0
                + fadhs * 985.0
                + fnh3l * 740.0;

            p_thermal[ir] = p_thermal[ir + 1]
                + 0.5 * (g_thermal[ir + 1] + g_thermal[ir]) * (r_thermal[ir + 1] - r_thermal[ir])
                    / KM2CM
                    * KM
                    * avg_dens;
        }

        writeln!(
            fout,
            "\nDensity (kg m-3) and pressure (MPa) profiles from the thermal code\n"
        )
        .ok()?;
        writeln!(fout, "Radius (km) \t Density (kg/m3) \t Pressure (MPa)").ok()?;
        for ir in 1..nr {
            let out_curr = &thermal_outputs[ir][time_step];
            let out_prev = &thermal_outputs[ir - 1][time_step];
            let r_curr_km = out_curr.radius_km / KM2CM;
            let r_prev_km = out_prev.radius_km / KM2CM;

            let rho_val = d_m_thermal[ir] * GRAM
                / (4.0 / 3.0 * PI * (r_curr_km.powi(3) - r_prev_km.powi(3)) * KM.powi(3));

            writeln!(
                fout,
                "{} \t {} \t {}",
                r_prev_km,
                rho_val,
                p_thermal[ir] / MPA
            )
            .ok()?;
        }

        Some(())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PlanMatEntry {
    pub db_index: usize,
    pub eos: i32,
    pub rho_0: f64,
    pub c: f64,
    pub nn: f64,
    pub ks_0: f64,
    pub ks_p: f64,
    pub v_0: f64,
    pub t_ref: f64,
    pub a: [f64; 2],
    pub b: [f64; 3],
}

impl PlanMatEntry {
    pub fn from(entry: &[&str]) -> Option<Self> {
        if entry.len() < 14 {
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

pub fn planmat(data_folder: &str, _n_comp: usize) -> Option<Vec<PlanMatEntry>> {
    let planmat_db = PathBuf::from(data_folder).join("Compression_planmat.txt");
    let planmat_db = fs::read_to_string(planmat_db).ok()?;
    Some(
        planmat_db
            .lines()
            .skip(4)
            .enumerate()
            .filter_map(|(n, s)| {
                if !n.is_multiple_of(2) {
                    // skip material lines
                    None
                } else {
                    PlanMatEntry::from(&s.split_whitespace().collect_vec())
                }
            })
            .sorted_by_key(|p| p.db_index)
            .collect(),
    )
}

pub fn planmat_index(planmat_db: &[PlanMatEntry], mat: usize) -> Option<usize> {
    // planmat_db.iter().position(|e| e.db_index == mat)
    planmat_db.binary_search_by_key(&mat, |e| e.db_index).ok()
}
