use std::{fs, path::Path};

use crate::{consts::*, input::IcyDwarfInput, planet_system::ZoneState};

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

pub struct CrackOutput {
    pub pressure: f64,
    pub brittle_strength: f64,
    pub crit_intensity: f64,
    pub stress_intensity: f64,
    /// Pore water pressure in MPa (named porosity to match output mappings)
    pub porosity: f64,
    /// Hydration stress/pressure in MPa (named deg_of_hydration to match output mappings)
    pub deg_of_hydration: f64,
    pub crack_size_hydr_old: f64,
    pub crack_size_diss_old: f64,
    pub crack_size: f64,

    // Updated states returned as part of the output
    pub crack: f64,
    pub p_pore: f64,
    pub p_hydr: f64,
    pub act: [f64; 3],
}

impl IcyDwarfInput {
    /// Calculation of the depth and profile of cracking over time.
    pub fn crack(
        &self,
        zone: &mut ZoneState,
        dtime: f64,
        data: &Data,
        circ: bool,
    ) -> CrackOutput {
        let warnings = self.housekeeping.warnings;
        let thermal_mismatch = self.core_crack.incl_therm_mismatch;
        let pore_water_expansion = self.core_crack.incl_pore;
        let hydration_dehydration = self.core_crack.incl_hydr;
        let dissolution_precipitation = self.core_crack.incl_dissol && circ;

        let rho_hydr = self.world_spec.rho_rock_hydr;
        let rho_rock = self.world_spec.rho_rock_dry;
        let (brittle_strength, _) = strain(zone.pressure, zone.x_hydr, zone.temp, zone.porosity);

        let mut dtdt = 0.0;
        let e_young = zone.x_hydr * E_YOUNG_SERP + (1.0 - zone.x_hydr) * E_YOUNG_OLIV;
        let nu_poisson = zone.x_hydr * NU_POISSON_SERP + (1.0 - zone.x_hydr) * NU_POISSON_OLIV;
        let k_ic = zone.x_hydr * K_IC_SERP + (1.0 - zone.x_hydr) * K_IC_OLIV;

        // Initialize mutable output
        let mut output = CrackOutput {
            pressure: zone.pressure / MPA,
            brittle_strength: brittle_strength / MPA,
            crit_intensity: k_ic,
            stress_intensity: 0.0,
            porosity: zone.p_pore / MPA,
            deg_of_hydration: zone.p_hydr / MPA,
            crack_size_hydr_old: 0.0,
            crack_size_diss_old: 0.0,
            crack_size: zone.crack_size,

            crack: zone.crack,
            p_pore: zone.p_pore,
            p_hydr: zone.p_hydr,
            act: zone.act,
        };

        //-------------------------------------------------------------------
        // Cracks open from thermal expansion / contraction mismatch
        // (Fredrich and Wong 1986, Vance et al. 2007)
        //-------------------------------------------------------------------
        if thermal_mismatch {
            dtdt = (zone.temp - zone.temp_old) / dtime;

            // Calculate T' in each layer over time, eq (2) of Vance et al. (2007)
            // T' is the temperature at zero stress from thermal mismatch
            if dtdt == 0.0 {
                dtdt = 1.0e-6 / dtime; // To ensure continuity of T', otherwise T'=0
            }

            let t_prime = QGBS
                / R_G
                / (12.0 * OMEGA * D0_DELTAB * e_young
                    / (3.0_f64.sqrt() * N_FIT * K_B * L_SIZE.powi(3) * dtdt.abs()))
                .ln();

            // Calculate the stress intensity K_I in each layer over time,
            // eq (4) of Vance et al. (2007)
            if t_prime != 0.0 {
                // Look up the right value of a(T,P) to use in eq(4)
                let delta_t_int = look_up(
                    (t_prime - zone.temp).abs(),
                    0.0,
                    DELTA_T_STEP,
                    SIZEA_TP,
                    warnings,
                );
                let p_int = look_up(zone.pressure, 0.0, P_STEP, SIZEA_TP, warnings);
                let integral_line = (data.atp[delta_t_int][p_int] / A_MIN) as usize; // Index in the integral table
                let atp_val = data.atp[delta_t_int][p_int];

                // Calculate K_I
                output.stress_intensity = (2.0 / (PI_GREEK * atp_val)).sqrt()
                    * data.integral[integral_line][1]
                    * e_young
                    * DELTA_ALPHA
                    / (2.0 * PI_GREEK * (1.0 - nu_poisson * nu_poisson))
                    * (t_prime - zone.temp).abs()
                    - zone.pressure * (PI_GREEK * atp_val).sqrt();
            }
        }

        //-------------------------------------------------------------------
        //               Cracks from hydration - dehydration
        //-------------------------------------------------------------------
        if hydration_dehydration && zone.x_hydr != zone.x_hydr_old {
            // Shrinking/widening of open cracks
            if output.crack > 0.0 {
                output.p_hydr = 0.0;
                // Initialize crack size
                if output.crack_size == 0.0 {
                    output.crack_size = SMALLEST_CRACK_SIZE;
                }
                output.crack_size_hydr_old = output.crack_size;
                let x_bar = (2.0 * 4.5e-5 * (-45.0e3 / (R_G * zone.temp)).exp() * dtime).sqrt();
                let num = zone.x_hydr_old * rho_hydr + (1.0 - zone.x_hydr_old) * rho_rock;
                let den = zone.x_hydr * rho_hydr + (1.0 - zone.x_hydr) * rho_rock;
                let d_crack_size = -2.0 * ((num / den).powf(0.333) - 1.0) * x_bar;
                if output.crack_size + d_crack_size < 0.0 {
                    output.p_hydr = e_young * (-d_crack_size - output.crack_size) / x_bar; // Residual rock swell builds up stresses
                    output.crack_size = 0.0; // Crack closes completely
                } else {
                    output.crack_size += d_crack_size;
                }
            } else {
                // Cracks may open if stresses develop as rock shrinks/swells
                let num = zone.x_hydr_old * rho_hydr + (1.0 - zone.x_hydr_old) * rho_rock;
                let den = zone.x_hydr * rho_hydr + (1.0 - zone.x_hydr) * rho_rock;
                output.p_hydr += 2.0 * e_young * ((num / den).powf(0.333) - 1.0);
            }
        }

        //-------------------------------------------------------------------
        //             Expansion of pore water as it is heated
        //            (Norton 1984, Le Ravalec and Guéguen 1994)
        //-------------------------------------------------------------------
        if pore_water_expansion {
            if zone.x_hydr >= 0.09 && zone.temp > zone.temp_old {
                // Look up the right value of alpha and beta, given P and T
                let tempk_int = look_up(zone.temp, TEMPK_MIN, DELTA_TEMPK, SIZEA_TP, warnings);
                let p_int = look_up(
                    zone.pressure / BAR,
                    P_BAR_MIN,
                    DELTA_P_BAR,
                    SIZEA_TP,
                    warnings,
                );
                // Calculate fluid overpressure from heating, including geometric effects (Le Ravalec & Guéguen 1994)
                output.p_pore += (1.0 + 2.0 * ASPECT_RATIO)
                    * data.alpha[tempk_int][p_int]
                    * (zone.temp - zone.temp_old)
                    / (data.beta[tempk_int][p_int] / BAR
                        + ASPECT_RATIO * 3.0 * (1.0 - 2.0 * nu_poisson) / e_young);
            }
        }

        //-------------------------------------------------------------------
        //          Dissolution / precipitation (Bolton et al. 1997)
        //-------------------------------------------------------------------
        if dissolution_precipitation {
            if output.crack > 0.0 {
                // Initialize crack size
                if output.crack_size == 0.0 {
                    output.crack_size = SMALLEST_CRACK_SIZE;
                }
                output.crack_size_diss_old = output.crack_size; // For output only
                let mut d_crack_size = 0.0;
                let surface_volume_ratio = 2.0 / output.crack_size; // Rimstidt and Barnes (1980) Fig. 6 for a cylinder/fracture

                // Use reaction constants at given T and P
                let tempk_int = look_up(
                    zone.temp,
                    TEMPK_MIN_SPECIES,
                    DELTA_TEMPK_SPECIES,
                    SIZEA_TP,
                    warnings,
                );
                let p_int = look_up(
                    zone.pressure / BAR,
                    P_BAR_MIN,
                    DELTA_P_BAR,
                    SIZEA_TP,
                    warnings,
                );

                let mut k_eq = [0.0; 3];
                k_eq[0] = 10.0_f64.powf(data.silica[tempk_int][p_int]);
                k_eq[1] = 10.0_f64.powf(data.chrysotile[tempk_int][p_int]);
                k_eq[2] = 10.0_f64.powf(data.magnesite[tempk_int][p_int]);

                let chem_time = 1.0e6;
                let itermax = 100;
                let dtime_chem = dtime / chem_time;

                let mu_xu = [MU_XU_SILICA, MU_XU_CHRYSOTILE, MU_XU_MAGNESITE];
                let nu_prod = [NU_PROD_SILICA, NU_PROD_CHRYSOTILE, NU_PROD_MAGNESITE];
                let ea_diss = [EA_SILICA, EA_CHRYSOTILE, EA_MAGNESITE];
                let molar_volume = [
                    MOLAR_VOLUME_SILICA,
                    MOLAR_VOLUME_CHRYSOTILE,
                    MOLAR_VOLUME_MAGNESITE,
                ];
                let crack_species = [
                    self.core_crack.dissol.of_silica,
                    self.core_crack.dissol.of_serp,
                    self.core_crack.dissol.of_carb,
                ];

                for i in 0..3 {
                    if crack_species[i] {
                        let mut iter = 0;
                        while iter < itermax {
                            iter += 1;

                            // (Act_prod in mol L-1 to scale with K, silica equation (i=0) assumes unit A/V).
                            // The Arrhenius term is equivalent to a dissociation rate constant kdiss in mol m-2 s-1.
                            let inner_term = (output.act[i] / RHO_H2OL).powf(nu_prod[i]) / k_eq[i];
                            let r_diss = surface_volume_ratio
                                * (-ea_diss[i] / (R_G * zone.temp)).exp()
                                * 1.0
                                * (1.0 - inner_term.powf(mu_xu[i]));

                            // Update crack size and update act[i]
                            if -nu_prod[i] * r_diss * dtime_chem > output.act[i] {
                                // Everything precipitates
                                d_crack_size -= output.act[i] / nu_prod[i] * molar_volume[i]
                                    / surface_volume_ratio;
                                output.act[i] = 0.0; // Can't have negative concentrations!
                                break;
                            } else {
                                if nu_prod[i] * r_diss * dtime_chem < 0.1 * output.act[i] {
                                    d_crack_size += r_diss * dtime / (iter as f64)
                                        * molar_volume[i]
                                        / surface_volume_ratio;
                                    output.act[i] += nu_prod[i] * r_diss * dtime / (iter as f64);
                                    break;
                                }
                                d_crack_size +=
                                    r_diss * dtime_chem * molar_volume[i] / surface_volume_ratio;
                                output.act[i] += nu_prod[i] * r_diss * dtime_chem; // We neglect the change in crack volume to calculate Act[i].
                            }
                        }
                    }
                }

                if output.crack_size + d_crack_size > 0.0 {
                    output.crack_size += d_crack_size;
                } else {
                    output.crack_size = 0.0; // Pore clogged
                    for a in output.act.iter_mut() {
                        *a = 0.0; // Reset old activity quotients
                    }
                }
            } else {
                // If the crack is closed, clear the old activity quotients
                for a in output.act.iter_mut() {
                    *a = 0.0;
                }
            }
        }

        // Record the updated pressure outputs before they potentially get reset by the following checks
        output.porosity = output.p_pore / MPA;
        output.deg_of_hydration = output.p_hydr / MPA;

        // Cases where cracks appear
        if thermal_mismatch {
            if output.stress_intensity >= k_ic && dtdt < 0.0 {
                output.crack = 1.0; // Cooling cracks
            }
            if output.stress_intensity >= k_ic && dtdt >= 0.0 {
                output.crack = 2.0; // Heating cracks
            }
        }
        if hydration_dehydration {
            if output.p_hydr.abs() > zone.pressure + brittle_strength {
                if output.p_hydr > 0.0 {
                    output.crack = 3.0; // Compressive hydration cracks
                } else {
                    output.crack = 4.0; // Dehydration cracks
                }
                output.p_hydr = 0.0;
            }
        }
        if pore_water_expansion {
            if output.p_pore > brittle_strength {
                output.crack = 5.0;
                output.p_pore = 0.0;
            }
            if output.crack.floor() == 5.0 {
                output.p_pore = 0.0;
            }
        }
        if dissolution_precipitation {
            if output.crack > 0.0
                && output.crack_size > output.crack_size_diss_old
                && (output.crack == output.crack.floor()
                    || output.crack == output.crack.floor() + 0.2)
            {
                output.crack = output.crack.floor() + 0.1; // Dissolution widened crack
            }
            if output.crack > 0.0
                && output.crack_size < output.crack_size_diss_old
                && (output.crack == output.crack.floor()
                    || output.crack == output.crack.floor() + 0.1)
            {
                output.crack = output.crack.floor() + 0.2; // Precipitation shrunk crack
            }
        }

        // Cases where cracks disappear
        if zone.mass_rock <= zone.mass_rock_init {
            output.crack = 0.0; // Trivial: not enough rock
        }
        if hydration_dehydration {
            if output.p_hydr > 0.0 && output.p_hydr <= zone.pressure + brittle_strength {
                output.crack = -2.0; // Crack closed because of hydration
            }
        }
        if dissolution_precipitation {
            if output.crack > 0.0 && output.crack_size <= 0.0 {
                output.crack = -1.0; // Crack closed after precipitation
            }
        }

        // Write the updated state back to zone
        zone.crack = output.crack;
        zone.crack_size = output.crack_size;
        zone.p_pore = output.p_pore;
        zone.p_hydr = output.p_hydr;
        zone.act = output.act;

        output
    }
}

pub fn read_data(dat_folder: &String) -> Option<Data> {
    Some(Data {
        atp: read_data_file(&format!("{}/Crack_aTP.txt", dat_folder))?,
        alpha: read_data_file(&format!("{}/Crack_alpha.txt", dat_folder))?,
        beta: read_data_file(&format!("{}/Crack_beta.txt", dat_folder))?,
        chrysotile: read_data_file(&format!("{}/Crack_chrysotile.txt", dat_folder))?,
        integral: read_data_file(&format!("{}/Crack_integral.txt", dat_folder))?,
        magnesite: read_data_file(&format!("{}/Crack_magnesite.txt", dat_folder))?,
        silica: read_data_file(&format!("{}/Crack_silica.txt", dat_folder))?,
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

/// Helper lookup function to return correct index in a table.
fn look_up(x: f64, mut x_var: f64, x_step: f64, size: usize, warnings: bool) -> usize {
    if x <= x_step {
        0
    } else if x > x_var + x_step * ((size - 1) as f64) {
        if warnings {
            println!(
                "IcyDwarf look_up: x={} above range, assuming x={}",
                x,
                x_step * ((size - 1) as f64)
            );
        }
        size - 1
    } else {
        let mut x_int = 0;
        for j in 0..size {
            let lower = x_var - 0.5 * x_step;
            let upper = x_var + 0.5 * x_step;
            if lower > 0.0 && x > lower && x < upper {
                x_int = j;
            }
            x_var += x_step;
        }
        x_int
    }
}

/// Calculates the creep rate in s-1 of ice, rock, or a mixture using
/// flow laws. Stresses are hydrostatic pressure/(1-porosity).
pub fn creep(t: f64, p: f64, x_ice: f64, porosity: f64, x_hydr: f64) -> f64 {
    let p_por = p / MPA / (1.0 - porosity);
    let eps_disl = 4.0e5 * p_por.powi(4) * (-60.0e3 / (R_G * t)).exp();
    let eps_basal = if t < 255.0 {
        3.9e-3 * (-49.0e3 / (R_G * t)).exp()
    } else {
        3.0e26 * (-192.0e3 / (R_G * t)).exp()
    } * D_FLOW_LAW.powf(-1.4)
        * p_por.powf(1.8);
    let eps_gbs = 5.5e7 * p_por.powf(2.4) * (-60.0e3 / (R_G * t)).exp();
    let eps_diff = 3.02e-14 * p_por * D_FLOW_LAW.powi(-2) * (-59.4e3 / (R_G * t)).exp();

    let creep_rate_ice = eps_diff + 1.0 / (1.0 / eps_basal + 1.0 / eps_gbs) + eps_disl;

    if x_ice > 0.3 {
        // The rock fragments are barely in contact and deformation is controlled entirely by the ice
        creep_rate_ice
    } else {
        // Deformation is controlled by both rock and ice properties
        let t_calc = t.max(140.0);
        let creep_rate_hydr =
            416869.38347 * p_por * D_FLOW_LAW.powi(-3) * (-240.0e3 / (R_G * t_calc)).exp();
        let creep_rate_dry = 177827.941004
            * p_por
            * D_FLOW_LAW.powf(-2.98)
            * ((-261.0e3 + p * 6.0e-6) / (R_G * t_calc)).exp();

        // Scaling from Roberts (2015)
        let term1 =
            (0.3 - x_ice) * (x_hydr * creep_rate_hydr + (1.0 - x_hydr) * creep_rate_dry).ln();
        let term2 = x_ice * creep_rate_ice.ln();
        ((term1 + term2) / 0.3).exp()
    }
}

/// Calculates the brittle strength in Pa and corresponding ductile
/// strain rate in s-1 of silicate rock.
pub fn strain(pressure: f64, x_hydr: f64, t: f64, porosity: f64) -> (f64, f64) {
    let hydr_strength = MU_F_SERP * pressure;
    let dry_strength = if pressure <= 200.0e6 {
        MU_F_BYERLEE_LOP * pressure
    } else {
        MU_F_BYERLEE_HIP * pressure + C_F_BYERLEE_HIP
    };
    let mut brittle_strength = x_hydr * hydr_strength + (1.0 - x_hydr) * dry_strength;
    brittle_strength /= 1.0 - porosity;

    let t_calc = t.max(140.0);
    let strain_rate = 416869.38347 * brittle_strength / MPA
        * D_FLOW_LAW.powi(-3)
        * (-240.0e3 / (R_G * t_calc)).exp();

    (brittle_strength, strain_rate)
}
