use crate::{
    consts::KM2CM,
    input::{Fracs, IcyDwarfInput},
    planet_system::{RHO_ADHS_TH, RHO_H2OL_TH, RHO_NH3L_TH},
};

/// This struct stores a snapshot of radial zone state loaded from thermal output files.
#[allow(dead_code)]
#[derive(Clone, Debug, Default)]
pub struct ThermalOut {
    /// Zone radius in centimeters.
    pub radius_km: f64,
    /// Temperature in Kelvin.
    pub temp_kelvin: f64,
    /// Rock mass in grams.
    pub mass_rock: f64,
    /// Water ice mass in grams.
    pub mass_ice: f64,
    /// Solid ammonia dihydrate mass in grams.
    pub mass_ammonia_solid: f64,
    /// Liquid water mass in grams.
    pub mass_water: f64,
    /// Liquid ammonia solution mass in grams.
    pub mass_ammonia_liquid: f64,
    /// Nusselt convection number.
    pub nusselt_num: f64,
    /// Amorphous ice fraction.
    pub ice_frac_amorphous: f64,
    /// Thermal conductivity.
    pub thermal_cond: f64,
    /// Degree of rock hydration.
    pub deg_of_hydr: f64,
    /// Matrix porosity fraction.
    pub porosity: f64,
    /// Core cracking flag.
    pub crack: bool,
    /// Tidal heating power rate in erg per second.
    pub tidal_heating_rate: f64,
}

impl ThermalOut {
    /// Parse a line of output text into a [`ThermalOut`] struct.
    ///
    /// # Parameters
    /// - `ln`: The line string to parse.
    ///
    /// # Returns
    /// An optional [`ThermalOut`] struct containing the parsed fields.
    pub fn from_line(ln: &str) -> Option<Self> {
        let parts = ln.split_whitespace().collect::<Vec<_>>();
        let radius_km = parts[0].parse::<f64>().ok()? * KM2CM;
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

    /// Calculate total zone mass in grams.
    ///
    /// # Returns
    /// The total mass sum.
    pub fn mass_total(&self) -> f64 {
        self.mass_rock
            + self.mass_ice
            + self.mass_ammonia_solid
            + self.mass_ammonia_liquid
            + self.mass_water
    }

    /// Calculate total zone volume and phase volumes.
    ///
    /// # Parameters
    /// - `input`: The [`IcyDwarfInput`] configuration reference.
    ///
    /// # Returns
    /// A tuple containing total volume and a [`Fracs`] struct of phase volumes.
    #[allow(dead_code)]
    pub fn vol(&self, input: &IcyDwarfInput) -> (f64, Fracs) {
        let vol_rock = self.mass_rock
            / (self.deg_of_hydr * input.world_spec.rho_hydr_th()
                + (1.0 - self.deg_of_hydr) * input.world_spec.rho_rock_th());
        let vol_ice = self.mass_ice / RHO_H2OL_TH;
        let vol_adhs = self.mass_ammonia_solid / RHO_ADHS_TH;
        let vol_water = self.mass_water / RHO_H2OL_TH;
        let vol_nh3l = self.mass_ammonia_liquid / RHO_NH3L_TH;
        (
            vol_rock + vol_adhs + vol_water + vol_nh3l,
            Fracs(vol_rock, vol_ice, vol_adhs, vol_water, vol_nh3l),
        )
    }

    /// Calculate phase mass fractions for the output zone.
    ///
    /// # Returns
    /// A [`Fracs`] struct containing phase mass fractions.
    pub fn fracs(&self) -> Fracs {
        let mass_total = self.mass_total();
        Fracs(
            self.mass_rock / mass_total,
            self.mass_ice / mass_total,
            self.mass_ammonia_solid / mass_total,
            self.mass_water / mass_total,
            self.mass_ammonia_liquid / mass_total,
        )
    }
}
