//! This module defines the input data structures and TOML parsing logic for the simulation.

use std::fs::{self};

use crate::consts::MYR2SEC;
use serde::{
    Deserialize,
    de::{self},
};
use serde_repr::Deserialize_repr;
pub mod recover;

/// This enum specifies the dissipation factor mode.
#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum QMode {
    /// Linear dissipation mode.
    #[default]
    Lin,
    /// Exponential decay mode.
    ExpDecay,
    /// Exponential change mode.
    ExpChange,
}

impl TryFrom<u8> for QMode {
    type Error = u8;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        if value < 3 {
            Ok(unsafe { std::mem::transmute::<u8, QMode>(value) })
        } else {
            Err(value)
        }
    }
}

impl<'a> TryFrom<&'a str> for QMode {
    type Error = String;

    fn try_from(value: &'a str) -> Result<Self, Self::Error> {
        match value.to_ascii_lowercase().trim() {
            "linear" | "lin" => Ok(Self::Lin),
            "decay" | "expdecay" | "exponential decay" => Ok(Self::ExpDecay),
            "expchange" | "exponential change" => Ok(Self::ExpChange),
            _ => Err(value.to_owned()),
        }
    }
}

fn deserialize_qmode<'de, D>(deserializer: D) -> Result<QMode, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum QModeInput {
        Int(u8),
        S(String),
    }

    match QModeInput::deserialize(deserializer)? {
        QModeInput::Int(n) => QMode::try_from(n).map_err(|v| {
            de::Error::custom(format!(
                "{} is not a valid dissipation factor model integer code",
                v
            ))
        }),
        QModeInput::S(s) => QMode::try_from(s.as_ref()).map_err(|s| {
            de::Error::custom(format!("{} is not a valid dissipation factor model", s))
        }),
    }
}

/// This enum specifies the orbital eccentricity model.
#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum EccModel {
    /// Standard second-order eccentricity model.
    #[default]
    E2,
    /// Coupled tenth-order eccentricity model.
    E10Cpl,
    /// Constant time-lag tenth-order eccentricity model.
    E10Ctl,
}

impl TryFrom<u8> for EccModel {
    type Error = u8;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        if value < 3 {
            Ok(unsafe { std::mem::transmute::<u8, EccModel>(value) })
        } else {
            Err(value)
        }
    }
}

impl<'a> TryFrom<&'a str> for EccModel {
    type Error = String;

    fn try_from(value: &'a str) -> Result<Self, Self::Error> {
        match value.to_ascii_lowercase().trim() {
            "second-order" | "e2" | "ecc2" => Ok(Self::E2),
            "coupled tenth-order" | "e10cpl" | "ecc10_cpl" => Ok(Self::E10Cpl),
            "time-lag tenth-order" | "e10ctl" | "ecc10_ctl" => Ok(Self::E10Ctl),
            _ => Err(value.to_owned()),
        }
    }
}

fn deserialize_ecc_model<'de, D>(deserializer: D) -> Result<EccModel, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum EccModelInput {
        Int(u8),
        S(String),
    }

    match EccModelInput::deserialize(deserializer)? {
        EccModelInput::Int(n) => EccModel::try_from(n).map_err(|v| {
            de::Error::custom(format!(
                "{} is not a valid eccentricity model integer code",
                v
            ))
        }),
        EccModelInput::S(s) => EccModel::try_from(s.as_ref())
            .map_err(|s| de::Error::custom(format!("{} is not a valid eccentricity model", s))),
    }
}

/// This struct holds configuration flags for simulation logging and recovery.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct Housekeeping {
    /// Flag to display simulation warnings.
    pub warnings: bool,
    /// Flag to recover simulation state from previous output files.
    pub recover: bool,
}

/// This struct defines the spatial and temporal grid settings for the simulation.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct Grid {
    /// The number of radial zones in each world.
    pub n_zones: usize,
    /// The length of a single time step in years.
    pub time_step: f64,
    /// The speedup multiplier factor.
    pub speedup: f64,
    /// The total duration of the simulation in megayears.
    pub time_total: f64,
    /// The frequency of output generation in megayears.
    pub output_every: f64,
}

impl Grid {
    /// Calculate the total number of output time steps in the simulation.
    ///
    /// # Returns
    /// The number of output time steps.
    pub fn output_time_step(&self) -> usize {
        (self.time_total / self.output_every) as usize
    }
}

/// This struct stores parameters for the tidal dissipation factor.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct TidalQ {
    /// The initial value of the tidal dissipation factor.
    pub init: f64,
    /// The current value of the tidal dissipation factor.
    pub today: f64,
    /// The mode of dissipation calculation.
    #[serde(deserialize_with = "deserialize_qmode")]
    pub mode: QMode,
}

/// This struct defines the physical properties of a planetary ring system.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct Ring {
    /// The total mass of the ring in grams.
    pub mass: f64,
    /// The inner radius of the ring in centimeters.
    pub inner: f64,
    /// The outer radius of the ring in centimeters.
    pub outer: f64,
}

/// This struct defines the physical properties of the primary central body.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct PrimaryWorld {
    /// The total mass of the primary body in grams.
    pub mass: f64,
    /// The volumetric radius of the primary body in centimeters.
    pub rad: f64,
    /// The dimensionless moment of inertia coefficient.
    pub moi_coef: f64,
    /// The tidal dissipation configuration for the primary body.
    pub tidal_q: TidalQ,
    /// The second Love number of the primary body.
    pub k2: f64,
    /// The second zonal harmonic coefficient.
    pub j2: f64,
    /// The fourth zonal harmonic coefficient.
    pub j4: f64,
    /// Flag to enable tidal resonant interactions.
    pub tidal_resonant: bool,
    /// The rotation period of the primary body in hours.
    pub spin_period: f64,
    /// The ring parameters associated with the primary body.
    pub ring: Ring,
}

/// This enum specifies the type of chondrite material.
#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr, PartialEq, Eq)]
pub enum ChondriteType {
    #[default]
    /// Carbonaceous Ivuna chondrite type.
    CI,
    /// Carbonaceous Ornans chondrite type.
    CO,
}

impl ChondriteType {
    fn from_string(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().trim() {
            "ci" | "chondr_i" | "ivuna" => Some(Self::CI),
            "co" | "chondr_o" | "orans" => Some(Self::CO),
            _ => None,
        }
    }
}

/// Helper to deserialize chondrite from boolean or integer in TOML files.
fn deserialize_chondrite<'de, D>(deserializer: D) -> Result<ChondriteType, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum ChondriteInput {
        Bool(bool),
        Int(u8),
        S(String),
    }

    match ChondriteInput::deserialize(deserializer)? {
        ChondriteInput::Bool(b) => Ok(if b {
            ChondriteType::CO
        } else {
            ChondriteType::CI
        }),
        ChondriteInput::Int(n) => {
            if n == 0 {
                Ok(ChondriteType::CI)
            } else if n == 1 || n == 2 {
                Ok(ChondriteType::CO)
            } else {
                Err(de::Error::custom(format!(
                    "{} is an invalid chondrite type integer code.",
                    n
                )))
            }
        }
        ChondriteInput::S(s) => {
            if let Some(c_type) = ChondriteType::from_string(&s) {
                Ok(c_type)
            } else {
                Err(de::Error::custom(format!(
                    "{} is an invalid chondrite type.",
                    &s
                )))
            }
        }
    }
}

/// This enum specifies the rheology model used for tidal calculations.
#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr, PartialEq, Eq)]
pub enum TidalModel {
    /// Maxwell viscoelastic rheology model.
    #[default]
    Maxwell = 2,
    /// Burgers viscoelastic rheology model.
    Burgers,
    /// Andrade viscoelastic rheology model.
    Andr,
    /// Sundberg-Cooper viscoelastic rheology model.
    SunCoop,
}

impl TidalModel {
    fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().trim() {
            "maxwell" => Some(Self::Maxwell),
            "burgers" => Some(Self::Burgers),
            "andr" | "andrade" => Some(Self::Andr),
            "suncoop" | "sundberg-cooper" => Some(Self::SunCoop),
            _ => None,
        }
    }
}

fn deserialize_tidal_model<'de, D>(deserializer: D) -> Result<TidalModel, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum TidalModelInput {
        Int(u8),
        S(String),
    }

    match TidalModelInput::deserialize(deserializer)? {
        TidalModelInput::Int(n) => {
            if (2..=5).contains(&n) {
                Ok(unsafe { std::mem::transmute::<u8, TidalModel>(n) }) // this is safe, as both have the same size.
            } else {
                Err(de::Error::custom(format!(
                    "{} is an invalid tidal model integer code.",
                    &n
                )))
            }
        }
        TidalModelInput::S(s) => {
            if let Some(model) = TidalModel::from_str(&s) {
                Ok(model)
            } else {
                Err(serde::de::Error::custom(format!(
                    "{} is an invalid tidal model.",
                    &s
                )))
            }
        }
    }
}

/// This struct defines physical parameters shared across planetary bodies.
#[allow(dead_code)]
#[derive(Default, Debug, Clone, Deserialize)]
pub struct WorldSpec {
    /// The density of dry rock material in grams per cubic centimeter.
    pub rho_rock_dry: f64,
    /// The density of hydrated rock material in grams per cubic centimeter.
    pub rho_rock_hydr: f64,
    /// Flag to indicate chondrite composition.
    #[serde(deserialize_with = "deserialize_chondrite")]
    pub chondrite: ChondriteType,
    /// The tidal rheology model.
    #[serde(deserialize_with = "deserialize_tidal_model")]
    pub rhelogy: TidalModel,
    /// The eccentricity model.
    #[serde(deserialize_with = "deserialize_ecc_model")]
    pub ecc_model: EccModel,
    /// Flag to enable tidal heating calculation.
    pub tidal_heating: bool,
    /// The lookup table array for material properties.
    pub lookup_tbl: Vec<f64>,
}

/// This struct defines physical and orbital properties for an icy moon.
#[allow(dead_code)]
#[derive(Default, Debug, Clone, Deserialize)]
pub struct IcyWorld {
    /// The name of the world.
    pub name: String,
    /// The total radius of the world in kilometers.
    pub planetary_rad: f64,
    /// The mean bulk density of the world in grams per cubic centimeter.
    pub planetary_dens: f64,
    /// The surface boundary temperature in Kelvin.
    pub temp_surf: f64,
    /// The initial interior temperature in Kelvin.
    pub temp_init: f64,
    /// The formation time in millions of years after solar system formation.
    pub t_form: f64,
    /// Flag indicating if the world formed from a planetary ring.
    pub from_ring: bool,
    /// The initial ammonia mass fraction.
    pub ammonia: f64,
    /// Flag indicating if the ocean contains salt.
    pub briny: bool,
    /// The initial degree of hydration.
    pub hydr_init: f64,
    /// Flag indicating if gas hydrates are present.
    pub hydrate: bool,
    /// The initial rock matrix porosity.
    pub por_init: f64,
    /// The mass fraction of rock in the body.
    pub rock_frac: f64,
    /// The initial water content ratio in rock.
    pub rock_h20: f64,
    /// Flag to start the body in a differentiated state.
    pub start_diff: bool,
    /// The initial semi-major axis of the orbit in kilometers.
    pub orb_a_init: f64,
    /// The initial eccentricity of the orbit.
    pub orb_e_init: f64,
    /// The initial orbital inclination in degrees.
    pub orb_i_init: f64,
    /// The initial obliquity in degrees.
    pub orb_o_init: f64,
    /// Flag to allow orbital migration.
    pub orb_can_change: bool,
    /// Flag indicating a retrograde orbit.
    pub retrograde: bool,
    /// The resonance locking duration in millions of years.
    pub t_reslock: f64,
}

/// This struct configures dissolution options for core cracking.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct CoreCrackDissol {
    /// Flag to include silica dissolution.
    pub of_silica: bool,
    /// Flag to include serpentine dissolution.
    pub of_serp: bool,
    /// Flag to include carbonate dissolution.
    pub of_carb: bool,
}

/// This struct configures mechanisms for core cracking simulation.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct CoreCrack {
    /// Flag to include thermal expansion mismatch stresses.
    pub incl_therm_mismatch: bool,
    /// Flag to include pore pressure stresses.
    pub incl_pore: bool,
    /// Flag to include hydration volume expansion stresses.
    pub incl_hydr: bool,
    /// Flag to include mineral dissolution stresses.
    pub incl_dissol: bool,
    /// Detailed mineral dissolution configuration.
    pub dissol: CoreCrackDissol,
}

/// This struct defines bounds and steps for geochemical grid searches.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct GeoInput {
    /// The minimum bound value.
    pub min: f64,
    /// The maximum bound value.
    pub max: f64,
    /// The step size for grid search.
    pub step: f64,
}

/// This struct configures geochemical subroutine parameters.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct SubroutinesGeo {
    /// Temperature search parameters.
    pub temp: GeoInput,
    /// Pressure search parameters.
    pub pressure: GeoInput,
    /// Electron activity search parameters.
    pub pe: GeoInput,
    /// Water-to-rock ratio search parameters.
    pub wr_ratio: GeoInput,
}

/// This struct configures cryolava subroutine parameters.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct SubroutinesCryo {
    /// The start time for cryolava calculations in steps.
    pub after: i32,
    /// The minimum temperature for CHNOSZ calculations in Kelvin.
    pub min_temp_chnosz: f64,
}

/// This struct holds flags to enable or disable simulation subroutines.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct Subroutines {
    /// Flag to run thermal evolution calculations.
    pub run_therm: bool,
    /// Flag to generate core cracking output.
    pub gen_crack_core: bool,
    /// Flag to generate water abundance output.
    pub gen_water_ab: bool,
    /// Flag to generate core crack stress intensity tables.
    pub gen_crack_sp: bool,
    /// Flag to run geochemical modeling routines.
    pub run_geo: bool,
    /// Flag to run core compression calculations.
    pub run_comp: bool,
    /// Flag to run cryolava modeling routines.
    pub run_cryo: bool,
    /// Geochemical routine configuration settings.
    pub geo: SubroutinesGeo,
    /// Cryolava routine configuration settings.
    pub cryo: SubroutinesCryo,
}

/// This struct contains the complete input parameters loaded from a configuration file.
#[derive(Default, Debug, Clone, Deserialize)]
pub struct IcyDwarfInput {
    /// Housekeeping configuration settings.
    pub housekeeping: Housekeeping,
    /// Spatial and temporal grid settings.
    pub grid: Grid,
    /// Primary central world properties.
    pub primary_world: PrimaryWorld,
    /// Physical specification options shared across worlds.
    pub world_spec: WorldSpec,
    /// List of secondary icy worlds in the simulation.
    pub worlds: Vec<IcyWorld>,
    /// Flags and settings for simulation subroutines.
    pub subroutines: Subroutines,
    /// Core cracking mechanism configuration.
    pub core_crack: CoreCrack,
}

impl IcyDwarfInput {
    /// Return the initial hydration vector for all zones across worlds.
    ///
    /// # Returns
    /// A nested vector containing initial hydration fractions.
    #[allow(dead_code)]
    pub fn x_hydr(&self) -> Vec<Vec<f64>> {
        self.worlds
            .iter()
            .map(|w| vec![w.hydr_init; self.grid.n_zones])
            .collect()
    }

    /// Calculate the step interval at which cryolava calculations run.
    ///
    /// # Returns
    /// The step interval for cryolava updates.
    pub fn t_cryo(&self) -> i32 {
        if self.grid.output_every > 0.0 {
            ((self.subroutines.cryo.after as f64 * MYR2SEC) / self.grid.output_every).floor() as i32
        } else {
            0
        }
    }

    /// Return the total number of secondary worlds in the simulation.
    ///
    /// # Returns
    /// The count of secondary worlds.
    pub fn n_moons(&self) -> usize {
        self.worlds.len()
    }
}

/// This struct holds phase mass fractions for rock, water ice, adhs, water, and ammonia liquid.
#[derive(Default, Debug, Clone)]
pub struct Fracs(pub f64, pub f64, pub f64, pub f64, pub f64);

/// Parse a TOML configuration file into an [`IcyDwarfInput`] struct.
///
/// # Parameters
/// - `toml_path`: The file system path to the input TOML file.
///
/// # Returns
/// An optional [`IcyDwarfInput`] containing the parsed data.
pub fn parse_toml(toml_path: &str) -> Option<IcyDwarfInput> {
    if !toml_path.ends_with(".toml") {
        println!("WARNING: File {} does not end with .toml", toml_path);
    }
    let Ok(toml_str) = fs::read_to_string(toml_path) else {
        eprintln!("ERROR: File {} does not exist!", toml_path);
        return None;
    };

    match toml::from_str(&toml_str) {
        Ok(val) => Some(val),
        Err(err) => {
            eprintln!("ERROR parsing TOML {}: {}", toml_path, err);
            None
        }
    }
}

#[cfg(test)]
mod test {
    use crate::input::parse_toml;

    #[test]
    fn test_input() {
        let parsed = parse_toml("inputs/input.toml");
        println!("{:#?}", parsed);
        assert!(parsed.is_some());
    }
}
