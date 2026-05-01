use std::{
    fs::{self, File},
    io::{self, BufRead},
    path::Path,
};

use serde::Deserialize;
use serde_repr::Deserialize_repr;

#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum QMode {
    #[default]
    Lin,
    ExpDecay,
    ExpChange,
}

#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum TidalModel {
    #[default]
    NA,
    Elastic,
    Maxwell,
    Burgers, // burgers ??? lmaooo
    Andrade,
}

#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum EccModel {
    #[default]
    E2,

    E10Cpl,
    E10Ctl,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct Housekeeping {
    pub warnings: bool,
    pub recover: bool,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct Grid {
    pub n_zones: usize,
    pub time_step: f64,
    pub speedup: f64,
    pub time_total: f64,
    pub output_every: f64,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct TidalQ {
    pub init: f64,
    pub today: f64,
    pub model: TidalModel,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct Ring {
    pub mass: f64,
    pub inner: f64,
    pub outer: f64,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct Saturn {
    pub mass: f64,
    pub rad: f64,
    pub moi_coef: f64,
    pub tidal_q: TidalQ,
    pub k2: f64,
    pub j2: f64,
    pub j4: f64,
    pub tidal_resonant: bool,
    pub spin_period: f64,
    pub ring: Ring,
}

#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum ChondriteType {
    #[default]
    CO,
    CI,
}

#[repr(u8)]
#[derive(Default, Debug, Clone, Deserialize_repr)]
pub enum TidalRhel {
    #[default]
    Maxwell = 2,
    Burgers,
    Andr,
    SunCoop,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct WorldSpec {
    pub rho_rock_dry: f64,
    pub rho_rock_hydr: f64,
    pub chondrite: bool,
    pub rhelogy: TidalRhel,
    pub ecc_model: EccModel,
    pub tidal_heating: bool,
    pub lookup_tbl: Vec<f64>,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct IcyWorld {
    pub name: String,
    pub r_p: f64,
    pub rho_p: f64,
    pub temp_surf: f64,
    pub temp_init: f64,
    pub t_form: f64,
    pub from_ring: bool,
    pub ammonia: f64,
    pub briny: bool,
    pub hydr_init: f64,
    pub hydrate: bool,
    pub por_init: f64,
    pub rock_frac: f64,
    pub rock_h20: f64,
    pub start_diff: bool,
    pub orb_a_init: f64,
    pub orb_e_init: f64,
    pub orb_i_init: f64,
    pub orb_o_init: f64,
    pub orb_can_change: bool,
    pub retrograde: bool,
    pub t_reslock: f64,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct CoreCrackDissol {
    pub of_silica: bool,
    pub of_serp: bool,
    pub of_carb: bool,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct CoreCrack {
    pub incl_mismatch: bool,
    pub incl_pore: bool,
    pub incl_hydr: bool,
    pub incl_dissol: bool,
    pub dissol: CoreCrackDissol,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct GeoInput {
    pub min: f64,
    pub max: f64,
    pub step: f64,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct SubroutinesGeo {
    pub temp: GeoInput,
    pub pressure: GeoInput,
    pub pe: GeoInput,
    pub wr_ratio: GeoInput,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct SubroutinesCryo {
    pub after: f64,
    pub min_temp_chnosz: f64,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct Subroutines {
    pub run_therm: bool,
    pub gen_crack_core: bool,
    pub gen_water_ab: bool,
    pub gen_crack_sp: bool,
    pub run_geo: bool,
    pub run_comp: bool,
    pub run_cryo: bool,
    pub geo: SubroutinesGeo,
    pub cryo: SubroutinesCryo,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct ParsedInput {
    pub housekeeping: Housekeeping,
    pub grid: Grid,
    pub saturn: Saturn,
    pub world_spec: WorldSpec,
    pub worlds: Vec<IcyWorld>,
    pub subroutines: Subroutines,
    pub core_crack: CoreCrack,
}

/// parses a .toml file to a [`ParsedInput`].
/// [TOML](https://toml.io/) is a configuration file format, which can be used
/// to also define the inputs for IcyDwarf. It is both human and machine-readable,
/// and follows the v25.x format very closely.
/// See inputs/input.toml for an example.
pub fn parse_toml(toml_path: &str) -> Option<ParsedInput> {
    if !toml_path.ends_with(".toml") {
        println!("WARNING: File {} does not end with .toml", toml_path);
    }
    let Ok(toml_str) = fs::read_to_string(toml_path) else {
        eprintln!("ERROR: File {} does not exist!", toml_path);
        return None;
    };

    toml::from_str(&toml_str).ok()
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
