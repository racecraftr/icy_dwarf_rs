use crate::{
    consts::*,
    create_output,
    input::{IcyDwarfInput, IcyWorld, Saturn, TidalQ, WorldSpec},
};

// pub fn planet_system(parsed: &ParsedInput) {
//     let base_mtx = vec![vec![0; parsed.grid.n_zones]; parsed.worlds.len()];
//     let base_vec = vec![0; parsed.worlds.len()];
// }

const TODAY: f64 = 4568.2 * MYR2SEC;
impl IcyDwarfInput {
    pub fn planet_sytem(&self, output_path: Option<String>) {
        // the input contains all the planets themselves. So not that hard.
        // if self.core_crack.incl_thermal_mismatch {}

        let real_time = 0.0;
        let q_prim = self.saturn.tidal_q.q_prim(&self.worlds, real_time);

        if self.housekeeping.recover {
            return;
        }

        //-------------------------------------------------------------------
        //                Normal init (from input file only)
        //-------------------------------------------------------------------
        for (idx, world) in self.worlds.iter().enumerate() {
            let name = &world.name;
            for s in ["Thermal", "Heats", "Crack_depth_WR", "Crack_stresses"] {
                match create_output(output_path.clone(), format!("{}_{}_{}.txt", idx, name, s)) {
                    Ok(()) => {}
                    Err(s) => eprintln!("{}", s),
                }
            }
            if self.saturn.mass > 0.0 {
                match create_output(output_path.clone(), format!("{}_{}_Orbit.txt", idx, name)) {
                    Ok(()) => {}
                    Err(s) => eprintln!("{}", s),
                }
            }
        }

        if self.saturn.mass > 0.0 {
            for s in ["Primary, Resonances, ResAcctFor, PCapture, icydwarf_outputs_1"] {
                match create_output(output_path.clone(), format!("{}.txt", s)) {
                    Ok(()) => {}
                    Err(s) => eprintln!("{}", s),
                }
            }
        }
    }
}

const RHO_H2OS_TH: f64 = RHO_H2OS * GRAM;
const RHO_ADHS_TH: f64 = RHO_ADHS * GRAM;

impl IcyWorld {
    pub fn mass(&self) -> f64 {
        self.planetary_dens * 4.0 / 3.0 * PI_GREEK * self.planetary_rad.powi(3)
    }

    pub fn rho_ice(&self) -> f64 {
        1.0 / (self.ammonia / XC) / (RHO_ADHS * GRAM)
    }

    pub fn phi(&self) -> f64 {
        todo!()
    }

    pub fn frac_rock_mass(&self, world_spec: &WorldSpec) -> f64 {
        let rho_ice = self.rho_ice();
        (1. - rho_ice / self.planetary_dens)
            / (1. - rho_ice / self.hydr_init * world_spec.rho_rock_th()
                + (1. - self.hydr_init) * world_spec.rho_hydr_th())
    }

    pub fn frac_rock_volume(&self, world_spec: &WorldSpec) -> f64 {
        let frac_rock_mass = self.frac_rock_mass(world_spec);
        frac_rock_mass * self.planetary_dens
            / (self.hydr_init * world_spec.rho_hydr_th()
                + (1. - self.hydr_init * world_spec.rho_rock_th()))
    }
}

impl WorldSpec {
    pub fn rho_hydr_th(&self) -> f64 {
        self.rho_rock_hydr * GRAM
    }

    pub fn rho_rock_th(&self) -> f64 {
        self.rho_rock_dry * GRAM
    }
}

impl TidalQ {
    pub fn q_prim(&self, planets: &Vec<IcyWorld>, real_time: f64) -> f64 {
        let t_form_min = planets.iter().fold(f64::MAX, |acc, e| acc.min(e.t_form));
        let scaling = (real_time - t_form_min) / (TODAY - t_form_min);
        match self.mode {
            crate::input::QMode::Lin => self.init + (self.today - self.init) * scaling,
            crate::input::QMode::ExpDecay => {
                self.init * ((self.today / self.init).ln() * scaling).exp()
            }
            crate::input::QMode::ExpChange => {
                self.init + 1. - self.init * ((self.today - self.init + 1.).ln() * scaling).exp()
            }
        }
    }
}
