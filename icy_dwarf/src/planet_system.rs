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

pub const RHO_H2OS_TH: f64 = RHO_H2OS * GRAM;
pub const RHO_H2OL_TH: f64 = RHO_H2OS_TH;
pub const RHO_ADHS_TH: f64 = RHO_ADHS * GRAM;
pub const RHO_NH3L_TH: f64 = XC / (1.0 / RHO_H2OL_TH) + (1.0 / RHO_ADHS_TH - 1.0 / RHO_H2OS_TH);

#[derive(Clone, Debug)]
pub struct ZoneState {
    pub radius: f64,
    pub dr: f64,
    pub temp: f64,
    pub temp_old: f64,
    pub mass_total: f64, // usually d_m in most cases
    pub mass_rock: f64,
    pub mass_rock_init: f64,
    pub mass_h2os: f64,
    pub mass_adhs: f64,
    pub mass_h2ol: f64,
    pub mass_nh3l: f64,
    pub energy_total: f64,
    pub porosity: f64,
    pub pressure: f64,
    pub x_hydr: f64,
    pub x_hydr_old: f64,
    pub kappa: f64,
    pub nusselt: f64,
}

#[derive(Clone, Debug)]
pub struct WorldState {
    pub name: String,
    pub a_orb: f64,
    pub e_orb: f64,
    pub i_orb: f64,
    pub obl: f64,
    pub n_orb: f64,
    pub lambda: f64,
    pub omega: f64,
    pub spin: f64,
    pub k2: f64,
    pub q_tide: f64,
    pub moi: f64,
    pub zones: Vec<ZoneState>,
    pub res_acct_for: Vec<f64>,
}

impl IcyDwarfInput {
    pub fn planet_system(&self, output_path: Option<String>) {
        let dtime = self.grid.time_step * 1.0e-6 * MYR2SEC;
        let n_time = (self.grid.time_total / self.grid.time_step) as usize;
        let n_steps = (self.grid.output_every / self.grid.time_step) as usize;

        let t_form_min = self
            .worlds
            .iter()
            .fold(f64::MAX, |acc, w| acc.min(w.t_form));

        //-------------------------------------------------------------------
        //                Normal init (from input file only)
        //-------------------------------------------------------------------
        let mut world_states: Vec<WorldState> = self
            .worlds
            .iter()
            .enumerate()
            .map(|(idx, world)| {
                let name = &world.name;
                for s in ["Thermal", "Heats", "Crack_depth_WR", "Crack_stresses"] {
                    let _ =
                        create_output(output_path.clone(), format!("{}_{}_{}.txt", idx, name, s));
                }
                if self.saturn.mass > 0.0 {
                    let _ =
                        create_output(output_path.clone(), format!("{}_{}_Orbit.txt", idx, name));
                }

                let rho_ice = 1.0
                    / ((world.ammonia / XC) / RHO_ADHS_TH
                        + (1.0 - world.ammonia / XC) / RHO_H2OS_TH);
                let rho_rock_hydr = world.hydr_init * self.world_spec.rho_hydr_th()
                    + (1.0 - world.hydr_init) * self.world_spec.rho_rock_th();
                let frock_pm =
                    (1.0 - rho_ice / world.planetary_dens) / (1.0 - rho_ice / rho_rock_hydr);

                let dr_grid = world.planetary_rad / (self.grid.n_zones as f64);
                let mut zones = Vec::with_capacity(self.grid.n_zones);
                let mut current_r = 0.0;

                for _ir in 0..self.grid.n_zones {
                    let next_r = current_r + dr_grid;
                    let d_vol = 4.0 / 3.0 * PI_GREEK * (next_r.powi(3) - current_r.powi(3));
                    let d_m = d_vol * world.planetary_dens;
                    let m_rock = d_m * frock_pm;
                    let m_h2os = d_m * (1.0 - frock_pm) * (1.0 - world.ammonia / XC);
                    let m_adhs = d_m * (1.0 - frock_pm) * (world.ammonia / XC);

                    let e_rock = m_rock * crate::heat_rock(world.temp_init);
                    let e_h2os = m_h2os * QH2O * world.temp_init * world.temp_init / 2.0;
                    let e_slush = m_adhs * QADH * world.temp_init * world.temp_init / 2.0;
                    let d_e = e_rock + e_h2os + e_slush;

                    zones.push(ZoneState {
                        radius: next_r,
                        dr: dr_grid,
                        temp: world.temp_init,
                        temp_old: world.temp_init,
                        mass_total: d_m,
                        mass_rock: m_rock,
                        mass_rock_init: m_rock,
                        mass_h2os: m_h2os,
                        mass_adhs: m_adhs,
                        mass_h2ol: 0.0,
                        mass_nh3l: 0.0,
                        energy_total: d_e,
                        porosity: world.por_init,
                        pressure: 0.0,
                        x_hydr: world.hydr_init,
                        x_hydr_old: world.hydr_init,
                        kappa: 0.0,
                        nusselt: 1.0,
                    });
                    current_r = next_r;
                }

                let n_orb = if world.orb_a_init > 0.0 {
                    (GCGS * self.saturn.mass / world.orb_a_init.powi(3)).sqrt()
                } else {
                    0.0
                };

                WorldState {
                    name: world.name.clone(),
                    a_orb: world.orb_a_init,
                    e_orb: world.orb_e_init,
                    i_orb: world.orb_i_init,
                    obl: world.orb_o_init,
                    n_orb,
                    lambda: 0.32 + (idx as f64) * 0.77,
                    omega: 0.58 + (idx as f64) * 0.27,
                    spin: n_orb,
                    k2: 0.001,
                    q_tide: 1.0e10,
                    moi: 0.4,
                    zones,
                    res_acct_for: vec![0.; self.n_moons()],
                }
            })
            .collect();

        if self.saturn.mass > 0.0 {
            for s in [
                "Primary",
                "Resonances",
                "ResAcctFor",
                "PCapture",
                "icydwarf_outputs_1",
            ] {
                let _ = create_output(output_path.clone(), format!("{}.txt", s));
            }
        }

        if self.housekeeping.recover {
            return;
        }
        //-------------------------------------------------------------------
        //                       Initialize time loop
        //-------------------------------------------------------------------
        let mut real_time = if !self.housekeeping.recover {
            t_form_min
        } else {
            0.0 // trecover would be here
        };
        real_time -= dtime;

        for _itime in 0..=n_time {
            real_time += dtime;
            let _q_prim = self.saturn.tidal_q.q_prim(&self.worlds, real_time);

            // TODO: Call Orbit and Thermal logic
        }
    }
}

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
