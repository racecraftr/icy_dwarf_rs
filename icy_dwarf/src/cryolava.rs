use crate::{
    consts::GRAM,
    input::{ParsedInput, ThermalOut},
};

pub const N_SPECIES: i32 = 10;

impl ParsedInput {
    pub fn cryolava(&self, thermal_out: &Vec<Vec<ThermalOut>>) {
        let r_p = self.worlds[0].r_p;
        let t = self.t_cryo();
        let chnoz_t_min = self.subroutines.cryo.min_temp_chnosz;
        let rho_dry = self.world_spec.rho_rock_dry;
        let rho_hydr = self.world_spec.rho_rock_hydr;
        let x_hydr: Vec<Vec<f64>> = vec![]; // TODO: initialize this variable (somehow)

        todo!()
    }

    fn calculate_mass_liquid(&self, thermal_out: &[Vec<ThermalOut>]) -> f64 {
        let t = self.grid.output_time_step();
        thermal_out.iter().map(|v| v[t].mass_water * GRAM).sum()
    }

    fn calculate_seafloor(&self, thermal_out: &Vec<Vec<ThermalOut>>) -> i32 {
        let t = self.t_cryo() as usize;
        if t >= self.grid.output_time_step() {
            println!("[[CALCULATE SEAFLOOR]] :: t >= NT");
        }
        let mut r = 0;
        while r < self.grid.n_zones {
            if thermal_out[r][t].mass_rock <= 0.0 {
                return r as i32 - 1;
            }
            r += 1;
        }
        self.grid.n_zones as i32 - 1
    }
}

fn f(p: f64, mliq: f64, abundances: &[f64], k_rxn: &[f64], x: f64) -> f64 {
    let mut f_x = p * mliq * k_rxn.iter().map(|&n| n + x).product::<f64>();
    for (i, &abundance) in abundances.iter().enumerate() {
        let lhs = k_rxn
            .iter()
            .enumerate()
            .filter_map(|(j, &n)| if j != i { Some(n + x) } else { None })
            .product::<f64>();
        f_x -= abundance * lhs;
    }
    f_x
}

fn f_prime(p: f64, mliq: f64, abundances: &[f64], k_rxn: &[f64], x: f64) -> Option<f64> {
    if abundances.len() != k_rxn.len() {
        return None;
    }
    let mut f_prime_x = p * mliq;
    let mut rhs = 0.0;
    for i in 0..k_rxn.len() {
        rhs += k_rxn
            .iter()
            .enumerate()
            .filter_map(|(j, &n)| if j != i { Some(n + x) } else { None })
            .product::<f64>();
    }

    f_prime_x *= rhs;
    let mut lhs = 0.0;
    for (i, abundance) in abundances.iter().enumerate() {
        let mut lhs2 = 0.0;
        for j in (0..k_rxn.len()).filter(|&j| j != i) {
            lhs2 += k_rxn
                .iter()
                .enumerate()
                .filter_map(|(u, &n)| if u != j && u != i { Some(n + x) } else { None })
                .product::<f64>();
        }
        lhs += lhs2 * abundance;
    }
    Some(f_prime_x - lhs)
}
