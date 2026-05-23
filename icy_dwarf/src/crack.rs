pub struct Data {
    pub atp: Vec<Vec<f64>>,
    pub alpha: Vec<Vec<f64>>,
    pub beta: Vec<Vec<f64>>,
    pub chrysopile: Vec<Vec<f64>>,
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

pub fn read_data(dat_folder: &String) {}
