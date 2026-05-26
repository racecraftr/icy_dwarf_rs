use crate::{input::IcyDwarfInput, thermal::ThermalOut};

const N_THERM: usize = 14;
const N_ORBIT: usize = 12;
const N_CRK_STRS: usize = 12;
const N_HEAT: usize = 10;
const N_REBOUND: usize = 7;

impl IcyDwarfInput {
    pub fn recover(&self, output_folder: String) {
        let nr = self.grid.n_zones;
        let mut thermal_out: ThermalOut = ThermalOut::default();
    }
}
