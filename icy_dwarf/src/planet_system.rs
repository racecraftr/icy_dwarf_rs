use crate::{create_output, input::ParsedInput};

// pub fn planet_system(parsed: &ParsedInput) {
//     let base_mtx = vec![vec![0; parsed.grid.n_zones]; parsed.worlds.len()];
//     let base_vec = vec![0; parsed.worlds.len()];
// }

impl ParsedInput {
    pub fn planet_sytem(&self, output_path: Option<String>) {
        // the input contains all the planets themselves. So not that hard.
        // if self.core_crack.incl_thermal_mismatch {}

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

    pub fn recover(&mut self) {}
}
