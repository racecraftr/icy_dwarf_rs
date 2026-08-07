# icy_dwarf_rs

Rust rewrite of Dr. Marc Neveu's (mars.f.neveu@nasa.gov) IcyDwarf program 
originally written in C. 

## Why?

Rust allows for a lot of things that C doesn't. 
- Cross-platform compatability (a main issue in the original codebase)
- Memory safety
- Better library logic - much less messier than C

## Roadmap 
- [ ] Implement multithreading for each planet when doing thermal calculations
- [ ] Create a gui for the program

## Installation & Execution

**Prerequisites**: 
- You need [`R`](https://cran.r-project.org/mirrors.html) installed on your machine with the [`CHNOSZ`](https://chnosz.net/) package installed globally. 
- The [`CHNOSZ`](https://chnosz.net/) package is written in fortran, so `libgfortran` must be installed. 
- The [`IPHREEQC`](https://www.usgs.gov/software/phreeqc-version-3) library installed in order to run water rock simulations.
- Given the project *is* written in rust, make sure that [`rustc` and `cargo`](https://rust-lang.org/tools/install/) are installed on your machine. 
- [`just`](https://just.systems/man/en/) should be installed to check if all the above packages are installed, as well as build and run the program.



Note that you 
may have to set the `$RHOME` variable on your machine differently to what it is configured to. 
If you are using **MacOS/Linux**, add the following 
to your `.bashrc`, `.zshrc`, or any shell configuration file you use:
  
```sh
export R_HOME=$(R RHOME)
```

For the `CHNOSZ` package, add this to your shell config file if needed:

```sh
alias gfortran=libgfortran
```

This will let you install the `CHNOSZ` package. run `R` in your terminal, and in the REPL, run 

```r
install.packages("CHNOSZ")
```

From here, installing and running the program is easy. 

**To compile and run:**:

```sh 
git clone https://www.github/racecraftr/icy_dwarf_rs.git
cd icy_dwarf_rs/icy_dwarf
just check # to check if all dependencies are installed 
just build # to build the program! 
cargo run [PATH-TO-INPUT-TOML]
```

<!--***Hold on, what is `build.rs`?***  
Because Rust uses external C/C++ libraries (IPHREEQC), a build script is used to compile them with the rust code. -->

## Input files

IcyDwarf uses TOML configuration files to specify simulation parameters. The root configuration object contains several primary tables mapping to physical properties, orbital parameters, and simulation settings.

### `[housekeeping]`
Configuration flags for logging and recovery.
* `warnings` (boolean): Display simulation warnings.
* `recover` (boolean): Recover simulation state from previous outputs.

### `[grid]`
Spatial and temporal simulation mesh configuration.
* `n_zones` (integer): Number of radial zones per world.
* `time_step` (float): Time step length in years.
* `speedup` (float): Speedup multiplier factor.
* `time_total` (float): Total simulation duration in megayears.
* `output_every` (float): Output generation frequency in megayears.

### `[primary_world]`
Physical properties of the central body.
* `mass` (float): Mass in grams.
* `rad` (float): Volumetric radius in centimeters.
* `moi_coef` (float): Dimensionless moment of inertia coefficient.
* `k2`, `j2`, `j4` (float): Second Love number and zonal harmonic coefficients.
* `tidal_resonant` (boolean): Enable tidal resonant interactions.
* `spin_period` (float): Rotation period in hours.
* **`[primary_world.tidal_q]`**: Tidal dissipation parameters.
  * `init`, `today` (float): Initial and current tidal dissipation factors.
  * `mode` (string/integer): Dissipation model (`"linear"`/0, `"decay"`/1, `"expchange"`/2).
* **`[primary_world.ring]`**: Ring system properties.
  * `mass` (float): Ring mass in grams.
  * `inner`, `outer` (float): Inner and outer ring radii in centimeters.

### `[world_spec]`
Shared physical material parameters across bodies.
* `rho_rock_dry`, `rho_rock_hydr` (float): Dry and hydrated rock densities (g/cm³).
* `chondrite` (string/integer/boolean): Chondrite composition type (`"CI"`/0/false or `"CO"`/1/true).
* `rhelogy` (string/integer): Tidal rheology model (`"maxwell"`/2, `"burgers"`/3, `"andr"`/4, `"suncoop"`/5).
* `ecc_model` (string/integer): Eccentricity model (`"e2"`/0, `"e10cpl"`/1, `"e10ctl"`/2).
* `tidal_heating` (boolean): Enable tidal heating calculation.
* `lookup_tbl` (array of floats): Lookup table for material properties.

### `[[worlds]]`
Array of tables defining individual secondary icy moons.
* `name` (string): World name.
* `planetary_rad` (float): Radius in kilometers.
* `planetary_dens` (float): Mean bulk density in g/cm³.
* `temp_surf`, `temp_init` (float): Surface and initial interior temperatures (K).
* `t_form` (float): Formation time post-solar system formation (Myr).
* `from_ring` (boolean): Formed from a planetary ring.
* `ammonia` (float): Initial ammonia mass fraction.
* `briny` (boolean): Indicates presence of salt in ocean.
* `hydr_init` (float): Initial degree of hydration.
* `hydrate` (boolean): Indicates gas hydrate presence.
* `por_init` (float): Initial rock matrix porosity.
* `rock_frac` (float): Rock mass fraction.
* `rock_h20` (float): Initial rock water content ratio.
* `start_diff` (boolean): Start body in a differentiated state.
* `orb_a_init` (float): Initial semi-major axis (km).
* `orb_e_init`, `orb_i_init`, `orb_o_init` (float): Initial eccentricity, inclination (degrees), and obliquity (degrees).
* `orb_can_change` (boolean): Allow orbital migration.
* `retrograde` (boolean): Indicates retrograde orbit.
* `t_reslock` (float): Resonance locking duration (Myr).

### `[core_crack]`
Configuration for core cracking mechanisms.
* `incl_therm_mismatch`, `incl_pore`, `incl_hydr`, `incl_dissol` (boolean): Include specific stress mechanisms (thermal mismatch, pore pressure, hydration volume expansion, mineral dissolution).
* **`[core_crack.dissol]`**: Dissolution configurations.
  * `of_silica`, `of_serp`, `of_carb` (boolean): Include silica, serpentine, or carbonate dissolution.

### `[subroutines]`
Toggles and parameters for various simulation subroutines.
* `run_therm`, `gen_crack_core`, `gen_water_ab`, `gen_crack_sp`, `run_geo`, `run_comp`, `run_cryo` (boolean): Flags to enable specific models and outputs.
* **`[subroutines.geo]`**: Geochemical routine bounds. Contains `temp`, `pressure`, `pe`, `wr_ratio` tables, each with `min`, `max`, and `step` (float).
* **`[subroutines.cryo]`**: Cryolava routine settings.
  * `after` (integer): Start step for cryolava calculations.
  * `min_temp_chnosz` (float): Minimum temperature for CHNOSZ in Kelvin.

## Libraries used

The libraries used for this project are listed below: 
- `serde`, `serde_repr`, and `toml` for the parsing of TOML 
- `faer` for linear algebra
- `clap` for command-line argument parsing 
- `nalgebra` and `num` for complex values and useful num traits 
- `itertools` for useful functions on iterators
- `extendr` for running `R` code in rust 

Thank you to these libraries for helping me out!
