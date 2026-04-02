mod input;
mod planet_system;
mod thermal;
mod tropf;

#[allow(dead_code)]
pub mod consts {
    // -----------------------------------------------------------------
    // PHYSICAL AND MATHEMATICAL CONSTANTS
    // -----------------------------------------------------------------

    /// Gravitational constant (SI)
    pub const G: f64 = 6.67e-11;

    /// Gravitational constant (cgs)
    pub const GCGS: f64 = 6.67e-8;

    /// Universal gas constant (J/(mol K))
    pub const R_G: f64 = 8.3145;

    /// Boltzmann's constant (J/K)
    pub const K_B: f64 = 1.3806502e-23;

    /// Ratio of the circumference of a circle to its radius
    pub const PI_GREEK: f64 = 3.14159265358979323846;

    /// Mass of the Earth (kg)
    pub const M_EARTH: f64 = 5.9721986e24;

    /// Radius of the Earth (m)
    pub const R_EARTH: f64 = 6.3675e6;

    /// AMU = 1/N_Avogadro
    pub const AMU: f64 = 1.66054e-24;

    // -----------------------------------------------------------------
    // UNIT CONVERSION FACTORS
    // -----------------------------------------------------------------

    /// km to m
    pub const KM: f64 = 1.0e3;

    /// cm to m
    pub const CM: f64 = 1.0e-2;

    /// km to cm
    pub const KM2CM: f64 = 1.0e5;

    /// g to kg
    pub const GRAM: f64 = 1.0e-3;

    /// bar to Pa
    pub const BAR: f64 = 1.0e5;

    /// Celsius to Kelvin
    pub const KELVIN: f64 = 273.15;

    /// Gyr to seconds
    pub const GYR2SEC: f64 = 3.15576e16;

    /// Myr to seconds
    pub const MYR2SEC: f64 = 3.15576e13;

    /// MeV to erg
    pub const MEV2ERG: f64 = 1.602e-6;

    /// Pa to barye (cgs)
    pub const PA2BA: f64 = 10.0;

    /// MPa to Pa
    pub const MPA: f64 = 1.0e6;

    // -----------------------------------------------------------------
    // GENERAL PARAMETERS
    // -----------------------------------------------------------------

    pub const RHO_H2OS: f64 = 935.0; // Density of H2O(s) kg/m³ (Feistel & Wagner 2006)
    pub const RHO_H2OL: f64 = 1000.0; // Density of H2O(l) kg/m³
    pub const RHO_ADHS: f64 = 985.0; // Density of ADH(s) kg/m³
    pub const RHO_NH3L: f64 = 740.0; // Density of NH3(l) kg/m³
    pub const XC: f64 = 0.321; // Ammonia content of eutectic H2O-NH3 mixture

    // -----------------------------------------------------------------
    // THERMAL PARAMETERS
    // -----------------------------------------------------------------

    /// Heat of hydration, erg/(g forsterite)
    pub const HHYDR: f64 = 5.75e9; 
    /// Heat capacity of rock below 275 K (cgs)
    pub const EROCK_A: f64 = 1.40e4; 
    /// Heat capacity of rock 275–1000 K, term 1 (cgs)
    pub const EROCK_C: f64 = 6.885e6; 
    /// Heat capacity of rock 275–1000 K, term 2 (cgs)
    pub const EROCK_D: f64 = 2.963636e3; 
    /// Heat capacity of rock above 1000 K (cgs)
    pub const EROCK_F: f64 = 1.20e7; 

    /// Heat capacity of water ice (erg/g/K)
    pub const QH2O: f64 = 7.73e4; 
    /// Heat capacity of ADH ice (erg/g/K)
    pub const QADH: f64 = 1.12e5; 
    /// Heat capacity of liquid water (erg/g/K)
    pub const CH2OL: f64 = 4.1885e7; 
    /// Heat capacity of liquid ammonia (cgs)
    pub const CNH3L: f64 = 4.7e7; 
    /// Latent heat of ADH melting (cgs)
    pub const LADH: f64 = 1.319e9; 
    /// Latent heat of H2O melting (cgs)
    pub const LH2O: f64 = 3.335e9; 
    /// Bulk permeability for D=1m cracks (m²)
    pub const PERMEABILITY: f64 = 1.0e-9; 
    /// Porosity from cracking (dimensionless)
    pub const CRACK_POROSITY: f64 = 0.01; 
    /// Temperature at which differentiation proceeds (K)
    pub const TDIFF: f64 = 140.0; 
    /// Temperature of full silicate hydration (K)
    pub const TDEHYDR_MIN: f64 = 700.0; 
    /// Temperature of full silicate dehydration (K)
    pub const TDEHYDR_MAX: f64 = 850.0; 
    /// Effective thermal conductivity, hydrothermal layer (cgs)
    pub const KAP_HYDRO: f64 = 100.0e5; 
    /// Effective thermal conductivity, convective slush (cgs)
    pub const KAP_SLUSH: f64 = 400.0e5;
    /// Effective thermal conductivity, convective ice (cgs) 
    pub const KAP_ICE_CV: f64 = 150.0e5; 
    /// Thermal conductivity, dry silicate rock (cgs)
    pub const KAPROCK: f64 = 4.2e5; 
    /// Thermal conductivity, hydrated silicate rock (cgs)
    pub const KAPHYDR: f64 = 1.0e5; 
    /// Thermal conductivity of ADH ice (cgs)
    pub const KAPADHS: f64 = 1.2e5; 
    /// Thermal conductivity of liquid water (cgs)    
    pub const KAPH2OL: f64 = 0.61e5; 
    /// Thermal conductivity of liquid ammonia (cgs)
    pub const KAPNH3L: f64 = 0.022e5; 
    /// Average expansivity of water (K⁻¹)
    pub const ALFH2OAVG: f64 = 1.0e-3; 
    /// Memory of old hydration state (0=none, 1=no change)
    pub const F_MEM: f64 = 0.75; 

    // -----------------------------------------------------------------
    // CRACKING PARAMETERS
    // -----------------------------------------------------------------

    pub const E_YOUNG_OLIV: f64 = 200.0e9; // Young's modulus, olivine (Pa) (Christensen 1966)
    pub const E_YOUNG_SERP: f64 = 35.0e9; // Young's modulus, serpentinite (Pa) (Christensen 1966)
    pub const NU_POISSON_OLIV: f64 = 0.25; // Poisson's ratio, olivine (Christensen 1966)
    pub const NU_POISSON_SERP: f64 = 0.35; // Poisson's ratio, serpentinite (Christensen 1966)
    pub const SMALLEST_CRACK_SIZE: f64 = 1.0e-2; // Smallest 1D or 2D crack size (m)

    pub const MU_F_SERP: f64 = 0.4; // Friction coefficient, serpentine (Escartin et al. 1997)
    pub const MU_F_BYERLEE_LOP: f64 = 0.85; // Friction coefficient, olivine <200 MPa (Byerlee 1978)
    pub const MU_F_BYERLEE_HIP: f64 = 0.6; // Friction coefficient, olivine 200–1700 MPa (Byerlee 1978)
    pub const C_F_BYERLEE_HIP: f64 = 50.0e6; // Frictional cohesive strength, olivine 200–1700 MPa (Pa)
    pub const D_FLOW_LAW: f64 = 500.0; // Grain size (microns)

    pub const K_IC_OLIV: f64 = 1.5e6; // Critical stress intensity, olivine (Pa·m^0.5)
    pub const K_IC_SERP: f64 = 0.4e6; // Critical stress intensity, serpentinite (Pa·m^0.5)
    pub const DELTA_ALPHA: f64 = 3.1e-6; // Thermal expansion anisotropy (K⁻¹)
    pub const QGBS: f64 = 3.75e5; // Activation enthalpy, grain boundary sliding (J/mol)
    pub const OMEGA: f64 = 1.23e-29; // Atomic volume (m³)
    pub const D0_DELTAB: f64 = 0.2377; // Grain boundary diffusion coefficient × width (m³/s)
    pub const N_FIT: f64 = 23.0; // Fitting parameter for diff eq (1)
    pub const L_SIZE: f64 = 0.25e-3; // Half grain size (m) (Vance et al. 2007)
    pub const A_VAR_MAX: f64 = 5.0e-5; // Max flaw size search upper bound (m)
    pub const A_MIN: f64 = 1.0e-7; // Minimum flaw size (m)

    pub const ASPECT_RATIO: f64 = 1.0e4; // Aspect ratio (width/length) of 2D water pores

    pub const NU_PROD_SILICA: f64 = 1.0; // Product stoichiometric coefficient, SiO2(s)
    pub const NU_PROD_CHRYSOTILE: f64 = 11.0; // Product stoichiometric coefficient, chrysotile
    pub const NU_PROD_MAGNESITE: f64 = 2.0; // Product stoichiometric coefficient, magnesite
    pub const MU_XU_SILICA: f64 = 1.0; // Q/K exponent, silica
    pub const MU_XU_CHRYSOTILE: f64 = 1.0; // Q/K exponent, chrysotile
    pub const MU_XU_MAGNESITE: f64 = 4.0; // Q/K exponent, magnesite (Pokrovski & Schott 1999)
    pub const EA_SILICA: f64 = 62.9e3; // Activation energy, silica reaction (J/mol)
    pub const EA_CHRYSOTILE: f64 = 70.0e3; // Activation energy, serpentine reaction (J/mol)
    pub const EA_MAGNESITE: f64 = 32.1e3; // Activation energy, carbonate reaction (J/mol)
    pub const MOLAR_VOLUME_SILICA: f64 = 29.0e-6; // Molar volume of silica (m³/mol)
    pub const MOLAR_VOLUME_CHRYSOTILE: f64 = 108.5e-6; // Molar volume of serpentine (m³/mol)
    pub const MOLAR_VOLUME_MAGNESITE: f64 = 28.018e-6; // Molar volume of carbonate (m³/mol)

    // Table sizes

    /// Integration steps
    pub const INT_STEPS: i32 = 10000;

    /// Data points in integral table
    pub const INT_SIZE: i32 = 1000;

    /// Size of square a(deltaT,P) table
    pub const SIZEATP: i32 = 100;

    /// deltaT intervals for a(deltaT,P) (K)
    pub const DELTA_T_STEP: f64 = 20.0;

    /// P intervals for a(deltaT,P) (Pa)
    pub const P_STEP: f64 = 2.5e6;

    /// Temperature step, 261–2241 K (K)
    pub const DELTA_TEMPK: f64 = 20.0;

    /// Pressure step, 0.1–2475.1 bar
    pub const DELTA_P_BAR: f64 = 25.0;

    /// Minimum temperature (K)
    pub const TEMPK_MIN: f64 = 261.0;

    /// Minimum pressure (bar)
    pub const P_BAR_MIN: f64 = 0.1;

    /// Minimum temperature for species (K)
    pub const TEMPK_MIN_SPECIES: f64 = 261.0;

    /// Temperature step for species (K)
    pub const DELTA_TEMPK_SPECIES: f64 = 7.0;

    // -----------------------------------------------------------------
    // WATER-ROCK PARAMETERS
    // -----------------------------------------------------------------

    pub const NVAR: i32 = 1024; // Geochemical variables per PHREEQC simulation
    pub const NAQ: i32 = 257; // Aqueous species (+ physical parameters)
    pub const NGASES: i32 = 15; // Gaseous species
    pub const NMINGAS: i32 = 389; // Minerals and gases
    pub const NELTS: i32 = 31; // 30 elements + 1 extra column

    // -----------------------------------------------------------------
    // ORBITAL EVOLUTION PARAMETERS
    // -----------------------------------------------------------------
    /// Max order to look for resonances
    pub const IJMAX: i32 = 5; 
    /// Minimum eccentricity
    pub const MIN_ECC: f64 = 1.0e-4; 
}

#[allow(dead_code)]
pub struct Thermalout {
    radius: f64,
    temp_k: f64,
    m_rock: f64,
    m_h2os: f64,
    m_nh31: f64,
    nu: f64,
    famor: f64,
    kappa: f64,
    xhydr: f64,
    pore: f64,
    crack: f64,
    w_tide: f64,
}

fn main() {}
