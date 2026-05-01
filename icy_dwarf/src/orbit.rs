pub fn laplace_coef(alpha: f64, j: f64, s: f64) -> f64 {
    let mut b_lap_j = 1_f64;
    let mut temp = 1_f64;
    for m in 0..200 {
        let m = m as f64;
        temp = temp * (s + m) / (1. + m) * (s + j + m) / (j + 1. + m) * alpha.powi(2);
        b_lap_j += temp * (j + m + 1.) * 2.;
        if temp < 1.0e-6 {
            break;
        }
    }
    b_lap_j
        * alpha.powf(j)
        * (0..200)
            .map(|m| (s + m as f64) / (1. + m as f64))
            .product::<f64>()
        * 2.
}

pub fn d_laplace_coef(alpha: f64, j: f64, s: f64) -> f64 {
    let mut db_lap_j = 1_f64;
    let mut temp = 1_f64;
    for m in 0..200 {
        let m = m as f64;
        temp = temp * (s + m) / (1. + m) * (s + j + m) / (j + 1. + m) * alpha.powi(2);
        db_lap_j += temp;
        if temp * (j + m + 1.) * 2. < 1.0e-6 {
            break;
        }
    }
    db_lap_j
        * alpha.powf(j - 1.)
        * (0..200)
            .map(|m| (s + m as f64) / (1. + m as f64))
            .product::<f64>()
        * 2.
}

pub fn d2_laplace_coef(alpha: f64, j: f64, s: f64) -> f64 {
    todo!(); //TODO: translate the code 
    // let mut d2b_lap_j = j * (j - 1.);
    // let mut temp = 1.;
    // for m in 0..200 {
    //     let m = m as f64;
    //     temp = temp * (s + m) / (1.0 + m) * (s + j + m) / (j + 1. + m) * alpha.powi(2);
    //     d2b_lap_j += temp * ((j + m + 1.) * 2. - 1.) * ((j + m + 1.) * 2. - 1.);
    // }
    // d2b_lap_j
}