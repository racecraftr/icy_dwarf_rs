use faer::{Mat, linalg::solvers::{DenseSolveCore, Svd}, prelude::Solve};
use num::{complex::Complex64};

use crate::{
    consts::{EROCK_A, EROCK_C, EROCK_D, EROCK_F},
    input::ParsedInput,
};

impl ParsedInput {
    pub fn thermal() {
        let tilt_scale: f64 = 1.0;
    }
}

fn heat_rock(T: f64) -> f64 {
    if T > 1000.0 {
        EROCK_A * 275.0 * 275.0
            + (1000.0 - 275.0) * (EROCK_C + EROCK_D * 1000.0)
            + (EROCK_F) * (T - 1000.0)
    } else if T > 275.0 {
        EROCK_A * 275.0 * 275.0 + (T - 275.0) * (EROCK_C + EROCK_C * T)
    } else {
        EROCK_A * T * T
    }
}

fn svd(mat: &Vec<Vec<f64>>) -> Option<Svd<f64>>{
    let rows = mat.len();
    let cols = mat[0].len();
    Mat::from_fn(
        rows, cols, |i, j| mat[i][j])
    .svd().ok()
}

/// Gauss jordan elimination on a complex matrix. 
/// Given a matrix and a vector, 
/// Uses full pivoting to solve Ax = b.
/// Returns both x and the inverse of A. 
fn gauss_jordan(mat: &Vec<Vec<Complex64>>, b: &Vec<Complex64>) -> Option<(Mat<Complex64>, Mat<Complex64>)> {
    let rows = mat.len();
    let cols = mat[0].len();
    if rows != cols || b.len() != cols {
        return None; // avoids unneccessary computation 
    }
    let A = Mat::from_fn(
        rows, cols, |i, j| mat[i][j]);
    if A.determinant() == Complex64::from(0.0) { // matrix is not invertible
        return None;
    }
    let b = Mat::from_fn(b.len(), 1, |i, _| b[i]);
    let lu = A.full_piv_lu();
    let x = lu.solve(b);
    let a_inv = lu.inverse();
    Some((x, a_inv))
}



#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_svd() {
        let mat: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 0.0, 0.0, 2.0],
            vec![0.0, 0.0, 3.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0, 0.0, 0.0],
        ];
        let Some(svd) = svd(&mat) else {
            assert!(false);
            return;
        };
        println!("S = {:?}", svd.S()); // S represents Sigma (the diagonal entries)
        println!("U = {:?}", svd.U()); // U represents U
        println!("V = {:?}", svd.V().transpose()); // V represents (V^*)^*
    }
}
