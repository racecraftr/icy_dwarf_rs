use std::{array, ops::Deref};

use faer::{
    Mat,
    col::generic::Col,
    prelude::Solve,
    sparse::{
        SparseColMat, Triplet,
        csc_numeric::{Own, generic::SparseColMat},
        linalg::solvers::{Lu, SymbolicLu},
    },
    traits::ComplexField,
};
use itertools::{Itertools, multizip};
use num::{
    complex::Complex64,
    traits::{ConstZero, Inv},
};

use crate::{input::IcyDwarfInput, traits::float_traits::FloatExt};
type C = Complex64;

const SH_TERMS: usize = 500;
const BASE_ARR: [Complex64; SH_TERMS] = [Complex64::ZERO; SH_TERMS];

impl IcyDwarfInput {
    pub fn tropf(
        &self,
        til_cesq: f64,
        til_t: f64,
        tilom: Complex64,
        s: usize,
        pn_fsf_amp: f64,
    ) -> (Complex64, f64, f64) {
        let diss_type = DissType::KinPE; // note that this is always true.
        let nf = 2_usize;
        let til_om = 1_f64;

        let tidal_pot_shc: [_; SH_TERMS] = array::from_fn(|i| {
            if i == nf - s {
                Complex64::I * 0.5 / pn_fsf_amp
            } else {
                Complex64::ZERO
            }
        });

        let vert_struct_shc = BASE_ARR;
        let sh_deg_arr = BASE_ARR;
        let divergence_sch = BASE_ARR;
        let curl_sch = BASE_ARR;
        let grav_potential = BASE_ARR;
        let potential_dissipation = til_t.inv();
        let slowness = Complex64::new(1., potential_dissipation / til_om) / til_cesq;

        let attn_hori = til_t.inv();
        let attn_vert = attn_hori;

        let n_vec: [_; SH_TERMS] = array::from_fn(|i| s + i);

        let diss_d = n_vec.map(|x| (attn_hori * (x + s) as f64).instead_of(0., f64::MIN_POSITIVE));
        let diss_r = n_vec.map(|x| (attn_vert * (x + s) as f64).instead_of(0., f64::MIN_POSITIVE));

        let l_alpha_d = arr_to_diag(&diss_d);
        let l_alpha_r = arr_to_diag(&diss_r);
        let lv_values = [slowness; SH_TERMS];
        let lv = arr_to_diag(&lv_values);

        let mut l_vec = n_vec.map(|n| {
            let n = n as f64;
            -n * (n + 1.)
        });
        if l_vec[0] == 0. {
            l_vec[0] = f64::MIN_POSITIVE
        }
        let ll = arr_to_diag(&l_vec);
        let la = arr_to_diag(&array::from_fn::<_, SH_TERMS, _>(|i| {
            (tilom + Complex64::I * diss_d[i]) * l_vec[i] - s as f64 * til_om
        }));

        let mut lc_values = [Complex64::ZERO; 2 * SH_TERMS - 2];
        lc_values[0] = Complex64::from({
            let n = n_vec[1] as f64;
            -(n * n - 1.) * (n + s as f64) / (2. * n + 1.)
        });
        for i in 0..lc_values.len() {
            let k = (i as f64 / 2.).ceil() as usize;
        }

        todo!()
    }
}

fn arr_to_diag<T, const N: usize>(arr: &[T; N]) -> SparseColMat<usize, T>
where
    T: ComplexField + Copy,
{
    let triplets: [_; N] = arr
        .iter()
        .enumerate()
        .map(|(i, &n)| Triplet::new(i, i, n))
        .collect_array()
        .unwrap();
    SparseColMat::try_new_from_triplets(N, N, &triplets).unwrap()
}

#[repr(u8)]
pub enum DissType {
    Kinetic,
    PE,
    KinPE,
}

fn ratio_factorials(n: usize, s: usize) -> f64 {
    ((n - s + 1)..=(n + s)).product::<usize>() as f64
}

/// uses the [`saer`] library to perform the Bi-conjugate Gradient Stabilized
/// method on a matrix to solve the equation Ax = b.
fn bicgstab(a: &[Vec<Complex64>], b: &[Complex64]) -> Vec<Complex64> {
    let rows = a.len();
    let cols = a[0].len();

    let triplets = a
        .iter()
        .enumerate()
        .flat_map(|(r, row)| {
            row.iter()
                .cloned()
                .enumerate()
                .filter(|(_, val)| *val != Complex64::from(0.0))
                .map(move |(c, val)| Triplet::new(r, c, val))
        })
        .collect::<Vec<_>>();

    let b = Col::from_fn(b.len(), |i| b[i]);
    let a_sparse = SparseColMat::try_new_from_triplets(rows, cols, &triplets).unwrap();

    let symb = SymbolicLu::try_new(a_sparse.symbolic()).unwrap();
    let lu = Lu::try_new_with_symbolic(symb, a_sparse.as_ref()).unwrap();
    let x = lu.solve(&b);

    x.iter().cloned().collect::<Vec<_>>()
}

/// Finds the complex eigenvalues of a real matrix.
fn eigen(mtx: &[Vec<f64>]) -> Option<Vec<Complex64>> {
    let rows = mtx.len();
    let cols = mtx[0].len();

    let mat = Mat::from_fn(rows, cols, |i, j| mtx[i][j]);
    mat.eigenvalues().ok() // that's it, really
}

fn globe_time_average(s_coefs: &[C], t_coefs: &[C], s: i32, n_vec: &[i32]) -> Vec<f64> {
    multizip((s_coefs, t_coefs, n_vec))
        .map(|(&sc, &tc, &n)| {
            let sc_c = sc.conjugate();
            let tc_c = tc.conjugate();
            (sc * tc_c + sc_c * tc).re / (2. * n as f64 + 1.)
                * ratio_factorials(n as usize, s as usize)
        })
        .collect()
}

#[cfg(test)]
mod eigen_tests {
    use super::*;
    #[test]
    fn eigen_test_1() {
        let diag_mtx = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, -2.0, 0.0],
            vec![0.0, 0.0, 5.0],
        ];
        let Some(mut eigenvalues) = eigen(&diag_mtx) else {
            panic!("Unable to find eigenvalues of matrix");
        };
        eigenvalues.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
        assert!(eigenvalues == vec![C::from(1.0), C::from(-2.0), C::from(5.0)])
    }

    #[test]
    fn eigen_test_2() {
        let diag_mtx = vec![
            vec![1.0, 3.0, -4.0],
            vec![0.0, -2.0, 9.5],
            vec![0.0, 0.0, 5.0],
        ];
        let Some(mut eigenvalues) = eigen(&diag_mtx) else {
            panic!("Unable to find eigenvalues of matrix");
        };
        eigenvalues.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
        assert!(eigenvalues == vec![C::from(1.0), C::from(-2.0), C::from(5.0)])
    }

    #[test]
    fn eigen_test_3() {
        // Block diagonal matrix with complex eigenvalues 1 +/- i and 2 +/- 3i
        let complex_mtx = vec![
            vec![1.0, -1.0, 0.0, 0.0],
            vec![1.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 2.0, -3.0],
            vec![0.0, 0.0, 3.0, 2.0],
        ];
        let Some(mut eigenvalues) = eigen(&complex_mtx) else {
            panic!("Unable to find eigenvalues of matrix");
        };
        // Sort by real part then imaginary part to ensure a deterministic order
        eigenvalues.sort_by(|a, b| a.re.total_cmp(&b.re).then(a.im.total_cmp(&b.im)));
        println!("eigenvalues are {:?}", &eigenvalues);
        let targets = [
            C::new(1.0, -1.0),
            C::new(1.0, 1.0),
            C::new(2.0, -3.0),
            C::new(2.0, 3.0),
        ];

        // check mismatch by a small amount
        for (a, b) in eigenvalues.iter().zip(targets.iter()) {
            assert!(
                (a.re - b.re).abs() < 1e-12,
                "Real part mismatch: {:?} != {:?}",
                a,
                b
            );
            assert!(
                (a.im - b.im).abs() < 1e-12,
                "Imaginary part mismatch: {:?} != {:?}",
                a,
                b
            );
        }
    }
}

#[cfg(test)]
mod bicgstab_tests {
    use super::*;

    #[test]
    fn test() {
        let a = [[2., 1., -1.], [-3., -1., 2.], [-2., 1., 2.]]
            .iter()
            .map(|v| v.iter().map(C::from).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let b = [-5., 9., 6.].iter().map(C::from).collect::<Vec<_>>();
        let x = bicgstab(&a, &b);
        println!("{:?}", &x);
        let res = x
            .iter()
            .zip([1., -2., 5.])
            .all(|(a, b)| (a.re - b).abs() < 1e-12);
        assert!(res);
    }
}

#[cfg(test)]
mod ratio_factorials_tests {
    use super::*;

    #[test]
    fn test_ratio_factorials() {
        // (2+2)! / (2-2)! = 24 / 1 = 24
        assert_eq!(ratio_factorials(2, 2), 24.0);
        // (3+1)! / (3-1)! = 24 / 2 = 12
        assert_eq!(ratio_factorials(3, 1), 12.0);
        // (3+2)! / (3-2)! = 120 / 1 = 120
        assert_eq!(ratio_factorials(3, 2), 120.0);
    }
}

#[cfg(test)]
mod globe_time_average_tests {
    use super::*;

    #[test]
    fn test_globe_time_average() {
        let s_coefs = [C::new(1.0, 1.0)];
        let t_coefs = [C::new(2.0, -1.0)];
        // sc = 1 + i, tc = 2 - i
        // sc_c = 1 - i, tc_c = 2 + i
        // sc * tc_c = (1 + i)(2 + i) = 2 + i + 2i - 1 = 1 + 3i
        // sc_c * tc = (1 - i)(2 - i) = 2 - i - 2i - 1 = 1 - 3i
        // sc * tc_c + sc_c * tc = 2.0 (real part)
        // For n = 2, s = 2:
        // (2.0) / (2 * 2 + 1) * ratio_factorials(2, 2)
        // = 2.0 / 5.0 * 24.0 = 9.6
        let res = globe_time_average(&s_coefs, &t_coefs, 2, &[2]);
        assert_eq!(res.len(), 1);
        assert!((res[0] - 9.6).abs() < 1e-12);
    }
}
