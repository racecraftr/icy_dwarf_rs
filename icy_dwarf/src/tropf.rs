use std::array;

use faer::{
    Col, Mat, Scale,
    linalg::zip::IntoView,
    prelude::Solve,
    sparse::{
        SparseColMat, Triplet,
        linalg::solvers::{Lu, SymbolicLu},
    },
    traits::ComplexField,
};
use itertools::multizip;
use num::{complex::Complex64, traits::Inv};

use crate::{
    input::IcyDwarfInput,
    traits::{float_traits::FloatExt, num_ext::NumExt},
};
type C = Complex64;

const C0: Complex64 = Complex64::ZERO;
const I: Complex64 = Complex64::I;

const SH_TERMS: usize = 500;
const BASE_ARR: [Complex64; SH_TERMS] = [Complex64::ZERO; SH_TERMS];

impl IcyDwarfInput {
    /// Port of the Tropf function from MatLab to Rust.
    /// Why? Well, MatLab ain't free...
    /// Combines the TROPF and tropf functions from the original C code.
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

        let source_sink_shc = BASE_ARR;
        let sh_deg_arr = BASE_ARR;
        let divergence_shc = BASE_ARR;
        let curl_shc = BASE_ARR;
        let grav_potential = BASE_ARR;
        let potential_dissipation = til_t.inv();
        let slowness = Complex64::new(1., potential_dissipation / til_om) / til_cesq;

        let attn_hori = til_t.inv();
        let attn_vert = attn_hori;

        let n_vec: [_; SH_TERMS] = array::from_fn(|i| s + i);

        let diss_d =
            n_vec.map(|x| (attn_hori * (x + s) as f64 + C0).instead_of(C0, f64::MIN_POSITIVE + C0));
        let diss_r =
            n_vec.map(|x| (attn_vert * (x + s) as f64 + C0).instead_of(C0, f64::MIN_POSITIVE + C0));

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
        let ll = arr_to_diag(&l_vec.map(|n| n + C0));
        let la = arr_to_diag(&array::from_fn::<_, SH_TERMS, _>(|i| {
            (tilom + Complex64::I * diss_d[i]) * l_vec[i] - s as f64 * til_om
        }));

        let mut lc_values = [Complex64::ZERO; 2 * SH_TERMS - 2];
        let mut tri_diag_indices = [(0_usize, 0_usize); 2 * SH_TERMS - 2];
        lc_values[0] = Complex64::from({
            let n = n_vec[1] as f64;
            -(n * n - 1.) * (n + s as f64) / (2. * n + 1.)
        });
        for i in (1..lc_values.len()).step_by(2) {
            let k = (i as f64 / 2.).ceil() as usize;
            let n = n_vec[k - 1] as f64;
            lc_values[i] = til_om * -n * (n + 1.) * (n + 2.) / (2. * n + 1.) + C0;
            tri_diag_indices[i] = (k, k - 1);
            if i < 2 * SH_TERMS - 3 {
                let n = n_vec[k + 1] as f64;
                lc_values[i + 1] = til_om * -n * (n + 1.) * (n + 1.) / (2. * n + 1.) + C0;
                tri_diag_indices[i + 1] = (k, k + 1);
            }
        }

        let triplets = lc_values
            .iter()
            .zip(tri_diag_indices)
            .map(|(v, (i, j))| Triplet::new(i, j, *v))
            .collect::<Vec<_>>();
        let lc = SparseColMat::try_new_from_triplets(2 * SH_TERMS - 2, 2 * SH_TERMS - 2, &triplets)
            .unwrap();

        let ld = arr_to_diag(&array::from_fn::<_, SH_TERMS, _>(|i| {
            (tilom + Complex64::I * diss_d[i]) * l_vec[i] - (s as f64) * tilom
                + 1. / tilom * l_vec[i] / lv_values[i] * l_vec[i]
        }));
        let l_vi = arr_to_diag(&[slowness.inv(); SH_TERMS]);
        let l_li = arr_to_diag(&l_vec.map(|l| l.inv() + C0));
        let l_bi = arr_to_diag(&array::from_fn::<_, SH_TERMS, _>(|i| {
            (tilom + diss_d[i] * I).inv() * l_vec[i] - s as f64 * til_om
        }));
        let eye = arr_to_diag(&[1. + C0; SH_TERMS]);

        let rho_ratio = 0.5;
        let iln_vals = arr_to_diag(&n_vec.map(|n| {
            let n = n as f64;
            -3. / (2. * n + 1.) * rho_ratio + C0
        }));

        // solvemethod is always 1, so we solve for pns.
        // Faer lets us do matrix multiplication really easily.

        let l_tilp =
            &l_li * (&la - &lc * &l_bi * &lc) * Scale::from_ref(&tilom) * &l_li * &lv + &eye;

        let tidal_pot_shc_col = arr_to_col(&tidal_pot_shc);
        let source_sink_shc_col = arr_to_col(&source_sink_shc);
        let divergence_shc_col = arr_to_col(&divergence_shc);
        let curl_shc_col = arr_to_col(&curl_shc);

        let q_tilp = (-&l_li) * (&la - &lc * &l_bi * &lc) * &l_li * &source_sink_shc_col
            + &l_li * &divergence_shc_col
            + &l_li * &lc * &l_bi * &curl_shc_col;

        let lhs_p = &q_tilp + &tidal_pot_shc_col;

        let pns = solve(&(&l_tilp + &iln_vals), &lhs_p);
        let dns = &l_li * &(Scale::from_ref(&tilom) * &lv * &pns + &source_sink_shc_col);

        let rns = &(-&l_bi) * (&lc * &dns + &curl_shc_col);

        let globe_time_average = |s_c: &Col<Complex64>, t_c: &Col<Complex64>| -> Col<Complex64> {
            Col::from_iter(
                globe_time_average(
                    s_c.into_view().as_slice(),
                    t_c.into_view().as_slice(),
                    s,
                    &n_vec,
                )
                .iter()
                .map(|n| Complex64::from(n)),
            )
        };

        let cal_wns = {
            let cal_wns_1 = Scale::from_ref(&(-I)) * &tidal_pot_shc_col;
            let cal_wns_2 = &pns * Scale::from_ref(&(-I * tilom));
            let cal_wns_3 = &lv * &cal_wns_2;

            let cal_wns_temp = globe_time_average(&cal_wns_1, &cal_wns_3);

            let cal_wns_1 = Scale::from_ref(&I) * (&tidal_pot_shc_col - &pns);

            cal_wns_temp + cal_wns_1
        };

        let cal_dns = Scale::from_ref(&-0.5.as_cplx())
            * (globe_time_average(&dns, &(&ll * &l_alpha_d * &dns))
                + globe_time_average(&(&l_alpha_d * &dns), &(&ll * &dns))
                + {
                    let rns = Scale::from_ref(&-I) * &rns;
                    globe_time_average(&rns, &(&ll * &l_alpha_r * &rns))
                });
        todo!()
    }
}

fn arr_to_col<T, const N: usize>(arr: &[T; N]) -> Col<T>
where
    T: ComplexField + Copy,
{
    Col::from_fn(N, |i| arr[i])
}

fn arr_to_diag<T, const N: usize>(arr: &[T; N]) -> SparseColMat<usize, T>
where
    T: ComplexField + Copy,
{
    let triplets = arr
        .iter()
        .enumerate()
        .map(|(i, &n)| Triplet::new(i, i, n))
        .collect::<Vec<_>>();
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

fn solve(mat: &SparseColMat<usize, Complex64>, b: &Col<Complex64>) -> Col<Complex64> {
    let symb = SymbolicLu::try_new(mat.symbolic()).unwrap();
    let lu = Lu::try_new_with_symbolic(symb, mat.as_ref()).unwrap();
    lu.solve(&b)
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

fn globe_time_average(s_coefs: &[C], t_coefs: &[C], s: usize, n_vec: &[usize]) -> Vec<f64> {
    multizip((s_coefs, t_coefs, n_vec))
        .map(|(&sc, &tc, &n)| {
            let sc_c = sc.conj();
            let tc_c = tc.conj();
            (sc * tc_c + sc_c * tc).re / (2. * n as f64 + 1.)
                * ratio_factorials(n as usize, s as usize)
        })
        .collect()
}

#[cfg(test)]
mod eigen_tests {
    use num::complex::ComplexFloat;

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
