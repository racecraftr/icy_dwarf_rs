use faer::{
    Mat,
    col::generic::Col,
    prelude::Solve,
    sparse::{
        Triplet,
        csc_numeric::generic::SparseColMat,
        linalg::solvers::{Lu, SymbolicLu},
    },
};
use num::complex::{Complex64, ComplexFloat};

#[repr(u8)]
pub enum DissType {
    Kinetic,
    PE,
    KinPE,
}

type C = Complex64;

pub fn tropf(
    tilceq: f64,
    til_T: f64,
    diss_type: DissType,
    til_om: C,
    nF: i32,
    s: i32,
    pn_fsf: f64,
    kn_fsf: C,
    kLovenF: C,
) {
    let til_omega = C::from(1.0);
}

fn dot_prod(v1: &[C], v2: &[C]) -> C {
    v1.iter().zip(v2.iter()).map(|(v1, v2)| v1 * v2).sum()
}

fn vec_norm(v: &[C]) -> f64 {
    v.iter()
        .map(|c| c * c.conj())
        .sum::<Complex64>()
        .abs()
        .sqrt()
}

fn ratio_factorials(n: usize, s: usize) -> f64 {
    ((n - s + 1)..(n + s)).product::<usize>() as f64
}

/// uses the [`saer`] library to perform the Bi-conjugate Gradient Stabilized
/// method on a matrix to solve the equation Ax = b.
/// Much of the work is done in conversion and setup,
/// but the actual map is being done by
fn bicgstab(a: &Vec<Vec<Complex64>>, b: &Vec<Complex64>) -> Vec<Complex64> {
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
fn eigen(mtx: &Vec<Vec<f64>>) -> Option<Vec<Complex64>> {
    let rows = mtx.len();
    let cols = mtx[0].len();

    let mat = Mat::from_fn(rows, cols, |i, j| mtx[i][j]);
    mat.eigenvalues().ok() // that's it, really
}

