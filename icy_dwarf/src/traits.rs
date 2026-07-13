pub mod float_traits {
    use num::complex::Complex64;

    pub trait FloatExt {
        fn unwrap_or_nan(self, fallback: Self) -> Self;
        fn max_assign(&mut self, other: Self);
        fn as_cplx(&self) -> Complex64;
    }

    impl FloatExt for f64 {
        fn unwrap_or_nan(self, fallback: Self) -> Self {
            if self.is_nan() { fallback } else { self }
        }

        fn max_assign(&mut self, other: Self) {
            *self = self.max(other);
        }

        fn as_cplx(&self) -> Complex64 {
            Complex64::from(self)
        }
    }
}

pub mod min_max {

    pub trait MinMax {
        fn min_max(self) -> Self;
    }

    impl MinMax for (f64, f64) {
        fn min_max(self) -> Self {
            let (a, b) = self;
            (a.min(b), a.max(b))
        }
    }
}

pub mod num_ext {
    use num::complex::Complex64;

    pub trait NumExt {
        fn instead_of(&self, cond_val: Self, instead: Self) -> Self;
    }

    impl NumExt for f64 {
        fn instead_of(&self, cond_val: Self, instead: Self) -> Self {
            if *self == cond_val { instead } else { *self }
        }
    }

    impl NumExt for Complex64 {
        fn instead_of(&self, cond_val: Self, instead: Self) -> Self {
            if *self == cond_val { instead } else { *self }
        }
    }
}
