pub mod float_traits {

    pub trait FloatExt {
        fn unwrap_or_nan(self, fallback: Self) -> Self;
        fn max_assign(&mut self, other: Self);
    }

    impl FloatExt for f64 {
        fn unwrap_or_nan(self, fallback: Self) -> Self {
            if self.is_nan() { fallback } else { self }
        }

        fn max_assign(&mut self, other: Self) {
            *self = self.max(other);
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
