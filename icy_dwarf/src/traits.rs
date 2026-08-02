//! This module provides helper traits for numerical operations.

/// This module defines helper methods for floating-point numbers.
pub mod float_traits {
    use num::complex::Complex64;

    /// Extension trait that adds utility methods to floating-point numbers.
    pub trait FloatExt {
        /// Return the value if it is not NaN. Otherwise, return the fallback value.
        ///
        /// # Parameters
        /// - `fallback`: The default value to return if `self` is NaN.
        fn unwrap_or_nan(self, fallback: Self) -> Self;

        /// Replace the current value with the maximum of the current value and another value.
        ///
        /// # Parameters
        /// - `other`: The value to compare against `self`.
        fn max_assign(&mut self, other: Self);

        /// Convert the floating-point number into a complex number.
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

/// This module defines helper methods for value substitution.
pub mod num_ext {
    use num::complex::Complex64;

    /// Trait to replace a specific value with another value.
    pub trait NumExt {
        /// Return a new value if the current value matches a target value.
        ///
        /// # Parameters
        /// - `cond_val`: The target value to check for equality.
        /// - `instead`: The replacement value to return if the condition matches.
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
