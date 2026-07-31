//! This build script configures Cargo to link the `iphreeqc` C++ library into the Rust binary executable.

/// Configure cargo native library search paths and dynamic linker flags.
fn main() {
    println!("cargo:rustc-link-search=native=/usr/local/lib");
    println!("cargo:rustc-link-lib=iphreeqc");
    println!("cargo:rustc-link-arg=-Wl,-rpath,/usr/local/lib");
}

