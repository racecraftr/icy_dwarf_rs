# icy_dwarf_rs

Rust rewrite of Dr. Marc Neveu's (mars.f.neveu@nasa.gov) IcyDwarf program 
originally written in C. 

## Why?

Rust allows for a lot of things that C doesn't. 
- Cross-platform compatability (a main issue in the original codebase)
- Memory safety
- Better library importing, such as those needed for linear algebra

## Installation & Execution

**Prerequisites**: 
- You need `R` installed on your machine with the `CHNOSZ` package installed globally. 
- The `CHNOSZ` package is written in fortran, so `libgfortran` must be installed. 
- The `IPHREEQC` library installed in order to 
- Given the project *is* written in rust, make sure `cargo` is installed on your machine. 

Note that you 
  may have to set the `$RHOME` variable on your machine differently to what it is configured to. 
  If you are using **MacOS/Linux**, add the following 
  to your `.bashrc`, `.zshrc`, or any shell configuration file you use:
  
```sh
export R_HOME=$(R RHOME)
```

For the `CHNOSZ` package, add this to your shell config file if needed:

```sh
alias gfortran=libgfortran
```

This will let you install the `CHNOSZ` package: 

```r
install.packages("CHNOSZ")
```

From here, installing and running the program is easy. 

**To compile and run:**:

```sh 
git clone https://www.github/racecraftr/icy_dwarf_rs.git
cd icy_dwarf_rs/icy_dwarf
cargo run []
```

**To install to your `$PATH$`**: 
```sh
git clone https://www.github/racecraftr/icy_dwarf_rs.git
cd icy_dwarf_rs/icy_dwarf
cargo install
```

***Hold on, what is `build.rs`?***  
Because Rust uses external C/C++ libraries (IPHREEQC), a build script is used to compile them with the rust code. 

## Running

Running the code is simple
```sh 
icy_dwarf
```

## Libraries used

The libraries used for this project are listed below: 
- `serde`, `serde_repr`, and `toml` for the parsing of TOML 
- `faer` for linear algebra
- `clap` for command-line argument parsing 
- `nalgebra` and `num` for complex values and useful num traits 
- `itertools` for useful functions on iterators
- `extendr` for running `R` code in rust 

Thank you to these libraries for helping me out!
