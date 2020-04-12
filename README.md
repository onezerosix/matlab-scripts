# misc-matlab-scripts
A collection of MATLAB scripts.

### Table of Contents
1. [Tridiagonal LU Decomposition Matrix Solver](#Tridiagonal-LU-Decomposition-Matrix-Solver)

### Tridiagonal LU Decomposition Matrix Solver
This script solves `Ax = d` for `x` using LU decomposition for the square tridiagonal matrix `A`.

A square tridiagonal matrix `A` has the following structure:\
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/tridiag_lu_decomp_A.png" width="40%" alt="matrix A">

The `a`, `b` and `c` variables are used as vector inputs for the script. The values must meet the following criteria for the script to work:\
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/tridiag_lu_decomp_A_conditions.png" width="50%" alt="matrix A conditions">

The output includes vector `x` and vector `z` such that `z = Ux`. It also includes vectors `alpha` and `beta` such that:\
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/tridiag_lu_decomp_LU.png" width="70%" alt="L and U matrices">
