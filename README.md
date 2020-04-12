# misc-matlab-scripts

A collection of MATLAB scripts.

### Table of Contents

1. [Tridiagonal LU Decomposition Matrix Solver](#Tridiagonal-LU-Decomposition-Matrix-Solver)  
  a. [Example](#Example)

### Tridiagonal LU Decomposition Matrix Solver

This script solves `Ax = d` for `x` using LU decomposition for the square tridiagonal matrix `A`.

A square tridiagonal matrix `A` has the following structure:  
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/tridiag_lu_decomp_A.png" width="40%" alt="matrix A">

The `a`, `b` and `c` variables are used as vector inputs for the script. The values must meet the following criteria for the script to work:  
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/tridiag_lu_decomp_A_conditions.png" width="50%" alt="matrix A conditions">

The output includes vector `x` and vector `z` such that `z = Ux`. It also includes vectors `alpha` and `beta` such that:  
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/tridiag_lu_decomp_LU.png" width="70%" alt="L and U matrices">

##### Example
```
>> [alpha, beta, z, x] = tridiag_lu_decomp([2;2;2;2;2], [1;1;1;1;1], [1;1;1;1;1], [3;4;4;4;3])

alpha =
    2.0000    1.5000    1.3333    1.2500    1.2000

beta =
         0    0.5000    0.6667    0.7500    0.8000
         
z =
    3.0000    2.5000    2.3333    2.2500    1.2000

x =
     1     1     1     1     1
```
