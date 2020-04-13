# misc-matlab-scripts

A collection of MATLAB scripts.

The input values for each function are generally not validated other than their lengths. If an error such as division by 0 occurs, it will propagate. Please refer to each function's help documentation for additional information.

### Table of Contents

1. [Tridiagonal LU Decomposition Matrix Solver](#Tridiagonal-LU-Decomposition-Matrix-Solver)  
  a. [Example](#Example~1)
2. [Least-Square Plot Using QR Decomposition](#Least-Square-Plot-Using-QR-Decomposition)  
  b. [Example](#Example~2)

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
### Least-Square Plot Using QR Decomposition

This script solves `Ax = y` for the n x 1 vector `x` in the least-square sense using QR decomposition. It also plots the least-square fit function along with the input coordinates for comparison.

The input n x 1 `basis` functions (with m x 1 input `t`) are used to form the m x n matrix `A`. For example, for n = 2 basis functions `b_1(t) = e^-t` and `b_2(t) = e^-2t` with m = 3 input coordinates given by `t = [0; 1; 2]` and `y = [0; 2; 1.06]`, we get the following 3 x 2 matrix:

```
    [ e^0   e^0  ]
A = [ e^-1  e^-2 ]
    [ e^-2  e^-4 ]
```

Two assumptions made for the script are: `A` has rank n and m >= n.

##### Example

```
>> [x] = lsquare_plot_with_qr([0;.5;1;1.3;2;3], [0;1.6;2;1.93;1.06;0.38], {@(t)exp(-t), @(t)exp(-2*t)})

x =
    8.3282
   -8.4245
```
<img src="https://raw.githubusercontent.com/onezerosix/misc-matlab-scripts/master/pictures/lsquare_plot_with_qr_example_graph.png" width="45%" alt="example plot">
