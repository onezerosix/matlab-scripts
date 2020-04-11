% tridiagonal matrix solver for x in Ax = d
function [alpha, beta, z, x] = tridiag_lu_decomp(n, a, b, c, d)
% INPUT
% n: size of matrix A
% a: diagonal values of A
% b: lower diagonal values of A (1st value will be ignored)
% c: upper diagonal values of A (nth value will be ignored)
% d: vector in Ax = d
% following criteria must be met for the script to work
% |a1| > |c1|, |an| > |bn|
% |ai| >= |bi| + |ci| and bici != 0 for i = 2,3,...,n-1

% OUTPUT
% alpha: diagonal values of U (upper diagonal values are c)
% beta: lower diagonal values of L (diagonal values are 1)
% z: vector values of z = Ux
% x: solution (vector values of x in Ax = d)

% validate input
assert(n > 2, 'n must be > 2') % tridiagonal
assert(length(a) == n, 'length of a should be n')
assert(length(b) == n, 'length of b should be n')
assert(length(c) == n, 'length of c should be n')
assert(length(d) == n, 'length of d should be n')
assert(abs(a(1)) > abs(c(1)), 'criteria not met: |a1| > |c1|')
assert(abs(a(n)) > abs(c(n)), 'criteria not met: |an| > |cn|')
for i=2:n-1
    assert(b(i)*c(i) ~= 0, 'criteria not met b%d * c%d != 0', i, i)
    assert(abs(a(i)) >= abs(b(i)) + abs(c(i)), ...
        'criteria not met: |a%d| >= |b%d| + |c%d|', i, i, i)
end

% create vectors for output
alpha = zeros(1,n);
beta = zeros(1,n);
z = zeros(1,n);
x = zeros(1,n);

% find the LU decomp of A
alpha(1) = a(1);
for i=2:n
    beta(i) = b(i) / alpha(i-1);
    alpha(i) = a(i) - (beta(i)*c(i-1));
end

% now have LUx = d
% (forward) substitution: z = Ux --> Lz = d
z(1)= d(1); % z1 = d1;
for i=2:n
    z(i) = d(i) - (beta(i)*z(i-1));
end

% (backward) substitution: z = Ux
x(n) = z(n) / alpha(n);
for i=n-1:-1:1
    x(i) = (z(i) - (c(i)*x(i+1)))/alpha(i);
end