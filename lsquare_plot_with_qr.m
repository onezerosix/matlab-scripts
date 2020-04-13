function [x] = lsquare_plot_with_qr(t, y, basis, plot_step)
% LSQUARE_PLOT_WITH_QR Solves Ax = y for x in the least-square sense and
% plots the least-square fit using QR decomposition.
%
% INPUT
% t & y: m paired input cordinates
% basis: cell array of n basis function handles that have t as input
%          these are used for forming A & plotting the least-square fit
% plot_step: (OPTIONAL) define the step for least-square fit plot,
%               defaults to 0.01
% OUTPUT
% a plot of t & y and a least-square fit function for the data
% x: least-square fit solution of Ax = y

% validate and possibly default input
m = length(t);
assert(m > 0, 'length of t must be > 0')
assert(length(y) == m, 'length of y must be same as t')
n = length(basis);
assert(n > 0, 'length of basis must be > 0')
if nargin == 3 || plot_step <= 0
    plot_step = 0.01;
end

% form matrix A
A = zeros(m,n);
for i = 1:m
    for j = 1:n
        A(i,j) = basis{j}(t(i));
    end
end

% find QR decomposition and w & U
[Q,R] = qr(A);
w = Q'*y;
w = w(1:n, :);
U = R(1:n, :);

% backward subs to solve x of Ux = w because U is upper triangular
x = zeros(n,1);
for i = n : -1 : 1
    other_terms = 0;
    for j = n : -1 : i+1
        other_terms = U(i,j)*x(j) + other_terms;
    end
    x(i) = (w(i) - other_terms) / U(i,i);
end

% make the curve from basis functions
lsquare_x = min(t) : plot_step:max(t);
lsquare_y = x(1) * basis{1}(lsquare_x);
for i = 2 : n
    lsquare_y = lsquare_y + x(i) * basis{i}(lsquare_x);
end

% plot the given data points and the curve
scatter(t,y);
hold on
plot(lsquare_x,lsquare_y);
hold off

end