function [m, a, b] = nat_clamped_cubic_spline(knot_x, knot_y, plot_step)
% NAT_CLAMPED_CUBIC_SPLINE Finds a natural clamped cubic spline for the
% input knots then plots the spline. The spline will have the form of
% S_i(x) = -(x - x_i+1)^3 * m_i / 6h_i + (x - x_i)^3 * m_i+1 +
%           a_i * (x_i+1 - x) + b_i * (x - x_i)
% where a_i = y_i / h_i - h_i * m_i / 6 for i = 1,2,...,n-1
% and b_i = y_i+1 / h_i - h_i * m_i+1 / 6 for i = 1,2,...,n-1
% see below for h and m.
%
% A tridiagonal matrix A with the equation Am = z will be used to solve for
% the m vector. The system of equations is
% h_i-1 * m_i-1 + 2(h_i-1 + h_i) * m_i + h_i * m_i+1 = 6(d_i - d_i-1)
% for i = 2, ..., n-1
% where h_i = x_i+1 - x_i and d_i = (y_i+1 - y_i) / h_i.
% 
% INPUT
% knot_x, knot_y: knots for the cubic spline with x in ascending order.
%                 Used for x_i, x_i+1, y_i and y_i+1 for the above
%                 equations.
% plot_step: (OPTIONAL) define the step for least-square fit plot,
%               defaults to 0.05.
%
% OUTPUT
% m: size n vector from Am = z (see above)
% a & b: size n vectors from S_i (see above)

% validate and possibly default input
n = length(knot_x);
assert(n > 0, 'length of knot_x must be > 0')
assert(length(knot_y) == n, 'length of knot_y must be same as knot_x')
if nargin == 2 || plot_step <= 0
    plot_step = 0.01;
end

% prepare variables
h = zeros(1, n-1); % knot x diffs
d = zeros(1, n-1); % knot slopes

for i = 1:(n-1)
	h(i) = knot_x(i+1) - knot_x(i);
	d(i) = (knot_y(i+1) - knot_y(i)) / h(i);
end

% create tri diag A and z of Am = z
diagonals = zeros(1, n-1);
lower_diagonals = zeros(1, n-1);
upper_diagonals = zeros(1, n-1);
z = zeros(1, n-1);

for i = 2:(n-1) % natural spline so i = 1 & n will need to have m = 0
	diagonals(i) = 2 * (h(i-1) + h(i));
	lower_diagonals(i) = h(i-1);
	upper_diagonals(i) = h(i);
	z(i) = 6 * (d(i) - d(i-1));
end

% first m should be 0 (already have necessary vars set to 0)
diagonals(1) = 1;

% solve Am = z for m
[~, ~, ~, m] = tridiag_lu_decomp(diagonals, lower_diagonals, ...
    upper_diagonals, z);

m = [m 0]; % set nth m to 0

% prepare more variables
a = zeros(1, n-1);
b = zeros(1, n-1);
 
% prepare for plot
num_points = 0;
for i = 1 : (n-1)
    num_points = num_points + 1 + (knot_x(i+1) - knot_x(i))/plot_step;
end
x_points = zeros(1, num_points);
y_points = zeros(1, num_points);

% solve a, b and get plot points
x_points_index = 1;
y_points_index = 1;
for i = 1 : (n-1)
    x = knot_x(i) : plot_step : knot_x(i+1);
    
    a(i) = (knot_y(i) / h(i)) - (h(i) * m(i) / 6);
    b(i) = (knot_y(i+1) / h(i)) - (h(i) * m(i+1)/6);
    
    for j = 1 : length(x)
        y_points(y_points_index) = ...
            -((x(j) - knot_x(i+1))^3 * m(i) / (6*h(i))) + ...
            ((x(j) - knot_x(i))^3 * m(i+1) / (6*h(i))) + ...
            (a(i) * (knot_x(i+1) - x(j))) + (b(i) * (x(j) - knot_x(i)));
        y_points_index = y_points_index + 1;
    end
        
    x_points(x_points_index : x_points_index + length(x) - 1) = x;
    x_points_index = x_points_index + length(x);
end

plot(x_points, y_points, knot_x, knot_y, 'o');
end