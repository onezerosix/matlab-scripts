function [x] = inverse_quad_interp(x, func)
% INVERSE_QUAD_INTERP Gives the result of interpolating the given function
% by the inverse quardatic method once and can continue interpolating until
% the user says to stop.
%
% INPUT
% x: the first three x values to use for interpolating
% func: the function to use
%
% OUTPUT
% x: the final three values of x (after user says to stop interpolating)

keep_looping = 1;

while keep_looping ~= 0
    y = [func(x(1)) func(x(2)) func(x(3))];

    % calculate gamma, alpha, beta
    gamma = (x(1) - x(3)) / (y(1) - y(3));
    alpha = (x(1) - x(2)) / (y(1) - y(2));
    beta = (alpha - gamma) / (y(2) - y(3));

    % compute next estimates
    nextX = x(1) - alpha*y(1) + beta*y(1)*y(2);
    x = [nextX, x(1), x(2)];
    
    disp(x);
    keep_looping = input('To stop interpolating, type 0, else 1: ');
end

end