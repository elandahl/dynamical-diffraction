% derivitive of X  to be used with ode45
% for algorithm of Larson and Barnhost, J.A.P. 51 (6), 1980.
% First written by Eric Landahl, 1/30/17. Last revised by EL 1/31/17.
% X = X(1) + i X(2)
% g and k are constants based on structure factors
% A is the depth at which the calculation is requested
% A_y is a vector of depths
% y is a vector of angular variables at different depths y(A_y)
% y and A_y are equally sized vectors
% yy is the value of y(A) interpolated from y(A_y)
% 

function dX = dX(A, X, y, A_y, g, k)

%yy = interp1 (A_y, y, A); % y(A)
yy = y; % overridden for now

dX(1) = k * (X(1)^2 - X(2)^2 + 1) + 2 * X(2) * (X(1) - yy) - 2 * g * X(1);
dX(2) = - (X(1)^2 - X(2)^2 + 1) + 2 * X(1) * (X(2) * k + yy) - 2 * g * X(2);
