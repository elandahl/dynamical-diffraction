% X0.m
% Function to determine rocking curve complex field X0 at a given angle y
% y is related to angle in radians
% EL 1/31/2017

function X0 = X0(y, g, k)

B = -(1+1i*k);
C = y + 1i*g;
xc = -C./B;
cond01 = (real(xc)>1);
cond02 = (abs(real(xc))<=1);
cond03 = (real(xc)<-1);
X1 = cond01.*(xc-sqrt(xc.^2 -1));
X2 = cond02.*(xc-1i*sqrt(1-xc.^2));
X3 = cond03.*(xc+sqrt(xc.^2 -1));
X0 = X1 + X2  + X3; 

end