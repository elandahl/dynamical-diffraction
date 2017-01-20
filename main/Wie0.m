% Wie0.m
% Function to determine unstrained crystal rocking curve
% using Wie's algorithim in J. Appl. Phys. 59, 3743 (1986)
% and erratum J. Appl. Phys. 70 (4), 2481 (1991)
% Vectorized as much as possible for computing speed
% suggestion to vectorize angle thanks to Dohn Arms
% Based on Wie.m by G. Jackson Williams 8/28/2009 and also his later cWie.c
% Written by Eric Landahl, 1/20/2017 
% Note sign on cond3 will be different than in Wiestep.m

function [X0] = Wie0 (theta, param)

%% Output:
%   X0 is the complex scattering amplitude of the unstrained crystal
%     scattering intensity is found by Intensity = X0.*conj(X0)
%% Inputs
%   theta is a vector of angles in radians to evaluate the rocking curve at
%   params are the crystal parmeters:
%     param(1) is complex psi0
%     param(2) is complex psiH
%     param(3) thetaB, is the Bragg angle in radians
%     param(4) phi, is the assymetry angle in radians (0 for symmetric)
%     param(5) delta, is the refractive angular offset (close to 0)
%     param(6) lambda, is the wavelength in meters

%% Unwrap parameters
psi_0r = real(param(1));
psi_0i = imag(param(1));
psi_Hr = real(param(2));
psi_Hi = imag(param(2));
thetaB = param(3);
phi = param(4);
delta = param(5);
lambda = param(6);

%% Constants (do not depend on strain or depth)
gamma_0 = sin(thetaB+phi);  % Eventually need to add Chi angle
gamma_H = sin(thetaB-phi);  % Eventually need to add Chi angle
b = abs(gamma_0/gamma_H);
g = ((1+b)*psi_0i)/(2*abs(psi_Hr)*sqrt(b));
k = psi_Hi/psi_Hr;
Lext = lambda*sqrt(abs(gamma_H)*gamma_0)/(pi*abs(psi_Hr));
alpha = -2*(theta - thetaB)*sin(2*thetaB); % contains no strain
y = ((1+b)*psi_0r - b*alpha)/(2*abs(psi_Hr)*sqrt(b));
r1 = abs(y.^2 - g^2 + k^2 - 1);         % APPENDIX FROM WIE et al.
r2 = 2*(y.*g - k);                       % Taking the square root of a complex number
r  = sqrt(r1.^2 + r2.^2);               % s = sqrt(C^2 - B^2)
s1 = sqrt((r+r1)/2); % typo as s2 in WIE et al.  
s2 = sqrt((r-r1)/2);
q1 = sqrt(1 - k^2 + g^2);
q2 = k/g;
% Using boolean indexing to perform "if/else" statements on vectors
% See: https://www.mathworks.com/matlabcentral/answers/129-how-does-logical-indexing-work
% Start by making arrays of condition tests that return a 1 if condition is true
cond1 = (y < - q1);
cond2a = (y <= q2);
cond2b = (y >= -q1);
cond2 = cond2a .* cond2b; % logical AND
cond3a = (y <= q1);
cond3b = (y >= q2);
cond3 = cond3a .* cond3b;  % logical AND
cond4 = (y > q1);
% Finally multiply each condition by the appropriate calculation and add all arrays
s = cond1.*(s1 + 1i*s2) + cond2.*(s2 + 1i*s1) + cond3.*(-s2 + 1i*s1) + cond4.*(-s1 + 1i*s2);
% Note above changed sign in front of cond3 to plus
% Continue to calculate X0 for all angles
B = -(1+1i*k);
C = y + 1i*g;
X0 = -B./(C - s);  
X0=X0/max(X0.*conj(X0));                  % Zero-strain solution

end
