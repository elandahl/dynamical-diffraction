% Function to take one depth step of size dz
% in propogating dynamical x-ray diffraction
% in a strained crystal up from the unstrained bulk
% using Wie's algorithim in J. Appl. Phys. 59, 3743 (1986)
% and erratum J. Appl. Phys. 70 (4), 2481 (1991)
% may be used with either fixed or adaptative step sizes
% and vectorized as much as possible for computing speed
% suggestion to vectorize angle thanks to Dohn Arms
% Based on Wie.m by G. Jackson Williams 8/28/2009 and also his later cWie.c
% Written by Eric Landahl, 11/1/2016 
% Last revised by EL 11/7/2016

function [X_out] = Wiestep (X_in, theta, z, dz, strain, param)

%% Output:
%   X_out is the complex scattering amplitude X at z + dz
%     will return X0, the unstrained X, if X_in is an array of zeros
%     scattering intensity is found by Intensity = X_out.*conj(X_out)
%% Inputs
%   X_in is the previous complex scattering intensity X at z (X_in=0 returns X0)
%   z is the depth where X_in was evaluated.
%   dz is the depth step size
%   theta is a vector of angles in radians to evaluate the rocking curve at
%   strain is a strain vector evaluated at z:
%     strain(1) is longitudinal
%     strain(2) is transverse
%     strain(3) is sheer
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
A01 = (z-dz)*(pi*abs(psi_Hr))/(lambda*sqrt(abs(gamma_H*gamma_0))); % Normalizing the crystal thuckness

%% Calculations that depend on depth and strain, but not angle
A01 = (z-dz)*(pi*abs(psi_Hr))/(lambda*sqrt(abs(gamma_H*gamma_0))); % Normalizing the crystal thuckness
A = (pi*z*abs(psi_Hr))/(lambda*sqrt(abs(gamma_H*gamma_0))); % Normalizing iterative crystal depth

%% Calculations that have an angular dependence
%c1 = (cos(phi)^2)*tan(thetaB) + sin(phi)*cos(phi);
%c2 = (sin(phi)^2)*tan(thetaB) - sin(phi)*cos(phi);
c1 = 2*sin(2*thetaB)*((cos(phi)^2)*tan(thetaB) + sin(phi)*cos(phi)); % from Erratum
c2 = 2*sin(2*thetaB)*((sin(phi)^2)*tan(thetaB) - sin(phi)*cos(phi)); % from Erratum
alpha = -2*(theta - thetaB)*sin(2*thetaB) - (c1*strain(1) + c2*strain(2));
y = ((1+b)*psi_0r - b*alpha)/(2*abs(psi_Hr)*sqrt(b));
r1 = abs(y.^2 - g^2 + k^2 - 1);         % APPENDIX FROM WIE et al.
r2 = 2*(y*g - k);                       % Taking the square root of a complex number
r  = sqrt(r1.^2 + r2.^2);               % s = sqrt(C^2 - B^2)
s1 = sqrt((r+r1)/2);
s2 = sqrt((r-r1)/2);
q1 = sqrt(1 - k^2 + g^2);
q2 = k/g;
r1 = abs(y.^2 - g^2 + k^2 - 1);          % APPENDIX FROM WIE et al.
% Using boolean indexing to perform "if/else" statements on vectors
% See: https://www.mathworks.com/matlabcentral/answers/129-how-does-logical-indexing-work
% Start by making arrays of condition tests that return a 1 if condition is true
cond1 = (y <= - q1);
cond2a = (y <= q2);
cond2b = (y >= -q1);
cond2 = cond2a .* cond2b; % logical AND
cond3a = (y <= q1);
cond3b = (y >= q2);
cond3 = cond3a .* cond3b;  % logical AND
cond4 = (y >= q1);
% Finally multiply each condition by the appropriate calculation and add all arrays
s = cond1.*(s1 + 1i*s2) + cond2.*(s2 + 1i*s1) + cond3.*(-s2 + 1i*s1) + cond4.*(-s1 + 1i*s2);
% Continue to calculate X_out for all angles
B = -(1+1i*k);
C = y + 1i*g;
X0 = -B./(C - s);                    % Zero-strain solution
if (max(abs(X_in))==0)
  X_out = X0; % Return zero strain solution if the input X array is zero
else
  X_out = (s.*X_in + 1i*(B+C.*X_in).*tan(s.*(A - A01)))./(s - 1i*(C+B.*X_in).*tan(s.*(A - A01))); 
end

endfunction
