% StrainCheck.m
% Program to check the arguments passed to WieAdapt.m
% Return warnings
% And fix any issues with the options (opts) array 
% And taper the strain function to zero at 5*Lext
% Written by Eric Landahl November 7, 2016
% Last revised by EL 12/7/16 to make compatible with standard MATLAB syntax

function [warnings Strain opts] = StrainCheck (Strain, z, theta, opts, params)

warnings = '';

%% Check number of arguments and length of arrays

if (nargin ~= 5)
  warnings = strcat(warnings, 'Incorrect number of input arguments to AdaptWie.m  ');
end

[M N] = size(Strain);
if (M ~= length(z))
  warnings = strcat(warnings, 'Strain and length not the same size.  ');
end

if (N ~= 3)
  warnings = strcat(warnings, 'Strain needs to have 3 elements at each depth: [long trans sheer].  ');
end

if (length(opts) ~= 4 )
  warnings = strcat(warnings, 'There need to be four elements in [opt].  Use 0 for defaults.  ');
end

if (length(params) ~= 6 )
  warnings = strcat(warnings, 'There need to be six elements in [params].  Use 0 for defaults.  ');
end

%% Unpack params
psi0 = params(1);
psiH = params(2);
thetaB = params(3);
phi = params(4);
delta = params(5);
lambda = params(6);

%% Unpack opts
tol = opts(1);
dz_min = opts(2);
dz_max = opts(3);
f = opts(4);

%% Check values of params

%if (~iscomplex(psi0) | ~iscomplex(psiH))
%  warnings = strcat(warnings, 'Psi are not complex.  ');
%end

if (thetaB <= 0 | thetaB >= pi/2) 
  warnings = strcat(warnings, 'The Bragg angle is out of range, or not in radians.  ');
end

if (abs(delta) >= pi/180)
  warnings = strcat(warnings, 'Refractive shift angle delta is out of range or not in radians.  ');
end

if (lambda <= 1e-12 | lambda >= 1e-8) 
  warnings = strcat(warnings, 'Lambda is out of range or not in meters.  ');
end

%% Check strain and depth matrices

if (max(z) >= 1e-3 )
  warnings = strcat(warnings, 'The maximum depth is out of range, or z is not in meters.  ');
end

if (length(z) <= 100)
  warnigns = strcat(warnings, 'The depth matrix z and Strain matrix should have at least 100 points.  ');
end

%% Condition to use fixed step size. 
if (tol <= 0)
 tol = 1e-10; % Set tolerance to 1 Angstrom
 f = 1; % and use fixed step size
 fprintf('Using fixed step size./n');
end

%% Using default minimum and maximum step sizes
if (dz_min <= 1e-10)
  dz_min = 1e-10; 
  fprintf('Depth steps have been set to a default minimum of 1 Angstrom.  ');
end
if (dz_max <= 1e-10 | dz_max >= 1.1e-6 )
  dz_max = 1e-6; 
  fprintf('Depth steps have been set to a default maximum size of 10,000 Angstrom.  ');
endif

%% Check if shift factor is reasonable
if (f <= 0)
  f = 1; 
  fprintf('Using fixed step size./n')
end

%% Check angle range
if (min(theta) <= 0 | max(theta) >= pi) 
  warnings = strcat(warnings, 'The angle range is too large, or not in radians.  ');
end
if (min(theta) >= thetaB | max(theta) <= thetaB)
  warnigns = strcat(warnings, 'The Bragg angle does not fall within the range of angles.  ');
end

%% Check Strain magnitude
if (max(max(abs(Strain))) >= 1e-2)
  warnings = strcat(warnings, 'The maximum strain value is greater than 1%.  ');
end

%% Check extinction depth
gamma0 = sin(thetaB+phi); 
gammaH = sin(thetaB-phi);
Lext = lambda*sqrt(abs(gammaH)*gamma0)/(pi*abs(real(psiH)));
if (5*Lext >= z(end))
  warnings = strcat(warnings, 'The depth array z and Strain do not go deep enough, to 5 x Lext.  ');
end

%% Strain tapering:  keeps strain approximately constant from 0 to 3*Lext, than tapers to ~0 at 5*Lext
taper = 0.5*(1-erf(2*z/Lext - 8));
Strain = [Strain(:,1).*taper' Strain(:,2).*taper' Strain(:,3).*taper'];

%% Update and pack up opts
opts(1) = tol;
opts(2) = dz_min;
opts(3) = dz_max;
opts(4) = f;

end
