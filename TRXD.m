% TRXD.m
% Program to calculate Time-Resolved X-Ray Diffraction 
% By Eric Landahl, DePaul University Physics Department, elandahl@depaul.edu
% First written December 13, 2016
% Based on previous work by Sooheyong Lee (KRISS) and G. Jackson Williams (LLNL)
%
% INPUTS:
%   model       strain model, selected from the following list:
%     thermal   Gaussian solution to prompt surface temperature rise
%   crystal     determines x-ray and strain properties, chosen from:
%     GaAs
%     Si
%     Ge
%     InSb
%   reflection  diffraction vector [h k l], i.e. [hkl] chosen from:
%     [0 0 4]
%   cut         surface normal vector [h k l], i.e. (hkl) chosen from:
%     [0 0 1]   
%   energy      x-ray energy in keV  
%   fluence     absorbed laser fluence in mJ/cm^2
%   angle       in degrees, one of the following:
%     a vector  angles to be calculated, relative to the Bragg angle
%     a value   a total range of angles, from which a vector is generated
%     0         a default is used
%   time        in seconds, one of the following:
%     a vector  times to be calculated, relative to time-zero
%     a value   a final time, from which a vector of times is generated
%     0         a default is used, based on the model
%
% OUTPUTS:
%   A           X-ray scattering amplitude array, returned as A(time, angle)             
%   time        a vector of times
%   angle       a vector of angles calculated, absolute, in degrees
%
% SAMPLE Usage:
% [A time angle] = TRXD ('thermal', 'Si', [0 0 4], [0 0 1], 10, 1, 0, 0)


function [A time angle] = TRXD (model, crystal, reflection, cut, energy, fluence, angle, time)

%% Some constants that can be changed to speed things up when using defaults
num_times = 50;
num_angles = 50;


%% Check inputs

if nargin ~= 8
  fprintf('Incorrect number of input arguments\n')
  nargin
  return
end

if (energy < 7) | (energy > 14)
  fprintf('Energy out of range\n')
  energy
  return
end

if reflection ~= [0 0 4]
  fprintf('Only [0 0 4] reflection is supported now.\n')
  reflection
end

if cut ~= [0 0 1]
  fprintf('Only [0 0 1] surface cut is supported now.\n')
  cut
  return
end

% Determine diffraction parameters

% Calculate assymetry angle; should be zero
phi = acos(dot(reflection/norm(reflection),cut/norm(cut)));
if phi > 1e-6
  fprintf ('The assymetry angle is not zero.\n')
  phi
  return
end

% Load X0h data
if strcmp(crystal,'GaAs')
  X0h = load ('GaAs400.dat');
  ID = 1; % ID number for crystal data
elseif strcmp(crystal,'InSb')
  X0h = load ('InSb400.dat');
  ID = 2;
elseif strcmp(crystal,'Si')
  X0h = load ('Si400.dat');
  ID = 3;
elseif strcmp(crystal,'Ge')
  X0h = load ('Ge400.dat');
  ID = 4;
else
  fprintf('Crystal needs to be either Si, GaAs, Ge, or InSb.\n')
  crystal
  return
end  

% Interpolate X0h data
% Note that X0h(:,1) are the energies in keV
p_0r = -abs(interp1(X0h(:,1), X0h(:,2), energy, 'spline', 'extrap'));   
p_0i = -abs(interp1(X0h(:,1), X0h(:,3), energy, 'spline', 'extrap'));   
p_Hr = -abs(interp1(X0h(:,1), X0h(:,4), energy, 'spline', 'extrap'));   
p_Hi = -abs(interp1(X0h(:,1), X0h(:,5), energy, 'spline', 'extrap'));   
delta = interp1(X0h(:,1), X0h(:,6), energy, 'spline', 'extrap');   
tB_deg = interp1(X0h(:,1), X0h(:,7), energy, 'spline', 'extrap'); 
thetaB= tB_deg*3.14159/180;  % Sergey Si, Bragg angle in radians
lambda= energy*1.23984E-11; % convert energy in keV to wavelength in meters

% Assemble params arrray containing parameters for the crystal diffraction
params(1) = p_0r + 1i*p_0i;
params(2) = p_Hr + 1i*p_Hi;
params(3) = thetaB;
params(4) = phi;
params(5) = delta;
params(6) = lambda; 

% Assemble opts array containing options for adaptative depth stepping
tol = 1e-6; % tolerance.  Higher value for more speed and less precision
dz_min = 1.1e-10; % Minimum step size in meters
dz_max = 1e-7; % Maximum step size in meters
f = 2; %Shift factor for convergence

opts(1) = tol;
opts(2) = dz_min;
opts(3) = dz_max;
opts(4) = f;

% Calculate extinction length
gamma0 = sin(thetaB+phi); 
gammaH = sin(thetaB-phi);
Lext = lambda*sqrt(abs(gammaH)*gamma0)/(pi*abs(p_Hr));

%% Calculate time array 
if time == 0
  time = 1e-8; % set default maximum time delay to 10 ns
end
if (length(time)==1) 
  time = time/num_times:time/num_times:time;
end

%% Calculate angular array
if (angle == 0)
  angle = 1e-2; % set default angular width to 10 mdeg
end
if (length(angle)==1)
  angle = (0:angle/num_angles:angle)-angle/2;
end

if strcmp(model,'thermal')
  % Load properties
  semi_data=csvread('Semiconductor_properties.csv');
  alpha_T = semi_data(7,ID+3); % 1/K
  kappa_T = semi_data(9,ID+3); % cm^2/s
  kappa_T = kappa_T / 10000; % convert from cm^2/s to m^2/s
  % Calculate initial temperature rise
  fluence = fluence*10; % Convert from mJ/cm^2 to J/m^2
  t_film = 70e-9; % Film thickness in m
  C_film = 904; %Specific heat of film in J/(kg K)
  rho_film = 2712; % Film density in kg/m^3
  T0 = fluence/(t_film * C_film * rho_film); % Initial temperature rise
  fprintf('A %d nm thick Aluminum film gives a temperature rise of %.1f K\n',t_film*1e9,T0)
  % Calculate teperature profile
  sigma = sqrt(2*kappa_T*time); % Width of temperature pulse
  dz = sigma(1); % Let depth step be the width at the first timepoint (FTCS Stability critereon)
  z = dz:dz:5*Lext;
  for m = 1:num_times % time loop
    for n = 1:length(z) % depth loop
      T(m,n) = dz*(T0/sigma(m))*(1/sqrt(2*pi))*exp(-(z(n)^2)/(2*sigma(m)^2));
    end
  end
  strain = alpha_T*T; % Strain is given by thermal expansion coefficient times temperature
else
  fprintf('Not a valid model\n')
  return
end

%%%%% 12/15/16: Next make "thermal" a separate function.  Then call WieAdapt.

%% Eventually, a strain propogation algorithim (with remeshing) goes here

 %%%for testing
A =strain;


end
