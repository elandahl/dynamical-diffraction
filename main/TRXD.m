% TRXD.m
% Program to calculate Time-Resolved X-Ray Diffraction 
% By Eric Landahl, DePaul University Physics Department, elandahl@depaul.edu
% First written December 13, 2016
% Last revised by EL 1.9.2017 to add unstrained amplitude output
% Based on previous work by Sooheyong Lee (KRISS) and G. Jackson Williams (LLNL)
%
% INPUTS:
%   model           strain model, selected from the following list:
%     thermalFilm   Gaussian solution to prompt surface temperature rise
%   crystal         determines x-ray and strain properties, chosen from:
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
%   A0          unstrained crystal scattering amplitude           
%   time        a vector of times
%   angle       a vector of angles calculated, absolute, in degrees
%
% SAMPLE Usage:
% [A A0 time angle Strain z] = TRXD ('thermalFilm', 'Si', [0 0 4], [0 0 1], 10, 1, 0, 0)


function [A A0 time angle Strain_save z] = TRXD (model, crystal, reflection, cut, energy, fluence, angle, time)
more off; % Turn on scrolling

%% Include subdirectories in path
addpath('main','include','strain_functions','data');

%% Some constants that can be changed to speed things up when using defaults
num_times = 20;
num_angles = 100;
time_f = 5e-7; % in seconds
angle_width = 1e-2; % in degrees



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
tol = 1e-1; % tolerance.  Higher value for more speed and less precision
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
  time = time_f; % set default maximum time delay to default
end
if (length(time)==1) 
  dt =time; % save for external file function
  time = time/num_times:time/num_times:time;
end

%% Calculate angular array
if (angle == 0)
  angle = angle_width; % set default angular width to 10 mdeg
end
if (length(angle)==1)
  angle = (0:angle/num_angles:angle)-0.9*angle/2;
end

angle = thetaB + angle*pi/180; % Convert angle to radians and add Bragg Angle

if strcmp(model,'thermalFilm')
  [longitudinal transverse sheer time_out z] = thermalFilm (crystal, fluence, time, 5.1*Lext);
elseif strcmp(model,'strainFile')
  filename = fluence;
  dz = 10e-9; % depth step in meters
  [longitudinal transverse sheer time_out z] = strainFile (filename, dz, dt);
else
  fprintf('Not a valid model\n')
  return
end
 
 % Need to add a single timepoint test for comparisson with Sergey's GID
 % Otherwise seems to work 12/28/16 EL
 
for m = 1: length(time)
  if time(m) == time_out(m) % if the time doesn't need to be remeshed
   st1(m,:) = longitudinal(m,:);
   st2(m,:) = transverse(m,:);
   st3(m,:) = sheer(m,:);
  else % remesh time if necessary
   st1(m,:) = interp1(time_out,longitudinal,time(m), 'spline', 'extrap');
   st2(m,:) = interp1(time_out,transverse,time(m), 'spline', 'extrap');
   st3(m,:) = interp1(time_out,sheer,time(m), 'spline', 'extrap');
  end
  Strain_in = [st1(m,:)' st2(m,:)' st3(m,:)']; % Evaluate strain at each time
  [X X0 err Steps Strain_out] = WieAdapt (Strain_in, z, angle, opts, params);
  A0 = X0.*conj(X0); % Unstrained crystal intensity
  A(m,:) = X.*conj(X); % Intensity Amplitude
  Strain_save(m,:,:) = Strain_out(:,:);
end

angle = (angle-thetaB-2*delta)*180/pi; % Convert angle back to degrees relative to Bragg

