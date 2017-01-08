% thermalFilm.m
% Simple strain model for a thin Aluminum film on a semiconductor
% Presumes that the temperature rise on the semiconductor surface
% is the same as the temperature jump of the film
% A more sophisticated model with better Physics (e.g. AMM, DMM) will be needed!
% First written by Eric Landahl, 12.28.2017
% Usually called by TRXD.m
%
%% INPUTS:
%   crystal     determines x-ray and strain properties, chosen from:
%     GaAs
%     Si
%     Ge
%     InSb
%   fluence     absorbed laser fluence in mJ/cm^2
%   time        a vector of times to be calculated in seconds
%   max_depth   usually 5*Lext, Lext is the x-ray extinction length in meters
%
%% OUTPUTS:
%   longitudinal    longitudinal strain, size = length(time_out) x length(z)
%   transverse      transverse strain, size = length(time_out) x length(z)
%   sheer           sheer strain, size = length(time_out) x length(z)
%   z               a vector of depths in meters
%   time_out        a vector of times returned in seconds, not the same as time_in
%
%% NOTE: The strain will be calculated out to time_in(end) however the time-steps
%        are chosen for convenience and accuracy in evaluating the strain. 
%        time_in is NOT usually equal to time_out
%        therefore the strain must be interpolated temporally after calling this
%        function
%
%% NOTE ALSO: Three strain components are returned, but usually only longitudinal
%             is non-zeros
%
%% TYPICAL USAGE
% 
% [st1 st2 st3 time_out z] = thermalFilm ('Si', 1, (1e-10:2e-10:1e-8), 1e-5);
%
function [longitudinal transverse sheer time_out z] = thermalFilm (crystal, fluence, time_in, max_depth)

% Remesh time
  time = time_in; 
  % In this simple model, there is no penalty for calculating many timepoints
  % so the given array of timepoints is used. 
  % Depth steps are chose below.

% Load properties
  load sample.dat; 
  ID = find(strcmp({sample.name}, crystal)==1);
  alpha_t = sample(ID).thermalExpansion.val; % 1/K
  D_t = sample(ID).thermalDiffusion.val; % cm^2/s
  D_t = D_t / 10000; % convert from cm^2/s to m^2/s
  
% Calculate initial temperature rise
  fluence = fluence*10; % Convert from mJ/cm^2 to J/m^2
  t_film = 70e-9; % Film thickness in m
  C_film = 904; %Specific heat of film in J/(kg K)
  rho_film = 2712; % Film density in kg/m^3
  T0 = fluence/(t_film * C_film * rho_film); % Initial temperature rise
  fprintf('A %d nm thick Aluminum film gives a temperature rise of %.1f K\n',t_film*1e9,T0)
  
% Calculate teperature profile
  sigma = sqrt(2*D_t*time); % Width of temperature pulse
  dz = min(sigma(1),max_depth/100); % Let depth step be the width at the first timepoint (FTCS Stability critereon)
  z = dz:dz:max_depth;
  for m = 1:length(time) % time loop
    for n = 1:length(z) % depth loop
      T(m,n) = dz*(T0/sigma(m))*(1/sqrt(2*pi))*exp(-(z(n)^2)/(2*sigma(m)^2));
    end
  end
  
% Calculate strains
  longitudinal = alpha_t*T; % Strain is given by thermal expansion coefficient times temperature
  transverse = 0*T; % No transverse strain
  sheer = 0*T; % No sheer strain
  time_out = time;

  
  