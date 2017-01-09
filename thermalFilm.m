% thermalFilm.m
% Classical strain model for a thin Aluminum film on a semiconductor
% Presumes that the temperature rise on the film is uniform and instantaneous
% and that the thermal contact is perfect
% A more sophisticated model with better Physics (e.g. AMM, DMM) will be needed!
% Also, does not include nanoscale phenomena (e.g. phonon mean free path)
% Also, no acoustic propogation is included.  This could be called later.
% Model: see p. 429 (Sec. 10.7) of Hahn, "Heat Conduction" 3rd edition
% Aluminum properties are hard coded
% Semiconductor properties looked up in sample.dat, which is required
% First written by Eric Landahl, 12.28.2016
% Revised by EL 1.9.2017
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

% Load substrate properties
  load sample.dat; 
  ID = find(strcmp({sample.name}, crystal)==1);
  alpha_t = sample(ID).thermalExpansion.val; % 1/K
  alpha2 = sample(ID).thermalDiffusion.val; % cm^2/s
  alpha2 = alpha2 / 10000; % convert from cm^2/s to m^2/s
  k2 = sample(ID).thermalConductivity.val; % W/(cm K)
  k2 = k2 * 100; % Convert from W/(cm K) to W/(m K)
  
% Aluminum film properties
  L = 70e-9; % Film thickness in m
  C_film = 904; %Specific heat of film in J/(kg K)
  rho_film = 2712; % Film density in kg/m^3
  k1 = 204; % Film thermal conductivity in W/(m K)
  alpha1 = 8.418E-5; % Film thermal diffusivity in m^2/s
  
% Calculate initial temperature rise
  fluence = fluence*10; % Convert from mJ/cm^2 to J/m^2
  T0 = fluence/(L * C_film * rho_film); % Initial temperature rise
  fprintf('A %d nm thick Aluminum film gives a temperature rise of %.1f K.\n',L*1e9,T0)
  
% Unitless parameters (see Hahn, "Thermal Conductivity", Eqs. 10-135 and 10-138)  
  mu = sqrt(alpha1/alpha2);
  beta = (k1/k2)/mu;
  gamma = (beta - 1)/(beta + 1);
  
% Spatial grid
  num_depths = 1000;  % number of depth points z to be calculated
  dz = max_depth/num_depths;
  z = dz:dz:max_depth;
  
% Meshgrid for calculation speed & ease
  [Time Z] = meshgrid(time,z); % Time and Z are 2D, time and z are 1D

% Calculate temperature profile
  max_n = 100; % number of terms in series expansion, default 100
  TT = 0.*Time.*Z; % 
  for n = 1: max_n  % Series expansion solution of heat equation
    T1 = erfc((2*n*L + mu*Z)./(2*(sqrt(alpha1*Time)))); % temporary
    T2 = erfc(((2*n + 2)*L + mu*Z)./(2*sqrt(alpha1*Time))); % temporary
    TT = TT + (gamma^n) * (T1 - T2); % temporary
  end
  T = T0 * (1/2) * (1 + gamma) * TT; % Temperature at all z and time  
  fprintf('The maximum bulk temperature rise is %.1f K.\n',max(max(T)));
  
% Calculate strains
  longitudinal = alpha_t*T; % Strain is given by thermal expansion coefficient times temperature
  transverse = 0*T; % No transverse strain
  sheer = 0*T; % No sheer strain
  time_out = time;

  
  