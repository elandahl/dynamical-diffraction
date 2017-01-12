% TRXD_Driver_File.m
% General TRXD master program. 
% To read in a 2D strain array from a previously genreated file.
% Prepared at Soo's request and tested with Soo1.dat
% By Eric Landahl, 1.10.2017
%
clear all; tic; more off;



%% Genreate fresh sample material properties data file
%sampledata; % creates file sample.dat database of material properties

%% Calculate TRXD for Si
model = 'strainFile';
crystal = 'Si';
reflection = [0 0 4];
cut = [0 0 1];
energy = 10; % in keV
filename = 'Soo1.txt'; % Note that this spot is usually for "fluence", but is filename for external files
angles = 0; % deg. relative, use 0 for default angles
times = 10e-12; % for external files, this is the timestep
fprintf('Starting TRXD calculation.\n')
[A A0 times angles Strain z] = TRXD (model, crystal, reflection, cut, energy, filename, angles, times);

fprintf('TRXD calculation done, preparing plots.\n')
%% Make plots, smooth data, calculate centroids
plot_opts = 'final'; % 'none' = no plots
ang_res = 2e-4; % angular resultion in degrees FWHM
[Intensity centroid FWHM] = TRXD_plots (A,A0,times,angles,Strain,z,ang_res,plot_opts);

toc