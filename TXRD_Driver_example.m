% TRXD_Driver_example.m
% General TRXD master program.  Use this program as a model.
% By Eric Landahl, 1.9.2017
%
clear all; tic; more off;

%% Include subdirectories in path
addpath('main','include','strain_functions','data');

%% Genreate fresh sample material properties data file
sampledata; % creates file sample.dat database of material properties
%% Calculate TRXD for Si
model = 'thermalFilm';
crystal = 'GaAs';
reflection = [0 0 4];
cut = [0 0 1];
energy = 10; % in keV
fluence = 1; % in mJ/cm^2
angles = 0; % deg. relative, use 0 for default angles
times = logspace(-3,3,30)*1e-9; % in seconds; use 0 for default times
fprintf('Starting TRXD calculation.\n')
[A A0 times angles Strain z] = TRXD (model, crystal, reflection, cut, energy, fluence, angles, times);

fprintf('TRXD calculation done, preparing plots.\n')
%% Make plots, smooth data, calculate centroids
plot_opts = 'final'; % 'none' = no plots
ang_res = 2e-4; % angular resultion in degrees FWHM
[Intensity centroid FWHM] = TRXD_plots (A,A0,times,angles,Strain,z,ang_res,plot_opts);

toc
