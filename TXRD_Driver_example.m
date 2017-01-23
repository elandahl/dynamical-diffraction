% TRXD_Driver_example.m
% General TRXD master program.  Use this program as a model.
% By Eric Landahl, 1.9.2017
% Updated to add benchmarking 1.16.2017
% Improved accuracy by fixing delta and final step interpolation 1/19/2017

clear all; tic; more off;

%% Include subdirectories in path
addpath('main','include','strain_functions','data','benchmarks');

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
toc;
%% Make plots, smooth data, calculate centroids
if ~strcmp(model,'benchmark') % Unless model is 'benchmark'
  fprintf('TRXD calculation done, preparing plots.\n')
  plot_opts = 'final'; % 'none' = no plots
  ang_res = 2e-4; % angular resultion in degrees FWHM
  [Intensity centroid FWHM] = TRXD_plots (A,A0,times,angles,Strain,z,ang_res,plot_opts);
elseif strcmp(model,'benchmark')
  % Loach benchmark data
  benchmark = load('benchmark.txt');
  Int_BM = benchmark(:,2);
  Ang_BM = benchmark(:,1);
  figure(20);clf;hold all;
    plot(angles*180/pi, A0,'-k')
    plot(angles*180/pi, A,'ob')
    plot(Ang_BM,Int_BM,'r','LineWidth',4)
    xlabel('Theta (degrees)')
    ylabel('Diffracted Intensity')
    legend('Unstrained','TRXD','GID')
    title([num2str(crystal) ' @ ' num2str(energy) ' keV,  ' num2str(reflection) ' and uniform strain of 1E-4, 1 um deep'])
  hold off;
  figure(21);clf;hold all;
    semilogy(angles*180/pi, A0,'-k')
    semilogy(angles*180/pi, A,'ob')
    semilogy(Ang_BM,Int_BM,'r','LineWidth',4)
    xlabel('Theta (degrees)')
    ylabel('Diffracted Intensity')
    legend('Unstrained','TRXD','GID')
    title([num2str(crystal) ' @ ' num2str(energy) ' keV,  ' num2str(reflection) ' and uniform strain of 1E-4, 1 um deep'])
    set(gca,'YScale','log')
  hold off;

end

