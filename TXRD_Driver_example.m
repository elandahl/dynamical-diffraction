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
model = 'strainFile1D';
crystal = 'GaAs';
reflection = [0 0 4];
cut = [0 0 1];
energy = 10; % in keV
fluence = 500; % in mJ/cm^2, or the timepoint for strainFile1D
angles = 0; % deg. relative, use 0 for default angles
times = logspace(-3,3,30)*1e-9; % in seconds; use 0 for default times, overridden for external files
fprintf('Starting TRXD calculation.\n')
[A A0 times angles Strain z] = TRXD (model, crystal, reflection, cut, energy, fluence, angles, times);
toc;
if strcmp(model,'benchmark')
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
    title([num2str(crystal) ' @ ' num2str(energy) ' keV,  ' num2str(reflection) ' and uniform strain of 1E-4, 2 um deep'])
  hold off;
  figure(21);clf;hold all;
    semilogy(angles*180/pi, A0,'-k')
    semilogy(angles*180/pi, A,'ob')
    semilogy(Ang_BM,Int_BM,'r','LineWidth',4)
    xlabel('Theta (degrees)')
    ylabel('Diffracted Intensity')
    legend('Unstrained','TRXD','GID')
    title([num2str(crystal) ' @ ' num2str(energy) ' keV,  ' num2str(reflection) ' and uniform strain of 1E-4, 2 um deep'])
    set(gca,'YScale','log')
  hold off;
elseif strcmp(model,'strainFile1D')
  figure(30);clf; hold on;
    plot(angles, A, '-b')
    xlabel('Angle (deg)')
    ylabel('Intensity')
    plot(angles, A0, '-k')
    legend('Strained','Unstrained')
    xlim([angles(floor(end/4)) angles(floor(3*end/4))])
    title([num2str(times*1e9) ' ns'])
  hold off;
  figure(31);clf; hold on;
    semilogy(angles, A, '-b')
    xlabel('Angle (deg)')
    ylabel('Intensity')
    semilogy(angles, A0, '-k')
    legend('Strained','Unstrained')
    title([num2str(times*1e9) ' ns'])
  hold off;
  figure(32);clf;hold on;
    plot(z*1e6,Strain(1,:,1));
    xlabel('Depth (um)')
    ylabel('Strain')
    title([num2str(times*1e9) ' ns'])
  hold off;
  dtheta =  sum(A.*angles)/sum(A) - sum(A0.*angles)/sum(A0); % in degrees
  strain_est = -(pi/180)*dtheta*cot(max(angles)*pi/180);
  fprintf('The centroid shift is %.3f mdeg or %.1e radians.\n',dtheta*1000,dtheta*pi/180)
  fprintf('For %s %d %d %d at %.1f keV the estimated average strain is %.1e .\n',crystal,reflection,energy,strain_est)
%% Make plots, smooth data, calculate centroids
else  % Unless model is 'benchmark'
  fprintf('TRXD calculation done, preparing plots.\n')
  plot_opts = 'final'; % 'none' = no plots
  ang_res = 2e-4; % angular resultion in degrees FWHM
  [Intensity centroid FWHM] = TRXD_plots (A,A0,times,angles,Strain,z,ang_res,plot_opts);

end

