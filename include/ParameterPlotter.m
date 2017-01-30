
% Load X0h data
% ParameterPlotter.m
% EL 1/26/17.  Plots x-ray crystal parameters
clear all; more off;
crystal = 'GaAs';

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

energy = X0h(:,1);
p_0r = X0h(:,2);
p_0i = X0h(:,3);
p_Hr =  X0h(:,4);
p_Hi = X0h(:,5);
delta = X0h(:,6);
tB_deg = X0h(:,7);

figure(1);hold on;
  plot(energy,p_0r,'-c')
  plot(energy,p_Hr,'-c')
  xlabel('Energy (keV)')
  ylabel('Re(Chi)')
  legend('0','H')
hold off

figure(2);hold on;
  plot(energy,p_0i,'-c')
  plot(energy,p_Hi,'-c')
  xlabel('Energy (keV)')
  ylabel('Im(Chi)')
  legend('0','H')
hold off

figure(3);hold on;
  plot(energy,p_0r,'-c')
  xlabel('Energy (keV)')
  ylabel('delta (rad)')
hold off

figure(4);hold on;
  plot(energy,p_0r,'-c')
  xlabel('Energy (keV)')
  ylabel('Bragg Angle (deg)')
hold off
