% LBM.m
% Simple 1D Lattice Boltzmann Model for Temperature
% Based on Escobar et al., Intl. J. of Heat and Mass Transfer 49, 97-107 (2006)
% First written by Eric Landahl, 12-21-2016
%
% Requries the functions:
%   Debye.m     finds the phonon energy density given a temperature and debye temperature
%   DebyeInv.m  find the temperature given a phonon energy density and debye temperature
%

clear all;
%
tau = 20; % In units of delta-T
t_max = 50; % Number of time steps
x_max = 500; % number of depth steps
%T0 = 100; % temperature jump
%T_room = 300; % equilibrium temperature
%Theta = 640; % Debye Temperature in K, Si is 640 K
W = 1/tau; % grid weight factor = delta-T/tau
b = 0;
%
%
%% use convention T(x)
%
%T = zeros(x_max) + T_room;
%T(floor(x_max/2)) = T0 + T_room;
%
%% Calculate Energy at all positions, in all directions
%E = Debye(T,Theta); % Eq. 6

E = zeros(x_max);
x0=floor(x_max/2);
E(x0-b:x0+b) = 1; % Heat spike

%plot(E)
%ylim([0 1])
%jnk = input('press enter to continue');
EL = E/2; 
ER = E/2;

for t = 1:t_max

% Divide the energy into left and right phonons


EL_new(1:x_max-1) = (1-W)*EL(2:x_max) + W*E(1:x_max-1)/2;
ER_new(2:x_max) = (1-W)*ER(1:x_max-1) + W*E(2:x_max)/2;

EL_new(x_max) = EL(x_max);
ER_new(1) = ER(1);

E = EL_new + ER_new;
%E(x0-b:x0+b) = E(x0-b:x0+b) + (t<t_max/10);
EL=EL_new;
ER=ER_new;

E_save(t,:)=E(:);
Esum(t) = sum(E); % Check for conservation of energy

pause(0.1)


plot(E)
title(num2str(t))

ylim([0 2])
%jnk = input('press enter to continue');

end


