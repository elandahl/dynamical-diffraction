% Debye.m
% Determines the thermal energy of a material given the temperature
% First written by Eric Landahl, 12-21-2016
% Possibly could be sped up using arrrayfun 
%
% INPUTS
%   T is the Temperature in Kelvin
%   Theta is the Debye Temperature of the materi in Kelvin
%   e.g. Theta = 640 for Silicon
% OUTPUTS
%    Energy is the Phonon Gas Energy Density, or Energy/(N*kB)

function [Energy] = Debye (T, Theta)
kB = 1.38064852E-23; % Boltzmann's Constant in J/K

for i = 1:length(T)
F = @(z)(z.^3)./(exp(z)-1);
Energy(i) = (9*T(i)).*(T(i)/Theta).^3*quad(F,0,Theta./T(i));
end


end
