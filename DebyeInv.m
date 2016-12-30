% DebyeInv.m
% Determines the temperature of a material given the Debye Energy
% Requires Debye.m
% First written by Eric Landahl, 12-21-2016
% Possibly could be sped up using arrrayfun 
%
% INPUTS
%   E is the Temperature/(kB*N) in Joules
%   Theta is the Debye Temperature of the materi in Kelvin
%   e.g. Theta = 640 for Silicon
% OUTPUTS
%    TT is temperature in K

function [TT] = DebyeInv (Energy, Theta)

for i = 1: length(Energy)
  myfun = @(T)(Debye(T,Theta)-Energy(i)); % Define a function that goes to zero at the desired energy
  TT(i) = fzero(myfun,Theta); % Find the temperature at that energy
end
  
end
