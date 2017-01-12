% strainFile.m
% Feeds strain from an external file into TRXD
% First written by Eric Landahl, 1.10.2017
% Usually called by TRXD.m
%
%% INPUTS:
%   filename    presumed to longitudinal only, strain(depth,time)
%   dz          depth step, in meters
%   dt          time step, in seconds
%% OUTPUTS:
%   longitudinal    longitudinal strain, size = length(time_out) x length(z)
%   transverse      transverse strain, size = length(time_out) x length(z)
%   sheer           sheer strain, size = length(time_out) x length(z)
%   z               a vector of depths in meters
%   time_out        a vector of times returned in seconds
%
%% NOTE ALSO: Three strain components are returned, but usually only longitudinal
%             is non-zeros
%
%% TYPICAL USAGE
% 
% [st1 st2 st3 time_out z] = strainFile ('Soo1.txt', 10e-9, 10e-12, 1e-5);
%
function [longitudinal transverse sheer time_out z] = strainFile(filename,dz,dt)

% Read in strain file
filename
longitudinal = dlmread(filename);
[numtimes numdepths] = size(longitudinal);
    
% Assemble depth array
z = dz*(1:numdepths);

% Assemble time array
time_out = dt*(1:numtimes);

% Set other strain components to zero
transverse =0.*longitudinal; % No transverse strain
sheer = 0.*longitudinal; % No sheer strain  
  