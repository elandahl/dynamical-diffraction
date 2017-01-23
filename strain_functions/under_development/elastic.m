% elastic.m
%
% For solving elastic wave propagation inside a semiconductor 
% driven by electronic strain
% First written by Eric Landahl 5/17/10, revised 5/19/10 as strain0.m
% Rewritten by EL 1/21/17 
%
% This solves the wave equation u_tt = c^2 u_xx in semi-invinit 1D geometry
% strain_in is the strain from thermal expansion and electron pressure
% strain_out adds in the elastic strain that is driven by strain_in
%
% Assumes linear response of the material at a single longitudinal sound speed
%% INPUTS:
%   crystal         determines x-ray and strain properties, chosen from:
%     GaAs
%     Si
%     Ge
%     InSb
%   t_in                a vector of times in seconds
%   z_in                a vector of depths in meters
%   strain_in     longitudinal strain, size = length(t_in) x length(z_in)
%% OUTPUTS:
%   strain_out    longitudinal strain, size = length(t_in) x length(z_in)
%
%% NOTE: Only longitudinal strain is implemented at this time
%
%% TYPICAL USAGE
% 
%       to be added

function [strain_out,z_out,t_out] = elastic(strain_in,z_in,t_in,crystal)

num_z = 300; % Number of depth points to be used in calculation

%% Load longitudinal sound speed.  
load ../../sample.dat;
ID = find(strcmp({sample.name}, crystal)==1);
v = sample(ID).soundLongitudinal.val; % in um/ns
v = v*1000; % Convert sound speed to m/s

%% Remesh
dz = z_in(end)/num_z;
dt = dz/v;
z=dz:dz:z_in(end); % spatial grid
t=(1:length(z))*dt;
num_t = length(t); % number of time points
s = zeros(num_t,num_z); % Grid for strain
fprintf('Original grid contained %i spatial and %i temporal points.\n', length(z_in),length(t_in))
fprintf('Remeshing grid for elastic calculation to %i spatial and %i temporal points.\n',num_z,num_t)
fprintf('The new grid has evenly spaced points of %.1f nm and %.1f ps.\n',dz*1e9,dt*1e12) 
[t_in_meshgrid,z_in_meshgrid] = meshgrid(t_in,z_in');
[t_meshgrid,z_meshgrid] = meshgrid(t,z');
strain_remesh = interp2(t_in_meshgrid,z_in_meshgrid,strain_in',t_meshgrid,z_meshgrid,'spline');
s0 = strain_remesh';

fprintf('Done remeshing strain input.  Strarting strain propagation.\n')
s (1,:) = s0(1,:);
s1 = zeros(num_t,num_z);
s2 = zeros(num_t,num_z);


for j = 1:num_t
  for i = 1:num_z
    if i-j < 0
      s1(j,i) = -s(1,j-i)/2;
    elseif i-j > 0
      s1 (j,i) = s(1,i-j)/2;
    else
      s1 (j,i) = 0;
    end
  end
end

for j = 1:num_t
  for i = 1: num_z
    if i + j > num_z  
      s2 (j,i) = 0;
    else
      s2 (j,i) = s(1,i+j)/2;
    end
  end
end



%% WORKS!!!
%for i = 1:num_z
%  for j = 1:num_t
%    if i-j < 0
%      s1(j,i) = -s(1,j-i)/2;
%    elseif i-j > 0
%      s1 (j,i) = s(1,i-j)/2;
%    else
%      s1 (j,i) = 0;
%    end
%  end
%end
%
%for i = 1: num_z
%  for j = 1:num_t
%    if i + j > num_z  
%      s2 (j,i) = 0;
%    else
%      s2 (j,i) = s(1,i+j)/2;
%    end
%  end
%end

 
   
s = s1 + s2;



%for i = 2: num_t
%  s1 (i,2:num_z) = s (i-1,1:num_z-1)/2; % Wave headed into bulk
%  s1 (i,1) = -s (i-1,1)/2; % Flips sign on reflection from surface
%  s2 (i,1:num_z-1) = s (i-1,2:num_z)/2; % Wave headed towards surface
%  s2 (i,num_z) = 0; % no wave in from bulk
%end
%fprintf('Strain propagation done.  Remeshing to original grid.\n')
%strain_remesh = interp2(t_in_meshgrid,z_in_meshgrid,strain_in',t_meshgrid,z_meshgrid,'spline');


%% IDEA:  USE LBM ?
%
%% NOTE:  Should vectorize this for speed
%for k = 1:num_t  % k is the time shift so driving strain is treated as IC
%    clear t0;
%    for m = 1:num_t-k+1      % This loop shifts the time array.  
%        t0(m) = t(k+m-1);     % Note that space (depth) is not shifted.
%    end
%    s0 = zeros(num_t-k+1,num_z); % Shrinks driving strain to match t
%    s0(1,:) = strain_in(k,:);  % Makes kth driving strain Initial Condition
%    [stemp] = strain0(s0,num_z,num_t); %2D strain wave from kth driving term
%    for m = 1:num_t-k  % Shift time back
%        sold(m+k,:)=s(m+k,:);
%        s(m+k,:)=stemp(m,:) + sold(m+k,:); % Add strain waves @ unshifted time
%    end
%end

strain_out = s;  %% CHANGE BACK TO s out after trouble shooting
z_out =z*1e6;
t_out =t*1e9;

end      