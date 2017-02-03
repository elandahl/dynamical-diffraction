% Larson.m
% for algorithm of Larson and Barnhost, J.A.P. 51 (6), 1980.
% EL, 1/31/17 

pkg load odepkg;
clear all;
more off; % Turn on scrolling

%% Include subdirectories in path
addpath('../include','../strain_functions','../data');

crystal = 'GaAs';
energy = 10; % keV
phi = 0; % assymetry angle in radians

% Make a dummy strain function.  Eventually have this passed in
z = 1e-9 * (1:3000); 
strain = 0 * exp(-z / 1e-6);


%% Load crystal parameters

% Load X0h data
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

% Interpolate X0h data
% Note that X0h(:,1) are the energies in keV
p_0r = -abs(interp1(X0h(:,1), X0h(:,2), energy, 'spline', 'extrap'))
p_0i = abs(interp1(X0h(:,1), X0h(:,3), energy, 'spline', 'extrap')) 
p_Hr = -abs(interp1(X0h(:,1), X0h(:,4), energy, 'spline', 'extrap'))   
p_Hi = -abs(interp1(X0h(:,1), X0h(:,5), energy, 'spline', 'extrap'))   
delta = interp1(X0h(:,1), X0h(:,6), energy, 'spline', 'extrap');   
tB_deg = interp1(X0h(:,1), X0h(:,7), energy, 'spline', 'extrap'); 
thetaB= tB_deg*3.14159/180;  % Sergey Si, Bragg angle in radians
lambda= energy*1.23984E-11; % convert energy in keV to wavelength in meters
% Assemble params arrray containing parameters for the crystal diffraction
params(1) = p_0r + 1i*p_0i;
params(2) = p_Hr + 1i*p_Hi;
params(3) = thetaB;
params(4) = phi;
params(5) = delta;
params(6) = lambda; 

g =  -p_0i / p_Hr
k = p_Hi / p_Hr

% Set angle theta (eventually loop)
tic;
num_angles = 50;
for i = 1:num_angles
theta = (i-num_angles/2)*thetaB/(2000*num_angles) + thetaB; % Angle in radians, near but not equal to Bragg angle

A_in = - p_Hr * pi * z / lambda; % normalized depth array
y0 = (sin(2*thetaB)/p_Hr) * (theta - thetaB + strain(end) * tan(thetaB)) - p_0r/p_Hr;

X0_out(i) = X0(y0, g, k);
theta_out(i) = theta;

y = (sin(2*thetaB)/p_Hr) * (theta - thetaB + strain * tan(thetaB)) - p_0r/p_Hr;

X_in(1) = real(X0_out(i));
X_in(2) = imag(X0_out(i));

A_span = [max(A_in) 0]; 
dA = A_in(2)-A_in(1);

opts = odeset('RelTol',1e-3,'AbsTol',1e-3,'InitialStep',dA,'MaxStep',dA*10000);

[A_out X_out] = ode45( @(A,X) dX(A, X, y0, A_in, g, k), A_span, X_in, opts);

X_save(i) = X_out(end);

end
toc;

figure(1);clf;hold on; 
  semilogy(theta_out*180/pi,X0_out.*conj(X0_out),'b')
  semilogy(theta_out*180/pi,X_save.*conj(X_save),'r')
hold off;



%
%
%
%
%Intensity = X_out(1).^2 + X_out(2).^2; % Magnitude squared of complex number
%
%z_out = - lambda * A_out / (p_Hr * pi); % Convert from normalized to depth in m
%
%figure(1);clf;hold on;
%  plot(z_out,Intensity)
%  xlabel('Depth (m)')
%  ylabel('Intensity')
%hold off;
%
%
%
%
%
