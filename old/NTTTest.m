% NTTTest.m
% Nanoscale Thermal Transport Tester
% Requires WieAdapt.m, Wiestep.m, StrainCheck.m, MFP.m
% First written November 9, 2016

clear all;
more off; % allows data output to continuously scroll during execution

for LMFP = 1:10

%% Constants for thermal transport / diffussion
T0 = 20; % Initial temperature rise over room temperature in Kelvin
alpha_T = 2.6e-6; % Thermal expansion of Si in 1/Kelvin
k0 = 0.0898; % total diffusivity in um^2/ns
Lambda_MFP = (LMFP-1)/5 + 1e-10; % Phonon Mean Free Path in um   

% Spatio-temporal grid for diffusion calculations
n = 100; % number of spatial grid points
m = 1000; % number of temporal points
T = zeros(n); % temperature initially set to zero everywhere
Length = 50; % Depth in um.  
h = Length/n; % spatial step in um
Duration = 1000; % Time duration in ns
tau = Duration/m; % spatial step in ns
fprintf('Diffusion timestep is %d ns and lengthstep is %d um \n',tau,h)
fprintf('Diffusion duration is %d ns and length is %d um \n',Duration,Length)
time = tau*(1:m);
x = h*(1:n);
fprintf('Calculating diffusion')
tic;
[T Tg TotalErr] = MFP (x, time, T0, k0, Lambda_MFP);
fprintf('\n')
TotalErr
toc

%% Testing with step function temperature
%n = 1000;
%m = 1;
%time=1:m;
%for nn = 1:n
%    x(nn)=0.02*nn;
%  for mm = 1:m
%    if x(nn) <= 1
%      T(nn,mm)=10;
%    else
%      T(nn,mm)=0;
%    endif
%  end
%end

%% Assemble params arrray containing parameters for the crystal diffraction
lambda= 1.23984E-10; % 10 keV
p_0r = -0.97631E-05;
p_0i = -0.14871E-06;
p_Hr = -0.49634E-05;
p_Hi = -0.13789E-06;
delta = 0.48816E-05; % Real part of refractive index
tB_deg = 27.167; % Bragg angle degrees
thetaB= tB_deg*3.14159/180;  % Sergey Si, Bragg angle in radians
phi = 0; % Assymetry angle in radians (zero for symmetric reflections)

params(1) = p_0r + 1i*p_0i;
params(2) = p_Hr + 1i*p_Hi;
params(3) = thetaB;
params(4) = phi;
params(5) = delta;
params(6) = lambda; 

%% Assemble opts array containing options for adaptative depth stepping
tol = 1e-6; % tolerance.  Higher value for more speed and less precision
dz_min = 1.1e-10; % Minimum step size in meters
dz_max = 1e-7; % Maximum step size in meters
f = 2; %Shift factor for convergence

opts(1) = tol;
opts(2) = dz_min;
opts(3) = dz_max;
opts(4) = f;

%% Assemble angle / depth grid
z = x*1e-6; % Depth in meters for WieAdapt.m
dtheta = 0.02*pi/180; % width of angular scan in radians
numpts = 200; % number of angular points
offset = -1e-3*pi/180; % shift center of calculation away from Braggg peak slightly
p = 1:numpts;
theta = thetaB - dtheta*((p-numpts/2)/numpts)+offset; % angles to calculate

X0 = Wiestep (0*theta, theta, z(end), dz_min, [0 0 0], params); % calculate unstrained amplitude
Int0 = X0.*conj(X0);
Cent0 = (180/pi)*sum(Int0.*theta)/sum(Int0);

%% Loop in time
fprintf('Calculating rocking curves')
figure(11);clf;
figure(13);clf;
figure(14);clf;
dAngle = delta/sin(2*thetaB)*(180/pi); % Correct for refractive index
tic;
for i = 1:m
  Strain = [alpha_T*T(:,i) 0*T(:,i) 0*T(:,i)];
  [X err Steps Strain_out] = WieAdapt (Strain, z, theta, opts, params);
  X_save(:,i)=X;
  if (i/10 == floor(i/10))
    fprintf('.')
  end
  Int_save(:,i) = X_save(:,i).*conj(X_save(:,i));
  Cent_save(i) = (180/pi)*sum(Int_save(:,i).*theta')/sum(Int_save(:,i));
  figure(11);hold on;
  semilogy(theta*180/pi,Int_save(:,i))
  xlabel('Angle (degrees)')
  ylabel('Diffraction Intensity')
  title([num2str(time(i)) ' ns'])
  hold off;
  figure(12);clf;hold on;
  plot(z*1e6,Strain_out(:,1))
  plot(z*1e6,Strain(:,1),':')
  xlabel('Depth (um)')
  ylabel('Longitudinal Strain')
  title([num2str(time(i)) ' ns'])
  hold off;
  figure(13);hold on;
  plot(time(i),(Cent_save(i)-Cent0)*1000)
  xlabel('Time (ns)')
  ylabel('Centroid shift (mdeg)')
  figure(14);hold on;
  plot(theta*180/pi,Int_save(:,i))
  xlabel('Angle (degrees)')
  ylabel('Diffraction Intensity')
  title([num2str(time(i)) ' ns'])
  hold off;
  pause(0.1)
end
fprintf('.\n')
toc

fname = [num2str(LMFP) '_' num2str(T0)]
save(fname)
clear -x LMFP;

end



