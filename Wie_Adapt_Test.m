% Wie_Adapt_Test.m
% Dynamical Diffraction code Tester
% Requires WieAdapt.m, Wiestep.m, StrainCheck.m, MFP.m
% First written November 9, 2016 by Eric Landahl
% Revised November 18, 2016 and used to benchmark against Sergey

clear all;
more off; % allows data output to continuously scroll during execution

%% Constants for thermal transport / diffussion
T0 = 10; % Initial temperature rise over room temperature in Kelvin
alpha_T = 2.6e-6; % Thermal expansion of Si in 1/Kelvin
k0 = 0.0898; % total diffusivity in um^2/ns
Lambda_MFP = 1e-6; % Phonon Mean Free Path in um   

%% Spatio-temporal grid for diffusion calculations
%n = 100; % number of spatial grid points
%m = 500; % number of temporal points
%T = zeros(n); % temperature initially set to zero everywhere
%Length = 50; % Depth in um.  
%h = Length/n; % spatial step in um
%Duration = 1000; % Time duration in ns
%tau = Duration/m; % spatial step in ns
%fprintf('Diffusion timestep is %d ns and lengthstep is %d um \n',tau,h)
%fprintf('Diffusion duration is %d ns and length is %d um \n',Duration,Length)
%time = tau*(1:m);
%x = h*(1:n);
%fprintf('Calculating diffusion')
%tic;
%[T Tg TotalErr] = MFP (x, time, T0, k0, Lambda_MFP);
%fprintf('\n')
%TotalErr
%toc

%% Testing with step function temperature
n = 1000;
m = 1;
time=1:m;
for nn = 1:n
    x(nn)=0.02*nn;
  for mm = 1:m
    if x(nn) <= 1
      T(nn,mm)=10;
    else
      T(nn,mm)=0;
    end
  end
end

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
tol = 1e-2; % tolerance.  Higher value for more speed and less precision
dz_min = 1.1e-10; % Minimum step size in meters
dz_max = 1e-7; % Maximum step size in meters
f = 10; %Shift factor for convergence

opts(1) = tol;
opts(2) = dz_min;
opts(3) = dz_max;
opts(4) = f;

%% Assemble angle / depth grid
z = x*1e-6; % Depth in meters for WieAdapt.m
dtheta = 0.02*pi/180; % width of angular scan in radians
numpts = 200; % number of angular points
offset = 1e-3*pi/180; % shift center of calculation away from Braggg peak slightly
p = 1:numpts;
theta = thetaB - dtheta*((p-numpts/2)/numpts)+offset; % angles to calculate

%% Loop in time
fprintf('Calculating rocking curves')

tic;
for i = 1:m
  Strain = [alpha_T*T(:,i) 0*T(:,i) 0*T(:,i)];
  [X err Steps Strain_out] = WieAdapt (Strain, z, theta, opts, params);
  X_save(:,i)=X;
  if (i/10 == floor(i/10))
    fprintf('.')
  end
%  figure(11);clf;hold on;
%  plot(theta*180/pi,X_save(:,i).*conj(X_save(:,i)))
%  xlabel('Angle (degrees)')
%  ylabel('Diffraction Intensity')
%  title([num2str(time(i)) ' ns'])
%  hold off;
%  figure(12);clf;hold on;
%  plot(z*1e6,Strain_out(:,1))
%  xlabel('Depth (um)')
%  ylabel('Longitudinal Strain')
%  title([num2str(time(i)) ' ns'])
%  hold off;
%  pause(0.1)
end
fprintf('.\n')
toc

dAngle = delta/sin(2*thetaB)*(180/pi); % Correct for refractive index

figure(5);clf;hold on;
Si_strained = load('Si_strained.dat'); % Data from Sergey's GID_SL server (see NOTES below)
angle1 = Si_strained(:,1); %i
Intensity1 = Si_strained(:,2);
p4=semilogy(angle1-dAngle, Intensity1,':k','LineWidth',4);
p5=semilogy(theta*180/pi,X.*conj(X),'-b','LineWidth',1);
xlabel('Theta (degrees)','FontSize',14)
ylabel('Diffracted intensity','FontSize',14)
AX=legend([p4 p5], 'strained GID-SL', 'strained');
xlim([theta(end)*180/pi theta(1)*180/pi])
ylim([1e-4 1.1])
set(gca,'fontsize',14)
set(AX,'FontSize',14);
hold off


