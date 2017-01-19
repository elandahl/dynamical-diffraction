% WieAdapt.m
% [X X0 err Steps] = WieAdapt (Strain, z, theta, opts, params)
%
% Function to take adaptive depth steps
% to propagate dynamical x-ray diffraction
% in a strained crystal up from the unstrained bulk
% using Wie's algorithim in J. Appl. Phys. 59, 3743 (1986)
% Written by Eric Landahl, 11/6/2016 
% Revised by EL 1/9/17 to also output unstrained scattering X0
% Last revised by EL 1/16/17 to improve convergence at high strain
% Improved accuracy by fixing delta and final step interpolation 1/19/2017
% For more details see help for Wiestep.m
%
% REQUIRES THESE FUNCTIONS:
%   Wiestep.m to take the actual individual steps
%   StrainCheck.m to make sure that the strain and parameters are OK
%
% INPUTS:
%   Strain  N x 3 matrix with N rows of [longitudinal transverse sheer]
%   z       N x 1 matrix with N rows of depths used in Strain, in meters
%   theta   M x 1 matrix containing angular points, in radians
%   opts    1 x 4 matrix of options [tol dz_min dz_max f]
%           tol     numerical tolerance used to set step sizes, e.g. 1e-3
%           dz_min  minimum depth step size, in meters, e.g. 1e-10
%           dz_max  maximum depth step size, in meters, e.g. 1e-7
%           f       shift factor, or step size change each iteration, e.g. 2
%   params  1 x 5 matrix of diffraction params: [psi0 psiH thetaB delta lambda]
%           psi0    complex
%           psiH    complex
%           thetaB  Bragg angle, in radians
%           phi     Assymetry angle, in radians
%           delta   real refractive index - 1 (~0, for small Bragg correction) 
%           lambda  x-ray wavelength, in meters
%
% OUTPUTS
%   X       M x 1 complex scattering amplitude. Intensity = X.*conj(X)
%   X0      Complex scattering amplitude of unstrained crystal
%   err     maximum error that occured during the calculation
%   Steps   number of steps required to calculate X
%
% NOTES
%   If f = 1 or any of the other opts are zero, the step size is fixed 
%   The initial step size is dz_min
%   Strain and z must agree, but do not need to be evenly spaced
%   z(end) must be > 5 X L_ext, and meet other criteria of StrainCheck.m
%   delta provides a small shift from dynam. diff. theory applied everywhere
%   zz decreases heading to the surface: zz + dz/2 is deeper than zz
%   zz is the incremented adapted depth; z is the input strain depth



function [X X0 err Steps Strain] = WieAdapt (Strain, z, theta, opts, params)

%% Unpack opts
tol = opts(1);
dz_min = opts(2);
dz_max = opts(3);
f = opts(4);

%% Unpack params
psi0 = params(1);
psiH = params(2);
thetaB = params(3);
phi = params(4);
delta = params(5);
lambda = params(6);

%% Check strain file and other inputs.  Stop if there is a problem. Taper Strain.
[warnings Strain opts] = StrainCheck (Strain, z, theta, opts, params);
if (length(warnings)~=0)
  warnings
  fprintf('\n  Ending program and returning null values.\n')
  X = 0; err = 0; Steps = 0;
  return
end

%% Set the maximum depth to 5 times the extinction length
gamma0 = sin(thetaB+phi); 
gammaH = sin(thetaB-phi);
Lext = lambda*sqrt(abs(gammaH)*gamma0)/(pi*abs(real(psiH)));
z0 = 4.9*Lext; % "infinite" crystal depth in meters (with a 0.05 safety gap)
zz = z0; % start at z0.  zz will be an updated value; z is the imported array

%% Start with unstrained crystal by putting a null array in for X_in
dz = dz_min; % begin with smallest possible step size
X_in = Wiestep (0*theta, theta, zz, dz, [0 0 0], params); % calculate unstrained amplitude
X0=X_in; % Save for later
intX0 = X0.*conj(X0); % Substrate intensity
max0 = max(intX0); % Maximium substrate intensity
max_errs=0;
X_half = X0; % Just to start

%% Adaptative stepping main loop
i = 1; % iteration variable 
while (zz>=0);
  st = interp1(z, Strain, zz,'spline', 'extrap');   % st is interpolation of Strain at zz
  X_out = Wiestep (X_in, theta, zz, dz, st, params); % take a full step at once
  zz_half = zz + dz/2; % half step depth (deeper by a half step)
  st_half = interp1(z, Strain, zz_half,'spline', 'extrap'); % st_half is interpolation of Strain at zz_half
  X_half = Wiestep (X_in, theta, zz_half, dz/2, st_half, params); % take a half step
  X_full = Wiestep (X_half, theta, zz, dz/2, st, params); % take a second half step
  intX_full = X_full.*conj(X_full); % intensity after two half steps
  intX_out = X_out.*conj(X_out); % intensity after one full step
  errs = sum(abs(intX_full-intX_out))/sum(intX_full); % error is the difference in intensities
  max_errs = max(errs,max_errs); % save maximum error so far
  zz_old =zz_half; % To use for corrections on final step
  if (errs<tol) && (dz<dz_max) % If the error is small then
    zz = zz - dz; % Move to shorter depth (closer to the surface)
    dz = dz*f; % increase dz if possible
    X_in = X_full; % Average full and two half-step results
  elseif  (errs<tol) && (dz>=dz_max) % If error is small but already at dz=dz_max
    zz = zz - dz; % Move to shorter depth (closer to the surface)
    X_in = X_full; % Average full and two half-step results
  elseif (errs>tol) && (dz<=dz_min) % If the error is large but dz=dz_min
    X_in = X_full; % Average full and two half-step results
    zz = zz - dz; % Move to shorter depth (closer to the surface)
  elseif (errs>tol) && (dz>dz_min) % If the error is large then
    dz = dz/(5*f); % shrink dz (don't record X, instead repeat)
  end  

% Optional plots to watch the algorithim work
%  err_save(i) = errs; % save error level
%  dz_save(i)=dz;
%  st_save(i)=st(1);
%  st_half_save(i)=st_half(1);
%  zz_save(i)=zz;
%%    figure(20);clf; hold on;
%%      plot(theta*180/pi,X_out.*conj(X_out),'-r','LineWidth',4)
%%      plot(theta*180/pi,X_full.*conj(X_full),'-b','LineWidth',1)
%%      plot(theta*180/pi,X0.*conj(X0),'-k')
%%      xlabel("Theta (degrees)")
%%      ylabel("Diffracted intensity")
%%      title(['Depth: ' num2str(zz*1e6) ' um, Step: ' num2str(dz*1e9) ' nm, Strain: ' num2str(st_half(1))])
%%    hold off;
%%    figure(21);clf; hold on;
%%      semilogy(theta*180/pi,X_out.*conj(X_out),'-r','LineWidth',4)
%%      semilogy(theta*180/pi,X_full.*conj(X_full),'-b','LineWidth',1)
%%      semilogy(theta*180/pi,X0.*conj(X0),'-k')
%%      xlabel("Theta (degrees)")
%%      ylabel("Diffracted intensity")
%%      title(['Depth: ' num2str(zz*1e6) ' um, Step: ' num2str(dz*1e9) ' nm, Strain: ' num2str(st_half(1))])
%%    hold off;
%%    figure(22);clf;hold on;
%%      plot(zz_save*1e6,st_save,'or','LineWidth',4)
%%      plot(zz_save*1e6,st_half_save,'xb','LineWidth',1)
%%      xlabel('Depth (um)')
%%      ylabel('Longitudinal Strain')
%%    hold off;
%%    figure(23);clf;hold on;
%%      plot(zz_save*1e6,dz_save*1e9)
%%      xlabel('Depth (um)')
%%      ylabel('Step size (nm)')
%%    hold off;
%%    figure(24);clf;hold on;
%%      plot(zz_save*1e6,err_save)
%%      xlabel('Depth (um)')
%%      ylabel('Error')
%%    hold off;
%%    pause(0.1)
%% End of optional plots
  i = i + 1; % increment up number of iterations
end

%% Return values
f1 = zz_old/(zz_old-zz); % Correction factor for final depth step
X =( f1*X_half + (1-f1)*X_full); % Correct final step
X = X/max(X.*conj(X)); % Normalize
err = max_errs; % return maximum error
Steps = i-1; % return number of steps

% Optional plots (after algorithm)
%    figure(20);clf; hold on;
%      plot(theta*180/pi,X_out.*conj(X_out),'-r','LineWidth',4)
%      plot(theta*180/pi,X_full.*conj(X_full),'-b','LineWidth',1)
%      plot(theta*180/pi,X0.*conj(X0),'-k')
%      xlabel("Theta (degrees)")
%      ylabel("Diffracted intensity")
%      title(['Depth: ' num2str(zz*1e6) ' um, Step: ' num2str(dz*1e9) ' nm, Strain: ' num2str(st_half(1))])
%    hold off;
%    figure(21);clf; hold on;
%      semilogy(theta*180/pi,X_out.*conj(X_out),'-r','LineWidth',4)
%      semilogy(theta*180/pi,X_full.*conj(X_full),'-b','LineWidth',1)
%      semilogy(theta*180/pi,X0.*conj(X0),'-k')
%      xlabel("Theta (degrees)")
%      ylabel("Diffracted intensity")
%      title(['Depth: ' num2str(zz*1e6) ' um, Step: ' num2str(dz*1e9) ' nm, Strain: ' num2str(st_half(1))])
%    hold off;
%    figure(22);clf;hold on;
%      plot(zz_save*1e6,st_save,'or','LineWidth',4)
%      plot(zz_save*1e6,st_half_save,'xb','LineWidth',1)
%      xlabel('Depth (um)')
%      ylabel('Longitudinal Strain')
%    hold off;
%    figure(23);clf;hold on;
%      plot(zz_save*1e6,dz_save*1e9)
%      xlabel('Depth (um)')
%      ylabel('Step size (nm)')
%    hold off;
%    figure(24);clf;hold on;
%      plot(zz_save*1e6,err_save)
%      xlabel('Depth (um)')
%      ylabel('Error')
%    hold off;
%    pause(0.1)
% End optional plots

end
