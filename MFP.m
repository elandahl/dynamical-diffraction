function [T Tg TotalErr] = MFP (x, time, T0, k0, Lambda)

%% MFP.m
% Crank Nicolson PDE solver
% For nanoscale thermal transport using simple Mean Free Path (MFP) model
% MFP model described in T.S. Fisher, "Thremal Energy at the Nanoscale" (2014)
% C-N algorithm:  LHS * T(n+1) = RHS, or T(n+1) = inv(LHS) * RHS
% Can optionally use tri_diag.m for matrix inversion, but seems inv() is faster
% Boundary condition is constant temperature deep in crystal
% Initial condition is zero everywhere except T0 initially at surface
% First written by Eric Landahl, October 5, 2016
% Revised and checked agains macroscopic case October 6, 2016
% Turned into a function SiMFP.m by EL on October 10, 2016
% Generalized, streamlined, and revised by EL  as MFP.m October 12, 2016

%% Inputs
% x is an equally spaced array of spatial gridpoints
% time is an equally spaced array of temporal gridpoints
% k0 is the bulk thermal conductivity, in units of length^2/time

%% Outputs
% T is the temperature evaluated on the x,time grid
% Tg is the analytic solution for the "Gaussian Kernel" (MFP = 0)
% Err is the sum total difference between T and Tg at all space and time
%     and Err should approach zero for Lambda=0

%% Spatio-temporal range
%  array convention is T(space,time) = T (n,m)
%  n => space, m => time
a = k0;
h = x(2)-x(1); % Uniform spatial steps
tau = time(2) - time(1); % Uniform temporal steps
n_max = length(x);
m_max = length(time);
T = zeros(n_max); % temperature initially set to zero everywhere
r = a*tau/(2*h^2); % coefficient in Crank-Nicholson matrix

T(1,1) = T0; % Initial Condition
r = (erf(x(1)/Lambda))*a*tau/(2*h^2); % coefficient in Crank-Nicholson matrix
%% Compute Left Hand Side matrix
LHS(1,1)=1+2*r;
LHS(1,2)=-2*r;
for n = 2:(n_max-1)
  r = (erf(x(n)/Lambda))*a*tau/(2*h^2); % coefficient in Crank-Nicholson matrix
  LHS(n,n-1) = -r;
  LHS(n,n) = 1+2*r;
  LHS(n,n+1) = -r;
end
LHS(n_max,n_max)=1+2*r;
LHS(n_max,n_max-1)=-2*r;

%% Loop over time
for m = 1:m_max-1 % time loop
  r = (erf(x(1)/Lambda))*a*tau/(2*h^2); % coefficient in Crank-Nicholson matrix
  RHS(1) = T(1,m) + r*(2*T(2,m) - 2* T(1,m));
  for n = 2:n_max-1 % Loop over space
    r = (erf(x(n)/Lambda))*a*tau/(2*h^2); % coefficient in Crank-Nicholson matrix
    RHS(n) = T(n,m) + r*(T(n+1,m) - 2*T(n,m) + T(n-1,m));
  end
  RHS(n_max) = T(n_max,m) + r*(-2*T(n_max,m) + 2*T(n_max-1,m));
  T(:,m+1) = inv(LHS)*RHS';  %Use this slower matrix inversion for non tri-diagonal matrices
  %T(:,m+1) = tri_diag(LHS,RHS); % Requires tri_diag.m fast matrix inversion for tri-diagonal matrices
  if (m/100 == floor(m/100))
    fprintf('.'); % Put a period on the output line every 10 percent
  end
end

%% Analytic solution
for m = 1:m_max
sigma(m) = sqrt(2*a*time(m)); % width of Gaussian heat pulse
  for n = 1:n_max
    % Tg, or Gaussian kernel solution is used 
    % Also, it is normalized against the step size, h
    Tg(n,m) = h*(T0/sigma(m))*(1/sqrt(2*pi))*exp(-(x(n)^2)/(2*sigma(m)^2));
   end
end
%
TotalErr=sum(sum(abs(Tg-T)))/sum(sum(Tg));

%% simple plots
%figure(1)
%clf
%hold on
%x=x/1000; % Convert to plotting in us
%p1=plot(x,T(:,1),'-r');
%p2=plot(x,T(:,floor(M/4)),'-m');
%p3=plot(x,T(:,floor(M/2)),'-g');
%p4=plot(x,T(:,floor(3*M/4)),'-b');
%p5=plot(x,T(:,M-1),'-k');
%plot(x,Tg(:,1),':r')
%plot(x,Tg(:,floor(M/4)),':m')
%plot(x,Tg(:,floor(M/2)),':g')
%plot(x,Tg(:,floor(3*M/4)),':b')
%plot(x,Tg(:,M-1),':k')
%legend([p1 p2 p3 p4 p5], [num2str(time(1)/1000) ' ns'], [num2str(time(floor(M/4))/1000) ' ns'], [num2str(time(floor(M/2))/1000) ' ns'], [num2str(time(floor(3*M/4))/1000) ' ns'], [num2str(time(M)/1000) ' ns'])
%xlim([0 x(end)])
%% ylim([-T(1,floor(M/4))*.05 T(1,floor(M/4))*1.5])
%xlabel('Depth (\mum)')
%ylabel('Temperature Change (K)')
%hold off
%
%
%%figure(2)
%%clf
%%surf(T)

% Check that total integrated temperature remains constant
%figure(3);clf;hold on;
%plot(time,sumT,'o')
%xlabel('Time (ps)')
%ylabel('Total integrated temperature')
%
%figure(4);
%plot(x(:),T(:,1:10))

endfunction
