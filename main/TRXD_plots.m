% TRXD_plots.m
% Produces plots of strain components and rocking curves
% May be called from either within of after TRXD.m
% By E.L.  Last revised 1/23/17 to workaround semilogy and MatLAB bug


function [Intensity centroid FWHM] = TRXD_plots (A,A0,time,angle,Strain,z,ang_res,plot_opts)

z = z*1e6; % Convert depth to um for plotting
time = time*1e9; % Convert time to ns for plotting
angle = angle*1e3; % Convert degrees to mdeg for plotting
ang_res = ang_res*1e3/(2*sqrt(2*log(2))); % Angular resolution for blurring

trans_status = max(max(Strain(:,:,2)))>1e-10; % 0 if no transverse strain
sheer_status = max(max(Strain(:,:,2)))>1e-10; % 0 if no sheer strain

if strcmp(plot_opts,'animate')==1
  fprintf('Animating plots at each timepoint.\n')
  figure(1);clf;
  figure(2);clf;
  figure(3);clf;
  figure(6);clf;
  figure(7);clf;
end

if strcmp(plot_opts,'animate')==1 && trans_status == 1
  figure(4);clf;
end

if strcmp(plot_opts,'animate')==1 && trans_status == 1
  figure(5);clf;
end

%% Unstrained intensity
  Int0 = A0.*conj(A0); % unstrained crystal x-ray intensity
  centroid0 = sum(Int0.*angle)/sum(Int0); % centroid of unstrained crystal

%% Plot at each time
for i = 1:length(time)

  longitudinal = Strain(i,:,1);
  transverse = Strain(i,:,2);
  sheer = Strain(i,:,3);
  Int = A(i,:).*conj(A(i,:)); % x-ray intensity
  centroid(i) = sum(Int.*angle)/sum(Int) - centroid0;
  f = (1/(ang_res*sqrt(2*pi)))*exp(-((angle-centroid0)/(sqrt(2)*ang_res)).^2);
  Int = conv(Int,f,'same');
  Int = Int/(max(Int));
  Intensity(i,:)=Int; % return convolved intensities
  FWHM(i) = 2*sqrt(2*log(2))*sqrt((1/(length(angle)-1))*sum(Int.*(angle-centroid(i)).^2));

  if strcmp(plot_opts,'animate')==1
  
  figure(1)
    plot(angle,Int);hold on;
    title([ num2str(time(i)) ' ns']);
    xlabel('Angle (deg)');
    ylabel('Diffraction Intensity')
  hold off; 
  
  figure(2)
    semilogy(angle,Int);hold on;
    title([ num2str(time(i)) ' ns']);
    xlabel('Angle (deg)');
    ylabel('Diffraction Intensity')
    set(gca,'YScale','log')
  hold off; 
  
  figure(3)
    plot(z,longitudinal);hold on;
    title([ num2str(time(i)) ' ns']);
    xlabel('Depth (um)');
    ylabel('Longituidnal Strain')
  hold off;  
  
  % Only plot transverse strain if it is nonzero
  if trans_status == 1 
    figure(4)
     plot(z,transverse);hold on;
     title([ num2str(time(i)) ' ns']);
     xlabel('Depth (um)');
     ylabel('Transverse Strain')
    hold off;  
  end
  
  % Only plot sheer strain if it is nonzero
  if sheer_status == 1
    figure(5)
      plot(z,sheer);hold on;
      title([ num2str(time(i)) ' ns']);
      xlabel('Depth (um)');
      ylabel('Sheer Strain')
    hold off;  
  end 
  
  figure(6);hold on;
    plot(time(i),centroid(i),'o');
    xlabel('Time (ns)');
    ylabel('Centroid (mdeg)');
  hold off

  figure(7);hold on;
    plot(time(i),FWHM(i),'o');
    xlabel('Time (ns)');
    ylabel('FWHM (mdeg)');
  hold off
  
  pause(0.5); % To allow plots to animate smoothly
  
  end
  
  end % End Time loop
  
  if strcmp(plot_opts,'none')~=1
  
  figure(10);clf; hold on;
  
  subplot(2,2,1)
    plot(time,centroid)
    xlabel('Time (ns)')
    ylabel('Centroid (mdeg)')
  
  subplot(2,2,3)
    plot(time,FWHM)
    xlabel('Time (ns)')
    ylabel('FWHM (mdeg)')
  
%  subplot(2,2,3)
%    surf(time,angle,(Intensity'),'LineStyle','none')
%    ylabel('Angle (mdeg)')
%    xlabel('Time (ns)')
%    title('Linear Scale')
%    view(2)
%    grid off
%    ylim([min(angle)/5 max(angle)/3])
%    xlim([time(1) time(end)])
%  
%  subplot(2,2,3)
%    surf(time,angle,(log(Intensity))','LineStyle','none')
%    ylabel('Angle (mdeg)')
%    xlabel('Time (ns)')
%    title('Log Scale')
%    view(2)
%    grid off
%    ylim([min(angle)/1 max(angle)/6])
%    xlim([time(1) time(end)])

ii(1)=1;  % Choose these timepoints
ii(2)=4;
ii(3)=10;
ii(4)=20;
ii(5)=25;
ii(6)=30;
time = time*1e-9; % convert to ns for plotting

subplot(2,2,2);hold all;
for j = 1:6
ti = ii(j);
plot(angle,Intensity(ii(j),:))
end
xlabel('Angle (mdeg)')
ylabel('X-ray Intensity')
xlim([-1 4.5])
ti=1;
lgd = sprintf('%.0e s', time(ti));
for idx=2:6, ti = ii(idx); lgd = strvcat(lgd, sprintf('%.0e s', time(ti))); end,
lgd=cellstr(lgd);
%for idx=1:6, disp(lgd{idx}), end
LEG=legend(lgd);
%set(gca, 'FontSize', 16)
set(LEG,'FontSize',8)


subplot(2,2,4);hold all;
for j = 1:6
ti = ii(j);
semilogy(angle,Intensity(ii(j),:))
set(gca,'YScale','log')
end
xlabel('Angle (mdeg)')
ylabel('X-ray Intensity')
ylim([1e-4 1.1])
xlim([-6 9])
set(gca,'YScale','log')


ti=1;
lgd = sprintf('%.0e s', time(ti));
for idx=2:6, ti = ii(idx); lgd = strvcat(lgd, sprintf('%.0e s', time(ti))); end,
lgd=cellstr(lgd);
%for idx=1:6, disp(lgd{idx}), end
LEG=legend(lgd);
%set(gca, 'FontSize', 16)
set(LEG,'FontSize',8)
set(gca,'YScale','log')

  end
  
end

