clear all;
close all;
%% load material spectral data
% 
% % gold material
% hagemann_gold =  importdata('Hagemann_gold.txt',' ',8);
% jc_gold = importdata('Johnson_Christy_gold.txt',' ',9);
% olmon_gold = importdata('Olmon-ev_gold.txt',' ',9);
% 
% % define epsilon function (micron input wavelength)
% eps_au = @(lambda) interp1(olmon_gold.data(:,1),(olmon_gold.data(:,2)+1i*olmon_gold.data(:,3)).^2,lambda,'pchip','extrap');
% 
% wavelength=1e-2:1e-2:50;
% 
% % and plot raw data n,k 
% figure(1);
% loglog(hagemann_gold.data(:,1),hagemann_gold.data(:,2),'.-',hagemann_gold.data(:,1),hagemann_gold.data(:,3),'.-', ...
%        jc_gold.data(:,1),jc_gold.data(:,2),'.-',jc_gold.data(:,1),jc_gold.data(:,3),'.-', ...
%        olmon_gold.data(:,1),olmon_gold.data(:,2),'.-',olmon_gold.data(:,1),olmon_gold.data(:,3),'.-', ...
%        wavelength,real(sqrt(eps_au(wavelength))),'-',wavelength,imag(sqrt(eps_au(wavelength))),'-');
% xlim([0.1 25]);

%% load file all units in microns

load('beamsplitter_scan_data_period__2.980000E+00__thickness__1.040000E-01__.mat');
% for now: manual paramter setting
period = 2.98;
gold_thickness = 0.119;

%%

titan_idx = 6;

T0p = THy_ZO(:,titan_idx,1,end);
Tp  = THy(:,titan_idx,1,end);
T0s = TEy_ZO(:,titan_idx,1,end);
Ts  = TEy(:,titan_idx,1,end);

% save(['./bs_scan_period__',int2str(period*1e+3),'__thickness__',int2str(gold_thickness*1e+3),'__titan_width__',num2str(au_ti_widths(titan_idx,2)*1e+3),'__transmission'],'wavelength','T0p','Tp','T0s','Ts')


%%

figure(1);
clf;

plot(wavelength,THy_ZO(:,titan_idx,1,end),'.-','LineWidth',2,'Markersize',20); 
hold on;
plot(wavelength,THy(:,titan_idx,1,end),'.-','LineWidth',2,'Markersize',20);
plot(wavelength,TEy_ZO(:,titan_idx,1,end),'.-','LineWidth',2,'Markersize',20);
plot(wavelength,TEy(:,titan_idx,1,end),'.-','LineWidth',2,'Markersize',20);

legend('max T_{0p}','max T_{p}','max T_{0s}','max T_{s}');
title(['Beamsplitter transmission w_{Ti} = ',num2str(au_ti_widths(titan_idx,2)),' \mum'],'fontsize',20);
xlabel('wavelength \lambda (\mum)','fontsize',20);
ylabel('Transmission T','fontsize',20);
set(gca,'fontsize',20,'linewidth',2);

xlim([2.5,25]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 20*3/4])
set(gcf,'papersize',[30,20]);

drawnow;
% print('-dpng',['./bs_scan_period__',int2str(period*1e+3),'__thickness__',int2str(gold_thickness*1e+3),'__titan_width__',num2str(au_ti_widths(titan_idx,2)*1e+3),'__transmission.png'],'-r100');
%print('-dpng',['./pictures/bs_scan_period__',int2str(period*1e+3),'__thickness__',int2str(thickness*1e+3),'__errors.png'],'-r100');


