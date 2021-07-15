clear all;
close all;

addpath('./fmm_source/');
%% helping functions
ZeroFunc = @(x) zeros(size(x));
OneFunc = @(x) ones(size(x));
BoxFunc = @(x,a,b) 0.5*(sign(x-a)-sign(x-b));

%% load material spectral data

% silver parameters Drude-Lorentz model
epsBgAg = 2.399296520158095;
dl1 = 0.44054073707611147;
wl1 = 0.0268977;
yl1 = 0.007019091816101157;
dl2 = 0.3642722518396299;
wl2 = 0.022091;
yl2 = 0.0021265004581278437;
wd1 = 0.0419425;
yd1 = 0.0005128576119935657;

% inline function calls for Drude and Lorentz model
% Drude model
chiDrude = @(w, wd, yd) wd.^2./(w.^2+1i.*w.*yd);
% Lorentz model
chiLorentz = @(w, dl, wl, yl) (wl.^2.*dl)./(wl.^2-w.^2-1i.*yl.*w);

% eps_Ag(omega) = eps_infty + chiDrude + chiLorentz_1 + chiLorentz_2
epsAg = @(w) epsBgAg + chiLorentz(w, dl1, wl1, yl1) + chiLorentz(w, dl2,wl2, yl2) - chiDrude(w, wd1, yd1);

% eps_SiO2
epsilonSiO2 = 2.18217;

%% potential omega loop here
% wavelength=10:1:25;
wavelength=100:5:700;
omega = 2*pi./wavelength;

% incAngle = 67.2;
incAngle = 60;
% incAngle = 0:2:89.99;
% incAngle = 0:2:89.9;
% incAngle = 29.5:0.01:30.5
% covarage_factor = 0.0:0.025:1.0;
covarage_factor = 0.5;
% sintheta = sin(pi.*incAngle./180);
fillling_factor = 0.7;
hmm_thickness = 20;
%%
time_in_smat = zeros(length(omega),length(incAngle));

TEy = zeros(length(omega),length(incAngle));
REy = zeros(length(omega),length(incAngle));
THy = zeros(length(omega),length(incAngle));
RHy = zeros(length(omega),length(incAngle));

TEy_Oth_order = zeros(length(omega),length(incAngle));
REy_Oth_order = zeros(length(omega),length(incAngle));
THy_Oth_order = zeros(length(omega),length(incAngle));
RHy_Oth_order = zeros(length(omega),length(incAngle));


for j=1:length(incAngle)
for i=1:length(omega)
% for j=1:length(incAngle)
%% initial data
layerstack.NumPts = 64;
layerstack.kappa = omega(i)*sin(pi.*incAngle(j)./180);
% layerstack.omega = 2*pi/500; 
layerstack.omega = omega(i);
layerstack.period = 500;
layerstack.kztolfactor = 1e-6;

% layerstack.layer_thicknesses = [500, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell ,hmm_unit_cell ,hmm_unit_cell ,hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, hmm_unit_cell, 500 ];
layerstack.layer_thicknesses  = [500 100 100 500]; 
% layerstack.layer_thicknesses = [100, 2,18, 2,18, 2,18, 2,18, 2,18, 2,18, 2,18, 2,18, 2,18, 2,18, 100 ];
% layerstack.layer_thicknesses = [100, 0,20, 0,20, 0,20, 0,20, 0,20, 0,20, 0,20, 0,20, 0,20, 0,20, 100 ];

xvec_sample = 0:5:1500;
zvec_sample = 0:3:1100;

%% input field functions

layerstack.UinfuncEy = @(x) OneFunc(x).*exp(1i.*layerstack.kappa.*x);
layerstack.UinfuncHy = @(x) OneFunc(x).*exp(1i.*layerstack.kappa.*x);
layerstack.DinfuncEy = @(x) ZeroFunc(x);
layerstack.DinfuncHy = @(x) ZeroFunc(x);

%% construction of material tensors
%helping fuctions
% specify needed epsilon values (air, diamond, gold)
% eps1 = [1.0,epsilonSiO2,epsAg(layerstack.omega)];
eps1 = [1.0,2.25];
eps1func = @(x) eps1(1).*OneFunc(x);
eps2func = @(x) eps1(1).*OneFunc(x) + (eps1(2) -eps1(1)).*(BoxFunc(mod(x,layerstack.period),0*layerstack.period,covarage_factor*layerstack.period));
eps3func = @(x) eps1(1).*OneFunc(x) + (eps1(3) -eps1(1)).*(BoxFunc(mod(x,layerstack.period),0*layerstack.period,covarage_factor*layerstack.period));
eps4func = @(x) eps1(2).*OneFunc(x);
mufunc = @(x) OneFunc(x);

layerstack.epsstruct(1).xx = eps1func;
layerstack.epsstruct(1).yy = eps1func;
layerstack.epsstruct(1).zz = eps1func;

layerstack.epsstruct(2).xx = eps2func;
layerstack.epsstruct(2).yy = eps2func;
layerstack.epsstruct(2).zz = eps2func;

layerstack.epsstruct(3).xx = eps4func;
layerstack.epsstruct(3).yy = eps4func;
layerstack.epsstruct(3).zz = eps4func;

layerstack.epsstruct(4).xx = eps1func;
layerstack.epsstruct(4).yy = eps1func;
layerstack.epsstruct(4).zz = eps1func;

% layerstack.epsstruct(4) = layerstack.epsstruct(2);
% layerstack.epsstruct(5) = layerstack.epsstruct(3);
% 
% layerstack.epsstruct(6) = layerstack.epsstruct(4);
% layerstack.epsstruct(7) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(8) = layerstack.epsstruct(4);
% layerstack.epsstruct(9) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(10) = layerstack.epsstruct(4);
% layerstack.epsstruct(11) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(12) = layerstack.epsstruct(2);
% layerstack.epsstruct(13) = layerstack.epsstruct(3);
% 
% layerstack.epsstruct(14) = layerstack.epsstruct(4);
% layerstack.epsstruct(15) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(16) = layerstack.epsstruct(4);
% layerstack.epsstruct(17) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(18) = layerstack.epsstruct(4);
% layerstack.epsstruct(19) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(20) = layerstack.epsstruct(4);
% layerstack.epsstruct(21) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(22) = layerstack.epsstruct(4);
% layerstack.epsstruct(23) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(24) = layerstack.epsstruct(4);
% layerstack.epsstruct(25) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(26) = layerstack.epsstruct(4);
% layerstack.epsstruct(27) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(28) = layerstack.epsstruct(4);
% layerstack.epsstruct(29) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(30) = layerstack.epsstruct(4);
% layerstack.epsstruct(31) = layerstack.epsstruct(5);
% 
% layerstack.epsstruct(32) = layerstack.epsstruct(1);


% layerstack.epsstruct(4) = layerstack.epsstruct(3);

layerstack.mustruct(1).xx = mufunc;
layerstack.mustruct(1).yy = mufunc;
layerstack.mustruct(1).zz = mufunc;

layerstack.mustruct(2) = layerstack.mustruct(1);
layerstack.mustruct(3) = layerstack.mustruct(1);
layerstack.mustruct(4) = layerstack.mustruct(1);
% layerstack.mustruct(5) = layerstack.mustruct(1);
% layerstack.mustruct(6) = layerstack.mustruct(1);
% layerstack.mustruct(7) = layerstack.mustruct(1);
% layerstack.mustruct(8) = layerstack.mustruct(1);
% layerstack.mustruct(9) = layerstack.mustruct(1);
% layerstack.mustruct(10) = layerstack.mustruct(1);
% layerstack.mustruct(11) = layerstack.mustruct(1);
% layerstack.mustruct(12) = layerstack.mustruct(1);
% layerstack.mustruct(13) = layerstack.mustruct(1);
% layerstack.mustruct(14) = layerstack.mustruct(1);
% layerstack.mustruct(15) = layerstack.mustruct(1);
% layerstack.mustruct(16) = layerstack.mustruct(1);
% layerstack.mustruct(17) = layerstack.mustruct(1);
% layerstack.mustruct(18) = layerstack.mustruct(1);
% layerstack.mustruct(19) = layerstack.mustruct(1);
% layerstack.mustruct(20) = layerstack.mustruct(1);
% layerstack.mustruct(21) = layerstack.mustruct(1);
% layerstack.mustruct(22) = layerstack.mustruct(1);
% layerstack.mustruct(23) = layerstack.mustruct(1);
% layerstack.mustruct(24) = layerstack.mustruct(1);
% layerstack.mustruct(25) = layerstack.mustruct(1);
% layerstack.mustruct(26) = layerstack.mustruct(1);
% layerstack.mustruct(27) = layerstack.mustruct(1);
% layerstack.mustruct(28) = layerstack.mustruct(1);
% layerstack.mustruct(29) = layerstack.mustruct(1);
% layerstack.mustruct(30) = layerstack.mustruct(1);
% layerstack.mustruct(31) = layerstack.mustruct(1);
% layerstack.mustruct(32) = layerstack.mustruct(1);

% layerstack.mustruct(4) = layerstack.mustruct(1);

%% calculate power contributions by smatrix algorithmn
tstartidx = tic;
 [ P0incEy,P0scatEy,PendincEy,PendscatEy,P0incHy,P0scatHy,PendincHy,PendscatHy, Ex,Ey,Ez,Hx,Hy,Hz,diffraction_order_data] = smat_with_modes( layerstack,xvec_sample,zvec_sample);
%  [ P0incEy,P0scatEy,PendincEy,PendscatEy,P0incHy,P0scatHy,PendincHy,PendscatHy,diffraction_order_data] = smat( layerstack);
time_in_smat(i,j) = toc(tstartidx);
toc(tstartidx);

TEy(i,j) = PendscatEy/P0incEy;
REy(i,j) = - P0scatEy/P0incEy;
THy(i,j) = PendscatHy/P0incHy;
RHy(i,j) = - P0scatHy/P0incHy;
% A(i) = 1 - T(i) -R(i);

%% diffraction order related spectra

% TEy_by_order = diffraction_order_data.PendscatEy./P0incEy;
% REy_by_order = - diffraction_order_data.P0scatEy./P0incEy;
% THy_by_order = diffraction_order_data.PendscatHy./P0incHy;
% RHy_by_order = - diffraction_order_data.P0scatHy./P0incHy;
% 
% figure(5);
% plot(diffraction_order_data.outAngleend,abs(TEy_by_order),'.-b',diffraction_order_data.outAngle0,abs(REy_by_order),'.-r');
% drawnow;
% figure(6);
% plot(diffraction_order_data.outAngleend,abs(THy_by_order),'.-b',diffraction_order_data.outAngle0,abs(RHy_by_order),'.-r');
% drawnow;

%

tmp_idx = find(diffraction_order_data.order_number == 0,1);

TEy_Oth_order(i,j) = diffraction_order_data.PendscatEy(tmp_idx)./P0incEy;
REy_Oth_order(i,j) = - diffraction_order_data.P0scatEy(tmp_idx)./P0incEy;
THy_Oth_order(i,j) = diffraction_order_data.PendscatHy(tmp_idx)./P0incHy;
RHy_Oth_order(i,j) = - diffraction_order_data.P0scatHy(tmp_idx)./P0incHy;

clear tmp_idx;

%%
figure(3);
imagesc(xvec_sample,zvec_sample,real(Ey).');
% title('|E_s|','fontsize',20);
% xlabel('position x (\mum)','fontsize',20);
% ylabel('position z (\mum)','fontsize',20);
colorbar;
% caxis([0,2]);

% set(gca,'fontsize',20,'linewidth',2);
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 20*3/4])
% set(gcf,'papersize',[30,20]);

drawnow;
% print('-dpng',['./pictures/bs_scan_period__2000__thickness__400__N__512__spol_mode_10mu.png'],'-r100');


figure(4);
% imagesc(xvec_sample,zvec_sample,(sqrt(real(Ex).^2+real(Ez).^2)).');
imagesc(xvec_sample,zvec_sample,real(Hy).');
% title('|E_p|','fontsize',20);
% xlabel('position x (\mum)','fontsize',20);
% ylabel('position z (\mum)','fontsize',20);
colorbar;
caxis([0,2]);
% set(gca,'fontsize',20,'linewidth',2);
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 20*3/4])
% set(gcf,'papersize',[30,20]);

drawnow;
% print('-dpng',['./pictures/bs_scan_period__2000__thickness__400__N__512__ppol_mode_10mu.png'],'-r100');
% 
% figure(5);
% imagesc(log10(abs(fftshift(fft2(Ey.'))))); colorbar;
% 
% drawnow;
% 
% figure(6);
% imagesc(log10(abs(fftshift(fft2(Hy.'))))); colorbar;
% drawnow;

end;
end;
%% post
% figure(1);
% plot(wavelength,T,'.-',wavelength,R,'.-',wavelength,A,'.-');
% figure(2);
% plot(omega,T,'.-',omega,R,'.-',omega,A,'.-');

figure(1);
plot(wavelength,TEy,'.-b',wavelength,REy,'.-r',wavelength,TEy_Oth_order,'x-b',wavelength,REy_Oth_order,'x-r');

figure(2);
plot(wavelength,THy,'.-b',wavelength,RHy,'.-r',wavelength,THy_Oth_order,'x-b',wavelength,RHy_Oth_order,'x-r');

figure(3);
plot(incAngle,TEy.','.-b',incAngle,REy.','.-r',incAngle,TEy_Oth_order.','o-b',incAngle,REy_Oth_order.','o-r');

figure(4);
plot(incAngle,THy.','.-b',incAngle,RHy.','.-r',incAngle,THy_Oth_order.','o-b',incAngle,RHy_Oth_order.','o-r');

figure(5);
subplot(1,3,1)
imagesc(incAngle,wavelength,TEy)
colorbar;
caxis([0,1]);
title('TEy')

subplot(1,3,2)
imagesc(incAngle,wavelength,REy)
colorbar;
caxis([0,1]);
title('REy')

subplot(1,3,3)
imagesc(incAngle,wavelength,1-TEy-REy)
colorbar;
caxis([0,1]);
title('AEy')

figure(6);
subplot(1,3,1)
imagesc(incAngle,wavelength,THy)
colorbar;
caxis([0,1]);
title('THy')

subplot(1,3,2)
imagesc(incAngle,wavelength,RHy)
colorbar;
caxis([0,1]);
title('RHy')

subplot(1,3,3)
imagesc(incAngle,wavelength,1-THy-RHy)
colorbar;
caxis([0,1]);
title('AHy')

% %%
% figure(3);
% imagesc(xvec_sample,zvec_sample,real(Ey).');
% colorbar;
% 
% figure(4);
% imagesc(xvec_sample,zvec_sample,real(Hy).');
% colorbar;

%
%% calculate layer data
% 
% [ layer_data(1) ] = process_single_layer( layerstack, 1 );
% [ layer_data(2) ] = process_single_layer( layerstack, 2 );
% [ layer_data(3) ] = process_single_layer( layerstack, 3 );

%%
% figure(13);
% plot(1:layerstack.NumPts,real(layer_data(1).kzvalsEy),'.',1:layerstack.NumPts,imag(layer_data(1).kzvalsEy),'.');

%% post processing
% x = -500:1:1500;
% x = layer_data(1).xvec;
% [ Ex,Ey,Ez,Hx,Hy,Hz ] = get_layer_modes( layer_data(1) , x);
% %%
% plot(x,real(Ey(:,1:7)),'.-')
