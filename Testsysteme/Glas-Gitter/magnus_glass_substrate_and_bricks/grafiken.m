M=load('magnus_glassystem.mat')
p = genpath('/users/stud/ledwon/Documents/Bachelorarbeit/Backups/Testsysteme/Glass-Gitter/glass_substrate_and_bricks/');
addpath(p)
S=load('glassystem.mat')
set(0, 'DefaultFigureRenderer', 'painters')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex')
k=12;
%%
figure(1)
subplot(1,2,1)

imagesc(100*M.wavelength,M.incAngle,log10(abs(M.MagnusTEy(:,:,5).'+M.MagnusREy(:,:,5).'-1)))

c=colorbar('northoutside','Direction','reverse');
title(c,'$\log _{10} |1-T_{TE}-R_{TE}| $ (Magnus)','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)

subplot(1,2,2)

imagesc(100*M.wavelength,M.incAngle,log10(abs(S.TEy(:,:,5).'+S.REy(:,:,5).'-1)))
c=colorbar('northoutside','Direction','reverse');
title(c,'$\log _{10} |1-T_{TE}-R_{TE}| $ ','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
%ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])


print('T_plus_R_TE','-dpdf')
%%
figure(8)
subplot(1,2,1)
imagesc(100*M.wavelength,M.incAngle,log10(abs(M.MagnusTHy(:,:,5).'+M.MagnusRHy(:,:,5).'-1)))
c=colorbar('northoutside','Direction','reverse');
title(c,'$\log _{10} |1-T_{TM}-R_{TM}| $ (Magnus)','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)

subplot(1,2,2)
imagesc(100*M.wavelength,M.incAngle,log10(abs(S.THy(:,:,5).'+S.RHy(:,:,5).'-1)))
c=colorbar('northoutside','Direction','reverse');
title(c,'$\log _{10} |1-T_{TM}-R_{TM}| $','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])
print('T_plus_R_TM','-dpdf')
%%
figure(3)
subplot(1,2,1)
imagesc(100*M.wavelength,M.incAngle,M.MagnusREy(:,:,5).')
c=colorbar('northoutside');
title(c,'$R_{TE}$ (Magnus)','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)

subplot(1,2,2)
imagesc(100*S.wavelength,S.incAngle,S.REy(:,:,5).')
c=colorbar('northoutside');
title(c,'$R_{TE}$','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])

print('R_TE','-dpdf')
%%
figure(4)
subplot(1,2,1)
imagesc(100*M.wavelength,M.incAngle,M.MagnusRHy(:,:,5).')
c=colorbar('northoutside');
title(c,'$R_{TM}$ (Magnus)','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)

subplot(1,2,2)
imagesc(100*S.wavelength,S.incAngle,S.RHy(:,:,5).')
c=colorbar('northoutside');
title(c,'$R_{TM}$','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])

print('R_TM','-dpdf')

%%
figure(7)
imagesc(100*M.wavelength,M.incAngle,log10(abs(M.MagnusTHy(:,:,5).'-S.THy(:,:,5).')))
c=colorbar('northoutside','Direction','reverse')
title(c,'$\log_{10}|\Delta R_{TM}|$','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])


print('R_TM_S_M','-dpdf')
%%
figure(10)
imagesc(100*M.wavelength,M.incAngle,log10(abs(M.MagnusTEy(:,:,5).'-S.TEy(:,:,5).')))
c=colorbar('northoutside','Direction','reverse')
title(c,'$\log_{10}|\Delta R_{TE}|$','Interpreter','LaTex','Fontsize', k)
c.FontSize = 10;

xlabel('Wellenl\"ange $\lambda$ [nm]','Interpreter','LaTex','Fontsize',k)
ylabel('Einfallswinkel $\varphi$ [$^\circ$]','Interpreter','LaTex','Fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])

print('R_TE_S_M','-dpdf')
%%
wavelength_idx = 49; % 49 enstspricht 560nm
incangle_idx = 15;
for i=1:length(M.numpts)
    konvTTM(i) = abs(M.MagnusTHy(wavelength_idx,incangle_idx,i)-M.MagnusTHy(wavelength_idx,incangle_idx,length(M.numpts)))./M.MagnusTHy(wavelength_idx,incangle_idx,length(M.numpts));
    konvRTM(i) = abs(M.MagnusRHy(wavelength_idx,incangle_idx,i)-M.MagnusRHy(wavelength_idx,incangle_idx,length(M.numpts)))./M.MagnusRHy(wavelength_idx,incangle_idx,length(M.numpts));
    konvTTE(i) = abs(M.MagnusTEy(wavelength_idx,incangle_idx,i)-M.MagnusTEy(wavelength_idx,incangle_idx,length(M.numpts)))./M.MagnusTEy(wavelength_idx,incangle_idx,length(M.numpts));
    konvRTE(i) = abs(M.MagnusREy(wavelength_idx,incangle_idx,i)-M.MagnusREy(wavelength_idx,incangle_idx,length(M.numpts)))./M.MagnusREy(wavelength_idx,incangle_idx,length(M.numpts));
   % konvM16TM(i) = abs(MNs16.THy(wavelength_idx,1,i)-MNs16.THy(wavelength_idx,1,length(MNs1.numpts)));
end
figure(11)
loglog(M.numpts,konvTTM,'x')
hold on
loglog(M.numpts,konvRTM,'o')
loglog(M.numpts,konvTTE,'d')
loglog(M.numpts,konvRTE,'.')
%loglog(MNs1.numpts,konvM16TM)
hold off
xlabel('Diskretisierungspunkte $N$','Interpreter','LaTex','FontSize',k)
ylabel('rel. Differenz $\Delta$','Interpreter','LaTex','Fontsize',k)
legend({'Transmission, TM', 'Reflexion, TM','Transmission, TE', 'Reflexion, TM'},'fontsize',k)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])
print('Konv_GG','-dpdf')


