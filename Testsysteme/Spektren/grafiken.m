k=12;

figure(1)
subplot(1,2,1)
plot(EV_TE_M3,'x','color','b')
xlabel('Realteil','Interpreter','LaTex','Fontsize',k)
ylabel('Imagin\"arteil','Interpreter','LaTex','Fontsize',k)
legend('EW TE-Pol.','Interpreter','LaTex')
axes('position',[.2 .25 .15 .15])
box on % put box around new pair of axes
indexOfInterest = (abs(real(EV_TE_M3)) < 30 & abs(real(EV_TE_M3)) < 80);
xlabel(' ','Interpreter','LaTex')
ylabel(' ','Interpreter','LaTex')
plot(EV_TE_M3(indexOfInterest),'x')
axis tight

subplot(1,2,2)
plot(EV_TM_M3,'x','color','b')
xlabel('Realteil','Interpreter','LaTex','Fontsize',k)
legend('EW TM-Pol.','Interpreter','LaTex')

axes('position',[.7 .25 .15 .15])
box on % put box around new pair of axes
indexOfInterest = (abs(real(EV_TM_M3)) < 10000 & abs(real(EV_TM_M3)) < 10000);
xlabel(' ','Interpreter','LaTex')
ylabel(' ','Interpreter','LaTex')
plot(EV_TM_M3(indexOfInterest),'x')
axis tight

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])

print('M3','-dpdf')
%%

figure(2)
subplot(1,2,1)
plot(EV_TE_M2,'x','color','b')
xlabel('Realteil','Interpreter','LaTex','Fontsize',k)
ylabel('Imagin\"arteil','Interpreter','LaTex','Fontsize',k)
legend('EW TE-Pol.','Interpreter','LaTex')
axes('position',[.2 .25 .15 .15])
box on % put box around new pair of axes
indexOfInterest = (abs(real(EV_TE_M2)) < 200 & abs(real(EV_TE_M2)) < 200);
xlabel(' ','Interpreter','LaTex')
ylabel(' ','Interpreter','LaTex')
plot(EV_TE_M2(indexOfInterest),'x')
axis tight


subplot(1,2,2)
plot(EV_TM_M2,'x','color','b')
xlabel('Realteil','Interpreter','LaTex','Fontsize',k)
legend('EW TE-Pol.','Interpreter','LaTex')

axes('position',[.7 .25 .15 .15])
box on % put box around new pair of axes
indexOfInterest = (abs(real(EV_TM_M2)) < 200 & abs(real(EV_TM_M2)) < 200);
xlabel(' ','Interpreter','LaTex')
ylabel(' ','Interpreter','LaTex')
plot(EV_TM_M2(indexOfInterest),'x')
axis tight

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])

print('M2','-dpdf')

%%
figure(3)
subplot(1,2,1)
plot(EV_TE_M1,'x','color','b')
xlabel('Realteil','Interpreter','LaTex','Fontsize',k)
ylabel('Imagin\"arteil','Interpreter','LaTex','Fontsize',k)
legend('EW TE-Pol.','Interpreter','LaTex')

subplot(1,2,2)
plot(EV_TM_M1,'x','color','b')
xlabel('Realteil','Interpreter','LaTex','Fontsize',k)
legend('EW TM-Pol.','Interpreter','LaTex')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])

print('M1','-dpdf')