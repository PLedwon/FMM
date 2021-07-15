numpts(1)
numpts(5)
numpts(8)
numpts(10)
k=12;
set(0, 'DefaultFigureRenderer', 'painters')

figure(1)
subplot(1,2,1)
plot(100*xvec,exMat(:,1),'linewidth',2),xlabel({'Position $x$ [nm]'},'Interpreter','latex','fontsize',k)
xlim([0, 400]), ylim([0.5, 1.1])
ylabel({'$E_x$'},'Interpreter','latex','fontsize',k);
legend({'$N=16$'},'Interpreter','latex','fontsize',k,'Location','south');



subplot(1,2,2)
plot(100*xvec,exMat(:,5),'linewidth',2),xlabel({'Position $x$ [nm]'},'Interpreter','latex','fontsize',k)
xlim([0, 400]), ylim([0.5, 1.1])

legend({'$N=64$'},'Interpreter','latex','fontsize',k,'Location','northwest');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])
print('N_low','-dpdf')

%%
figure(2)
subplot(1,2,1)
plot(100*xvec,exMat(:,8),'linewidth',2),xlabel({'Position $x$ [nm]'},'Interpreter','latex','fontsize',k)
xlim([0, 400]), ylim([0.5, 1.25])
ylabel({'$E_x$'},'Interpreter','latex','fontsize',k);
legend({'$N=256$'},'Interpreter','latex','fontsize',k,'Location','south');


subplot(1,2,2)
plot(100*xvec,exMat(:,10),'linewidth',2),xlabel({'Position $x$ [nm]'},'Interpreter','latex','fontsize',k)
xlim([0, 400]), ylim([0.5, 1.25])

legend({'$N=512$'},'Interpreter','latex','fontsize',k,'Location','south');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[6, 3],'Position',[0 0 6 3])
print('N_high','-dpdf')
