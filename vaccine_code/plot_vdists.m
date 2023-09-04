% script to plot vaccine arrival date distributions and generate Figure 2
% from the report

clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

vtimes = [180:60:1080];

% generate distributions using discretenormal.m
vdist_pois1 = discretenormal([0:length(vtimes)-1], 3,   1.5);
vdist_pois2 = discretenormal([0:length(vtimes)-1], 7.5, 1.5);
vdist_pois3 = discretenormal([0:length(vtimes)-1], 12,  1.5);
vdist_pois4 = discretenormal([0:length(vtimes)-1], 7.5, 1);
vdist_pois5 = discretenormal([0:length(vtimes)-1], 7.5, 1.5);
vdist_pois6 = discretenormal([0:length(vtimes)-1], 7.5, 2);


f = figure(1);
f.Position = [600 600 600 600];
%sgtitle('Vaccine arrival date distributions','FontSize',18)

% mu varying, sigma fixed
subplot(2,1,1)
bar(vtimes,[vdist_pois1; vdist_pois2; vdist_pois3],1,"grouped",'FaceAlpha',0.8)
xticks([180:120:1080])
xtickangle(0)
axis([120 1140 0 0.3])
%xlabel('Vaccination arrival time (days)')
ylabel('Probability')
title('Variance fixed, $\sigma=1.5$')
legend('$\mu=3$','$\mu=7.5$','$\mu=12$','Interpreter','latex','Location','eastoutside')
grid on

% mu fixed, varying sigma
subplot(2,1,2)
bar(vtimes,[vdist_pois4; vdist_pois5; vdist_pois6],1,"grouped",'FaceAlpha',0.8)
xticks([180:120:1080])
xtickangle(0)
axis([120 1140 0 0.45])
xlabel('Vaccination arrival time (days)')
ylabel('Probability')
title('Mean fixed, $\mu=7.5$')
legend('$\sigma=1$','$\sigma=1.5$','$\sigma=2$','Interpreter','latex','Location','eastoutside')
grid on

saveas(gcf,'./vacc_images/vaccdists.png')
