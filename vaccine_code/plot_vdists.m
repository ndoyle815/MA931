clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

vtimes = [180:60:1080];
%vdist_pois1 = poisspdf([0:length(vtimes)-1],4);  % poisson distribution
%vdist_pois2 = poisspdf([0:length(vtimes)-1],8);  % poisson distribution
%vdist_pois3 = poisspdf([0:length(vtimes)-1],12);  % poisson distribution
vdist_pois1 = discretenormal([0:length(vtimes)-1],3,1.5);  % poisson distribution
vdist_pois2 = discretenormal([0:length(vtimes)-1],7.5,1.5);  % poisson distribution
vdist_pois3 = discretenormal([0:length(vtimes)-1],12,1.5);  % poisson distribution
vdist_pois4 = discretenormal([0:length(vtimes)-1],7.5,1);  % poisson distribution
vdist_pois5 = discretenormal([0:length(vtimes)-1],7.5,1.5);  % poisson distribution
vdist_pois6 = discretenormal([0:length(vtimes)-1],7.5,2);  % poisson distribution
vdist_unif = 1/length(vtimes).*ones(1,length(vtimes));  % uniform distribution


f = figure(1);
f.Position = [600 600 600 600];
%sgtitle('Vaccine arrival date distributions','FontSize',18)

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

%lambda = 1/575;
%x = [0:1200];

%figure(2)
%plot(x,1-exp(-lambda*x),'Color','b')
%hold on
%scatter(22*7,1-0.762,50,'k','filled')
%xline(1/lambda,'color','r')
%axis([min(x) max(x) 0 1.1])
%grid on