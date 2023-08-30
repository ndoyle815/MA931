clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',22)

A1 = load('StratRankings_1500_500.mat');
A2 = load('StratRankings_1500_1000.mat');
A3 = load('StratRankings_1500_1500.mat');
A4 = load('StratRankings_1075_1000.mat');

Ranks(:,:,1) = A1.rank;
Ranks(:,:,2) = A2.rank;
Ranks(:,:,3) = A3.rank;
Ranks(:,:,4) = A4.rank;

weights = [0.6:0.01:1];

cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

figure('Position',[600 600 800 1200])
%[1+2,4+5,7+8,10+11]
%[1:3,5:7,9:11,13:15]
for k = 1:4
    %subplot(4,3,[3*k-2,3*k-1])
    %subplot(4,4,[4*k-3,4*k-1])
    subplot(4,1,k)
    hold on
    for strat = 1:4
        plot(weights,Ranks(strat,:,k),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
    end
    hold off
    set(gca, 'YDir','reverse')
    axis([min(weights) max(weights) 0.5 4.5])
    yticks([1,2,3,4])
    ylabel('Ranking')
    
    if k == 4
        xlabel('$w_1 \; (w_2 = 1 - w_1)$')
        legend({'S1 (Cautious Easing)','S2 (Suppression)','S3 (Slow Control)','S4 (Rapid Control)'},'Fontsize',14,'Interpreter','Latex','Location','east');
    end
    %title(strcat('$H_c$ = ',{' '},num2str(Hc),', ',' maxtime =  ',{' '},num2str(maxtime)))
    grid on
end

% add a bit space to the figure
%fig = gcf;
%fig.Position(3) = fig.Position(3) + 100;
% add legend
%Lgnd = legend({'S1 (Cautious Easing)','S2 (Suppression)','S3 (Slow Control)','S4 (Rapid Control)'},'Interpreter','Latex','Location','eastoutside','FontSize',20);
%Lgnd.Position(1) = 0.72;
%Lgnd.Position(2) = 0.5;

saveas(gcf,'../Costoutputs.png')