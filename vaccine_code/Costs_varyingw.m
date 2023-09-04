% script to run control strategies without vaccination and rank according
% to the objective function, producing subplots for Figure 4 in the report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% load default parameters
para0 = load('./mats/Parameters.mat');

% vaccination start times
vstart_times = [180:60:1080];
vstarts = [2*max(vstart_times), 360];

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 1000;  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 150 100 700; 50 150 100 200; 150 350 500 650; 275 350 425 500];
strategies = [1:length(thresholds)];

% add control thresholds defined by strategy
para = para0;
para.init = 1;
para.maxtime = maxtime;
para.Hmax = 1075;        % modify hospital capacity

% define functional weights
weights = [0.6:0.01:1];  % varying w1, w2 = 1 - w1
w3 = 2;

% stores cost function outputs
ns = length(strategies);
nw = length(weights);
fs = zeros(ns,nw);

tic
for strat = strategies
    % set switching thresholds
    para.T10 = thresholds(strat,1);
    para.T01 = thresholds(strat,2);
    para.T21 = thresholds(strat,3);
    para.T12 = thresholds(strat,4);

    % run preliminary simulation to get ICs
    [Prelim, Prelim_ICs] = Get_ICs(para0);

    % Run main simulation
    [Classes] = ODEmodel(para, Prelim_ICs);

    % use post-processor to compute metrics of interest
    [~, ~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, ~, ~, ~, ~, ~] = PostProcessor(Classes);
    
    for w = 1:nw
        % evaluate cost function
        fs(strat,w) = CostFunction([weights(w), 1-weights(w), w3], para, Peak_hospital, FinalHospital, Days_lockdown, Days_Tier2);
    end
end
toc

% rank strategies
[~,idx] = min(fs);
[sorted, ii] = sort(fs);
[~,rank] = sort(ii);

% Plotting
cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

figure('Position',[600 600 1200 300])
hold on
for strat = strategies
    plot(weights,rank(strat,:),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
end
hold off
set(gca, 'YDir','reverse')
axis([min(weights) max(weights) 0.5 4.5])
xlabel('$w_1 \; (w_2 = 1 - w_1)$')
ylabel('Strategy Ranking')
legend({'S1 (Cautious Easing)','S2 (Suppression)','S3 (Slow Control)','S4 (Rapid Control)'},'Interpreter','Latex','Location','eastoutside','FontSize',16)
title(strcat('$H_c$ = ',{' '},num2str(para.Hmax),', ',' maxtime =  ',{' '},num2str(maxtime)))
grid on

saveas(gcf,strcat('./vacc_images/CostFunction_',num2str(para.Hmax),'_',num2str(maxtime),'.png'))

save(strcat('./mats/StratRankings_',num2str(para.Hmax),'_',num2str(maxtime),'.mat'),"rank")
