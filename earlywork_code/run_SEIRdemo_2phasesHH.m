% first attempt - SEIR deterministic with a simple "social distancing"
% control
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

load ../mats/Distributions.mat
load ../mats/Probabilities.mat

% model parameters
gamma = 1/7;                     % infectious period
sigma = 1/5;                     % latency period
omega = 0*1/400;                 % recovered period
tau = 0.25;                      % relative infectiousness of asymptomatic
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
N = 2.*[50000; 125000; 40000];   % population structure
n = size(N,1);                   % number of age classes
atrisk_prop = N(end)/sum(N);     % at-risk proportion of population
bedsper1000 = 2.5;               % UK hospital beds per 1,000 population
Hmax = bedsper1000*sum(N)/1000;  % capacity

% transmission matrix
beta = 0.7.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 600;  % main simulation

% define strategy numbers and thresholds
strategies = [1:6];
thresholds = [50 200 125 700; 50 150 100 200; 100 300 400 550; 150 250 350 450; 50 300 400 500; 300 600 400 700];

% arrays to store lockdown metrics
[Days_lockdown, Peak_incidence, Peak_hospital, Last_lockdown, FinalSize, FinalHospital, FinalDeaths, NPhases] = deal([]);

figure('Position',[400 400 1000 1000])

tic
for strat = strategies
    % Define model parameters as a structure
    para0 = struct('beta',beta,'gamma',gamma,'sigma',sigma,'omega',omega,'tau',tau, ...
                  'da',da,'N',N,'n',n,'strategy',strat,'init',0,'maxtime',t_init, ...
                  'tgap',18,'tdelay',3,'tdiff',7,'hosp_rates',[0.1; 0.15; 0.3], ...
                  'epsilon',1/8,'delta',1/10,'rho',0.1);

    % dummy thresholds to allow infections to build with no intervention
    para0.U12 = 20000;
    para0.U01 = 20000;
    para0.L10 = 20000;
    para0.L21 = 20000;    

    % Define initial conditions as a structure
    E0 = 6e-4;  % initial exposed
    ICs = struct('S',(1-E0).*para0.N, 'E1',E0.*para0.N, 'E2',zeros(n,1), 'E3',zeros(n,1), ...
                 'IA',zeros(n,1), 'IS',zeros(n,1), 'IPH',zeros(n,1), 'IH',zeros(n,1), ...
                 'R',zeros(n,1), 'Cases',zeros(n,1), 'Hosp',zeros(n,1));

    % run preliminary simulation
    [Prelim] = SEIR_demo_2phasesHH(para0,ICs);

    % add control thresholds defined by strategy
    para = para0;
    para.maxtime = maxtime;
    para.L10 = thresholds(strat,1);
    para.U01 = thresholds(strat,2);
    para.L21 = thresholds(strat,3);
    para.U12 = thresholds(strat,4);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.U12
        para.init = 1;
    else
        para.init = 2;
    end

    % Get new initial conditions to begin in lockdown
    ICs = struct('S',Prelim.S(end,:), 'E1',Prelim.E1(end,:), 'E2',Prelim.E2(end,:), ...
                 'E3',Prelim.E3(end,:), 'IA',Prelim.IA(end,:), 'IS',Prelim.IS(end,:), ...
                 'IPH',Prelim.IPH(end,:), 'IH',Prelim.IH(end,:), 'R',Prelim.R(end,:), ...
                 'Cases',Prelim.Cases(end,:), 'Hosp',Prelim.Hosp(end,:));

    % Run model
    [Classes] = SEIR_demo_2phasesHH(para,ICs);

    % Find times where social distancing is enforced
    nx = size(Classes.SD,1);
    ix1 = find(Classes.SD(:,2)==1);
    ix2 = find(Classes.SD(:,2)==2);

    % Compute daily new infections, hospitalisations and deaths
%     hosp_rates = [0.1; 0.15; 0.3];
%     death_rates = [0.001; 0.01; 0.1];
%     
%     Classes.DailyInf = [0 0 0; Prelim.Cases(2:end,:) - Prelim.Cases(1:end-1,:); Classes.Cases(2:end,:) - Classes.Cases(1:end-1,:)];
%     HH = Classes.DailyInf*hosp_rates;
%     newHH = zeros(maxtime+t_init+1,1);
%     lag = [0:length(Distribution_Symptoms_to_Hospital)-1];
%     for t = t_init+1:length(newHH)
%         newHH(t) = newHH(t) + Distribution_Symptoms_to_Hospital*HH(t-lag);
%     end
%     Classes.Hospital = newHH(t_init+1:end);
%     Classes.Deaths = Classes.DailyInf*death_rates;
% 
%     DL = sum(Classes.SD(ix1+2)-Classes.SD(ix1)) + sum(Classes.SD(ix2+2)-Classes.SD(ix2));
    
    % append SD for plotting if we end in a restriction
    if Classes.SD(end,1) ~= 0
        Classes.SD(end+2,:) = [para.maxtime 0];
    end

    % plotting
    subplot(3,2,strat)
    for i = ix1'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'y', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    for i = ix2'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    plot(Classes.t, sum(Classes.IH,2), 'k', 'LineWidth', 2.5)
    %hold on
    %plot(Classes.t, sum(Classes.IS,2), 'b', 'LineWidth', 2.5)
    yline(para.U12,':','U12','LabelVerticalAlignment','middle','LineWidth',1.5,'FontSize',8)
    yline(para.U01,':','U01','LabelVerticalAlignment','middle','LineWidth',1.5,'FontSize',8)
    yline(para.L21,':','L21','LabelVerticalAlignment','middle','LineWidth',1.5,'FontSize',8)
    yline(para.L10,':','L10','LabelVerticalAlignment','middle','LineWidth',1.5,'FontSize',8)

    if strat > 4
        xlabel('Time (days)')
    end
    if mod(strat,2) == 1
        ylabel('Number in Hospital')
    end
    axis([0 maxtime 0 max([Hmax,1.15*max(sum(Classes.IH,2))])])
    title(strcat("Strategy S",num2str(strat)))
    grid on

    %round(max(sum(Classes.IH,2)))
    % Compute metrics
%     Peak_incidence = [Peak_incidence round(max(sum(Classes.IS,2)))];
%     Peak_hospital = [Peak_hospital round(max(Classes.Hospital))];
%     FinalSize = [FinalSize round(sum(Classes.R(end,:)))];
%     FinalHospital = [FinalHospital round(sum(Classes.Hospital))];
%     FinalDeaths = [FinalDeaths round(Classes.Cases(end,:)*death_rates)];
%     Last_lockdown = [Last_lockdown round(Classes.SD(end,1))];
%     Days_lockdown = [Days_lockdown round(DL)];
%     NPhases = [NPhases ceil(nx/2)];
    
end
toc

%save figure
% saveas(gcf,strcat('../images/waningimmunity_t600.png'))
% 
% % save metrics to table
% save("../mats/metrics2phases.mat","Days_lockdown","Last_lockdown","NPhases","Peak_incidence","Peak_hospital", ...
%     "FinalSize","FinalHospital","FinalDeaths","strategies")