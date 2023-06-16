% first attempt - SEIR deterministic with a simple "social distancing"
% control
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% model parameters
gamma = 1/7;                     % infectious period
sigma = 1/5;                     % latency period
tau = 0.25;                      % relative infectiousness of asymptomatic
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
N = 3.*[50000; 125000; 40000];   % population structure
n = size(N,1);                   % number of age classes
atrisk_prop = N(end)/sum(N);     % at-risk proportion of population

% transmission matrix
%load('./mats/contacts.mat')
beta = 1.0.*[1.35, 2.41, 0.3; 0.93, 3.01, 0.54; 0.35, 1.64, 1.21].*gamma;

% Define time to run model for
maxtime = 600;

% define strategy numbers and thresholds
strategies = [1:6];
thresholds = [1500 6000; 1000 2000; 2000 4000; 4000 6000; 1000 8000; 4000 8000];

% arrays to store lockdown metrics
[Days_lockdown, Peak_incidence, Peak_hospital, Last_lockdown, FinalSize, FinalHospital, FinalDeaths, NPhases] = deal([]);

figure('Position',[400 400 1000 1000])

for strat = 1:6
    % Define model parameters as a structure
    para = struct('beta',beta,'gamma',gamma,'sigma',sigma,'tau',tau,'da',da, ...
                'N',N,'n',n,'strategy',strat);

    % dummy thresholds to allow infections to build with no intervention
    para.Imin = 20000;
    para.Imax = 20000;
    para.Imin_risk = para.Imin*atrisk_prop;
    para.Imax_risk = para.Imax*atrisk_prop;

    % Define initial conditions as a structure
    E0 = 0.01;  % initial exposed
    ICs = struct('S',(1-E0).*para.N, 'E1',E0.*para.N, 'E2',zeros(n,1), 'E3',zeros(n,1), ...
                 'IS',zeros(n,1), 'IA',zeros(n,1), 'R',zeros(n,1), 'Cases',zeros(n,1));

    % Preliminary run - allows us to get ICs to begin in lockdown
    %[Prelim] = SEIR_demo_NBD(para,ICs,30);
    t_init = 30;
    [Prelim] = SEIR_demo(para,ICs,t_init);

    % add control thresholds defined by strategy
    para.Imin = thresholds(strat,1);
    para.Imax = thresholds(strat,2);
    para.Imin_risk = para.Imin*atrisk_prop;
    para.Imax_risk = para.Imax*atrisk_prop;

    % Get new initial conditions to begin in lockdown
    ICs = struct('S',Prelim.S(end,:), 'E1',Prelim.E1(end,:), 'E2',Prelim.E2(end,:), ...
                 'E3',Prelim.E3(end,:), 'IS',Prelim.IS(end,:), 'IA',Prelim.IA(end,:), ...
                 'R',Prelim.R(end,:),'Cases',Prelim.Cases(end,:));

    % Run model
    %[Classes] = SEIR_demo_NBD(para,ICs,maxtime);
    [Classes] = SEIR_demo(para,ICs,maxtime);

    % Find times where social distancing is enforced
    nx = size(Classes.SD,1);

    % Compute daily new infections, hospitalisations and deaths
    hosp_rates = [0.1; 0.15; 0.3];
    death_rates = [0.001; 0.01; 0.1];
    Classes.DailyInf = [Prelim.Cases(2:end,:) - Prelim.Cases(1:end-1,:); Classes.Cases(2:end,:) - Classes.Cases(1:end-1,:)];
    Classes.Hospital = Classes.DailyInf*hosp_rates;
    Classes.Deaths = Classes.DailyInf*death_rates;

    DL = 0;

    % Plot total infections
    %figure('Position',[400 400 600 400])
    %clf
    subplot(3,2,strat)
    for i = 1:2:nx
        DL = DL + Classes.SD(i+1,1) - Classes.SD(i,1);
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+1,1) Classes.SD(i+1,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    plot(Classes.t, sum(Classes.IS,2), 'k', 'LineWidth', 2.5)
    hold on
    yline(para.Imax,':','Upper Threshold')
    hold on
    yline(para.Imin,':','Lower Threshold')
    if strat > 4
        xlabel('Time (days)')
    end
    if mod(strat,2) == 1
        ylabel('Infected Individuals')
    end
    title(strcat("Strategy ",num2str(strat)))
    axis([0 maxtime 0 max([10000,1.15*max(sum(Classes.IS,2))])])
    grid on

    % Compute metrics
    Peak_incidence = [Peak_incidence round(max(sum(Classes.IS,2)))];
    Peak_hospital = [Peak_hospital round(max(Classes.Hospital))];
    FinalSize = [FinalSize round(sum(Classes.R(end,:)))];
    FinalHospital = [FinalHospital round(Classes.Cases(end,:)*hosp_rates)];
    FinalDeaths = [FinalDeaths round(Classes.Cases(end,:)*death_rates)];
    %Days_over_x
    %Days_over_y
    Last_lockdown = [Last_lockdown floor(Classes.SD(end,1))];
    Days_lockdown = [Days_lockdown round(DL)];
    NPhases = [NPhases nx];

    %save figure
    %saveas(gcf,strcat('./images/Strat_',num2str(strat),'.png'))
end

% I have been using this figure to check the total age cohort sizes remain
% constant
%figure(2)
%for j = 1:3
%    plot(Classes.t, Classes.S(:,j) + Classes.E1(:,j) + Classes.E1(:,j) + Classes.E1(:,j) + Classes.I(:,j) + Classes.R(:,j))
%    hold on
%end
%legend({'0-19','20-64','65+'})
%axis([0 maxtime para.N(3)-10 para.N(3)+10])
%grid on

%save figure
saveas(gcf,strcat('./images/Strats.png'))

%figure(2)
%plot(Prelim.t, sum(Prelim.I,2))
%grid on

% save metrics to table
save("./mats/metrics.mat","Days_lockdown","Last_lockdown","NPhases","Peak_incidence","Peak_hospital", ...
    "FinalSize","FinalHospital","FinalDeaths","strategies")