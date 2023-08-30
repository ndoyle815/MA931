% first attempt - SEIR deterministic with a simple "social distancing"
% control
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

load mats/Distributions.mat
load mats/Probabilities.mat

% model parameters
gamma = 1/7;                     % infectious period
sigma = 1/5.28;                  % latency period
omega = 0*1/200;                 % recovered period
tau = 0.25;                      % relative infectiousness of asymptomatic
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
N = 2.*[50000; 125000; 40000];   % population structure
n = size(N,1);                   % number of age classes
atrisk_prop = N(end)/sum(N);     % at-risk proportion of population

% transmission matrix
%load('./mats/contacts.mat')
%beta = 1.0.*[1.35, 2.41, 0.3; 0.93, 3.01, 0.54; 0.35, 1.64, 1.21].*gamma;
beta = 0.7.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 500;  % main simulation

% define strategy numbers and thresholds
strategies = [1:6];
%thresholds = [1500 6000; 1000 2000; 2000 4000; 4000 6000; 1000 8000; 4000 8000];
thresholds = [1000 6000; 500 1500; 1000 3000; 2500 4500; 500 7000; 3000 7000];

% arrays to store lockdown metrics
[Days_lockdown, Peak_incidence, Peak_hospital, Last_lockdown, FinalSize, FinalHospital, FinalDeaths, NPhases] = deal([]);

figure('Position',[400 400 1000 1000])

tic
for strat = strategies
    % Define model parameters as a structure
    para = struct('beta',beta,'gamma',gamma,'sigma',sigma,'omega',omega,'tau',tau, ...
                  'da',da,'N',N,'n',n,'strategy',strat,'init',0,'maxtime',t_init, ...
                  'tgap',16,'tdelay',5);

    % dummy thresholds to allow infections to build with no intervention
    para.Imin = 20000;
    para.Imax = 20000;
    para.Imin_risk = para.Imin*atrisk_prop;
    para.Imax_risk = para.Imax*atrisk_prop;

    % Define initial conditions as a structure
    E0 = 0.005;  % initial exposed
    ICs = struct('S',(1-E0).*para.N, 'E1',E0.*para.N, 'E2',zeros(n,1), 'E3',zeros(n,1), ...
                 'IS',zeros(n,1), 'IA',zeros(n,1), 'R',zeros(n,1), 'Cases',zeros(n,1));

    % Preliminary run - allows us to get ICs to begin in lockdown
    %[Prelim] = SEIR_demo_NBD(para,ICs,30);
    [Prelim] = SEIR_demo(para,ICs);
    
    % add control thresholds defined by strategy
    para.init = 1;
    para.maxtime = maxtime;
    para.Imin = thresholds(strat,1);
    para.Imax = thresholds(strat,2);
    para.Imin_risk = para.Imin*atrisk_prop;
    para.Imax_risk = para.Imax*atrisk_prop;

    % Get new initial conditions to begin in lockdown
    ICs = struct('S',Prelim.S(end,:), 'E1',Prelim.E1(end,:), 'E2',Prelim.E2(end,:), ...
                 'E3',Prelim.E3(end,:), 'IS',Prelim.IS(end,:), 'IA',Prelim.IA(end,:), ...
                 'R',Prelim.R(end,:),'Cases',Prelim.Cases(end,:));

    % Run model
    %[Classes] = SEIR_demo_discretised(para,ICs);
    [Classes] = SEIR_demo(para,ICs);

    % Find times where social distancing is enforced
    nx = size(Classes.SD,1);

    % Compute daily new infections, hospitalisations and deaths
    hosp_rates = [0.1; 0.15; 0.3];
    death_rates = [0.001; 0.01; 0.1];
    %hosp_lag = 5;
    
    Classes.DailyInf = [0 0 0; Prelim.Cases(2:end,:) - Prelim.Cases(1:end-1,:); Classes.Cases(2:end,:) - Classes.Cases(1:end-1,:)];
    HH = Classes.DailyInf*hosp_rates;
    newHH = zeros(maxtime+t_init+1,1);
    lag = [0:length(Distribution_Symptoms_to_Hospital)-1];
    for t = t_init+1:length(newHH)
        newHH(t) = newHH(t) + Distribution_Symptoms_to_Hospital*HH(t-lag);
    end
    Classes.Hospital = newHH(t_init+1:end);
    Classes.Deaths = Classes.DailyInf*death_rates;

    DL = 0;

    % Plot total infections
    %figure('Position',[400 400 600 400])
    %clf
    subplot(3,2,strat)
    for i = 1:4:nx
        DL = DL + Classes.SD(i+2,1) - Classes.SD(i,1);
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    plot(Classes.t, sum(Classes.IS,2), 'k', 'LineWidth', 2.5)
    %hold on
    %plot(Classes.t, 10.*Classes.Hospital, 'b', 'LineWidth', 2.5)
    %plot(Classes.t, Classes.IS(:,3)./sum(Classes.IS,2).*10000, 'b', 'LineWidth', 2.5)
    %hold on
    %plot(Classes.t, Classes.IA(:,3)./sum(Classes.IA,2).*10000, 'g', 'LineWidth', 2.5)
    hold on
    yline(para.Imax,':','Upper','LineWidth',1.5)
    hold on
    yline(para.Imin,':','Lower','LineWidth',1.5)
    if strat > 4
        xlabel('Time (days)')
    end
    if mod(strat,2) == 1
        ylabel('Infected Individuals')
    end
    title(strcat("Strategy S",num2str(strat)))
    axis([0 maxtime 0 max([10000,1.1*max(sum(Classes.IS,2))])])
    grid on

    % Compute metrics
    Peak_incidence = [Peak_incidence round(max(sum(Classes.IS,2)))];
    Peak_hospital = [Peak_hospital round(max(Classes.Hospital))];
    FinalSize = [FinalSize round(sum(Classes.R(end,:)))];
    FinalHospital = [FinalHospital round(sum(Classes.Hospital))];
    FinalDeaths = [FinalDeaths round(Classes.Cases(end,:)*death_rates)];
    %Days_over_x
    %Days_over_y
    Last_lockdown = [Last_lockdown round(Classes.SD(end,1))];
    Days_lockdown = [Days_lockdown round(DL)];
    NPhases = [NPhases ceil(nx/2)];
    
    %save figure
    %saveas(gcf,strcat('./images/Strat_',num2str(strat),'.png'))
end
toc

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