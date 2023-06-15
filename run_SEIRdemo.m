% first attempt - SEIR deterministic with a simple "social distancing"
% control
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% model parameters
gamma = 1/7;
sigma = 1/3.5;
N = [50000; 125000; 40000];   % population structure
n = size(N,1);                % number of age classes
atrisk_prop = N(end)/sum(N);  % at-risk proportion of population

% transmission matrix
load('./mats/contacts.mat')
%beta = 0.5*gamma.*ones(n);
%beta(1:n-1,1:n-1) = 2*gamma;
beta = contacts;

% Define time to run model for
maxtime = 500;

% define strategy numbers and thresholds
strategies = [1:6];
thresholds = [1600 7500; 1200 2000; 2500 4500; 4000 7000; 600 9000; 4000 9000];

% arrays to store lockdown metrics
[Days_lockdown, Peak_incidence, Last_lockdown, NPhases] = deal([]);

figure('Position',[400 400 1000 1000])

for strat = 1:6
    % Define model parameters as a structure
    para = struct('beta',beta,'gamma',gamma,'sigma',sigma,'N',N,'n',n,'strategy',strat,'init',1);

    % dummy thresholds to allow infections to build with no intervention
    para.Imin = 20000;
    para.Imax = 20000;
    para.Imin_risk = para.Imin*atrisk_prop;
    para.Imax_risk = para.Imax*atrisk_prop;

    % Define initial conditions as a structure
    E0 = 0.0015;  % initial exposed
    ICs = struct('S',(1-E0).*para.N, 'E1',E0.*para.N, 'E2',zeros(n,1), 'E3',zeros(n,1), ...
                 'I',zeros(n,1), 'R',zeros(n,1));

    % Preliminary run - allows us to get ICs to begin in lockdown
    %[Prelim] = SEIR_demo_NBD(para,ICs,30);
    [Prelim] = SEIR_demo(para,ICs,25);

    % add control thresholds defined by strategy
    para.Imin = thresholds(strat,1);
    para.Imax = thresholds(strat,2);
    para.Imin_risk = para.Imin*atrisk_prop;
    para.Imax_risk = para.Imax*atrisk_prop;
    para.init = 0;

    % Get new initial conditions to begin in lockdown
    ICs = struct('S',Prelim.S(end,:), 'E1',Prelim.E1(end,:), 'E2',Prelim.E2(end,:), ...
                 'E3',Prelim.E3(end,:), 'I',Prelim.I(end,:), 'R',Prelim.R(end,:));

    % Run model
    %[Classes] = SEIR_demo_NBD(para,ICs,maxtime);
    [Classes] = SEIR_demo(para,ICs,maxtime);

    % Find times where social distancing is enforced
    nx = size(Classes.SD,1);

    % Compute daily new infections and "hospitalisations"
    %Tcases = sum(Classes.Cases,2);
    %Classes.NewInf = Tcases(2:end) - Tcases(1:end-1);
    %p_hosp = [0.01, 0.04, 0.20];
    %Classes.H = p_hosp.*Classes.I;

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
    plot(Classes.t, sum(Classes.I,2), 'k', 'LineWidth', 2.5)
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
    axis([0 maxtime 0 1.15*max([10000,max(sum(Classes.I,2))])])
    grid on

    % Compute metrics
    Peak_incidence = [Peak_incidence round(max(sum(Classes.I,2)))];
    %Days_over_x
    %Days_over_y
    Last_lockdown = [Last_lockdown Classes.SD(end,1)-1];
    Days_lockdown = [Days_lockdown DL];
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
save("./mats/metrics.mat","Days_lockdown","Last_lockdown","NPhases","Peak_incidence","strategies")