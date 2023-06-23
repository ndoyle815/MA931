% script to run a single simulation
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% model parameters
gamma = 1/7;                     % infectious period
sigma = 1/5;                     % latency period
omega = 1/200;                   % recovered period
tau = 0.25;                      % relative infectiousness of asymptomatic
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
N = 3.*[50000; 125000; 40000];   % population structure
n = size(N,1);                   % number of age classes
atrisk_prop = N(end)/sum(N);     % at-risk proportion of population
strat = 1;

% transmission matrix
%load('./mats/contacts.mat')
beta = 1.0.*[1.35, 2.41, 0.3; 0.93, 3.01, 0.54; 0.35, 1.64, 1.21].*gamma;

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 600;  % main simulation

% Define model parameters as a structure
para0 = struct('beta',beta,'gamma',gamma,'sigma',sigma,'omega',omega,'tau',tau, ...
               'da',da,'N',N,'n',n,'strategy',strat,'init',0,'maxtime',t_init, ...
               'tgap',16,'tdelay',5);

% dummy thresholds to allow infections to build with no intervention
para0.Imin = 20000;
para0.Imax = 20000;
para0.Imin_risk = para0.Imin*atrisk_prop;
para0.Imax_risk = para0.Imax*atrisk_prop;

% run preliminary simulation to get ICs
[Prelim, ICs] = Get_ICs(para0);

% DEFINE RELAXATION AND REINTRODUCTION THRESHOLDS HERE
thresholds = [1000,2000];

% add control thresholds defined by strategy
para = para0;
para.init = 1;
para.maxtime = maxtime;
para.Imin = thresholds(1);
para.Imax = thresholds(2);
para.Imin_risk = para.Imin*atrisk_prop;
para.Imax_risk = para.Imax*atrisk_prop;

% Run main simulation
[Classes] = SEIR_demo(para,ICs);

% use post-processor to compute metrics of interest
[Peak_incidence, Peak_hospital, FinalSize, FinalHospital, Last_lockdown, Days_lockdown, NPhases, nx] = PostProcessor(Prelim, Classes, para, t_init);

% plotting
figure('Position',[400 400 800 600])
for i = 1:4:nx
    patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
    hold on
end
plot(Classes.t, sum(Classes.IS,2), 'k', 'LineWidth', 2.5)
hold on
yline(para.Imax,':','Upper')
hold on
yline(para.Imin,':','Lower')

xlabel('Time (days)')
ylabel('Infected Individuals')
title(strcat("Strategy S",num2str(strat)))
axis([0 maxtime 0 max([10000,1.15*max(sum(Classes.IS,2))])])
grid on

%save figure
saveas(gcf,strcat('./images/Strat_',num2str(thresholds(1)),'_',num2str(thresholds(2)),'.png'))

% evaluate cost function
weights = [1,1];
F = CostFunction(weights, para0, para)