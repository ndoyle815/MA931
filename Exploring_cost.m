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
beta = 1.0.*[1.35, 2.41, 0.3; 0.93, 3.01, 0.54; 0.35, 1.64, 1.21].*gamma;

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 400;  % main simulation

% Define model parameters as a structure
para0 = struct('beta',beta,'gamma',gamma,'sigma',sigma,'omega',omega,'tau',tau, ...
               'da',da,'N',N,'n',n,'strategy',strat,'init',0,'maxtime',t_init, ...
               'tgap',16,'tdelay',5);

% dummy thresholds to allow infections to build with no intervention
para0.Imin = 20000;
para0.Imax = 20000;
para0.Imin_risk = para0.Imin*atrisk_prop;
para0.Imax_risk = para0.Imax*atrisk_prop;

% add control thresholds defined by strategy
para = para0;
para.init = 1;
para.maxtime = maxtime;

% define strategy numbers and thresholds
strategies = [1:6];
thresholds = [1500 6000; 1000 2000; 2000 4000; 4000 6000; 1000 8000; 4000 8000];

% define functional weights
weights = [0.5, 0.5; 0.8, 0.2; 1, 0];

% stores cost function outputs
fs = zeros(length(strategies),size(weights,1));

for w = 1:size(weights,1)
    for strat = strategies
        para.Imin = thresholds(strat,1);
        para.Imax = thresholds(strat,2);
        para.Imin_risk = para.Imin*atrisk_prop;
        para.Imax_risk = para.Imax*atrisk_prop;

        % evaluate cost function
        fs(strat,w) = CostFunction(weights(w,:), para0, para);
    end
end
fs
figure(1)
plot(fs(:,1),'b')
hold on
plot(fs(:,2),'r')
hold on
plot(fs(:,3),'g')
