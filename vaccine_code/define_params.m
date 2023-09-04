% model parameters
gamma = 1/7;                     % infectious period
epsilon = 1/5.28;                % incubation period
omega = 1/800;                   % recovered period (omega_S)
red = 4/3;                       % omega_A = red*omega
zeta = 1/8;                      % period between symptoms and hospital
delta = 1/10;                    % hospitalisation period
tau = 0.25;                      % relative infectiousness of asymptomatic
rho = 0.1;                       % relative infectiousness of hospitalised
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
hosp_rates = [0.1; 0.15; 0.3];   % probability of hospitalisation given symptomatic infection
N = [100000; 250000; 80000];     % population structure
n = size(N,1);                   % number of age classes
bedsper1000 = 2.5;               % UK hospital beds per 1,000 population
Hmax = bedsper1000*sum(N)/1000;  % capacity

% transmission matrix
beta = 0.7.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];

% initial exposed
E0 = 6e-4;

% control parameters
strat = 1;                       % default strategy
init = 0;                        % default control state (no control)
tgap = 18;                       % decision-making gap (relaxation)
tdelay = 3;                      % natural delay in implementing decision
tdiff = 7;                       % tgap-tdiff = decision-making gap (reintroduction)

% default switching thresholds
T01 = 10000;                     % No Control -> Intermediate Control
T10 = 10000;                     % Intermediate Control -> No Control
T12 = 10000;                     % Intermediate Control -> Lockdown
T21 = 10000;                     % Lockdown -> Intermediate Control

% Define time to run model for
maxtime = 30;                    % preliminary run

% generate fixed vaccine parameters
% parameters varying sigmoid curve (ie. rollout speed)
tc = 100;                        % time to complete 50% vaccination
kappa = 0.05;                    % logistic shaping parameter
stagger = tc/2;                  % time between vaccination commencement of age groups
nu_a = [1;1;1];                  % vaccine transmission coefficient
efficacy = 0.9;                  % default vaccine efficacy

vstart_times = [180:60:1080];
vstart = 2*max(vstart_times);    % default vaccine arrival date

save("./mats/Parameters.mat","beta","gamma","epsilon","omega","red","zeta","delta",...
    "tau","rho","da","hosp_rates","N","n","bedsper1000","Hmax","E0","strat","init",...
    "tgap","tdelay","tdiff","T01","T10","T12","T21","maxtime","tc","kappa","stagger",...
    "nu_a","efficacy","vstart",'-mat')
