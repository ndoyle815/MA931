% basic cost function to evaluate

% weights: (w1,w2) determine how we weigh disease burden and lockdown
% stringency component
% para0: ODE parameter set for preliminary run
% para: ODE parameter set for main simulation
% NB: thresholds defined as para.Imin, para.Imax in para

function [F, FinalHospital, Days_lockdown] = CostFunction_HH(weights, para0, para, Hc)

% run preliminary simulation to get ICs
[Prelim, ICs] = Get_ICs_HH(para0);

% Run main simulation
[Classes] = SEIR_demo_2phasesHH(para,ICs);

% use post-processor to compute metrics of interest
[Peak_incidence, Peak_hospital, FinalSize, FinalHospital, Last_lockdown, Days_lockdown, NPhases, nx] = PostProcessor_HH(Classes);

Burden = FinalHospital/sum(para.N);
Stringency = 0.5*(0.7*Days_lockdown/para.maxtime)^2;

[Burden, Stringency, Peak_hospital];

F = weights(1)*Burden + weights(2)*Stringency + exp(weights(3)*(Peak_hospital-Hc));