% basic cost function to evaluate

% weights: (w1,w2) determine how we weigh disease burden and lockdown
% stringency components
% para0: ODE parameter set for preliminary run
% para: ODE parameter set for main simulation
% NB: thresholds defined as para.Imin, para.Imax in para

function F = CostFunction(weights, para0, para)

% run preliminary simulation to get ICs
[Prelim, ICs] = Get_ICs(para0);

% Run main simulation
[Classes] = SEIR_demo(para,ICs);

% use post-processor to compute metrics of interest
[Peak_incidence, Peak_hospital, FinalSize, FinalHospital, Last_lockdown, Days_lockdown, NPhases, nx] = PostProcessor(Prelim, Classes, para, para0.maxtime);

Burden = weights(1)*Peak_hospital;
Stringency = weights(2)*Days_lockdown;

[Burden, Stringency]

F = Burden + Stringency;