% basic cost function to evaluate

% weights: (w1,w2) determine how we weigh disease burden and lockdown
% stringency component
% para0: ODE parameter set for preliminary run
% para: ODE parameter set for main simulation
% NB: thresholds defined as para.Imin, para.Imax in para

function F = CostFunction_HH(weights, para, Hc, Peak_hospital, FinalHospital, Days_lockdown, Days_Tier2, DL, DT2)


Burden = FinalHospital/sum(para.N);
Stringency = 0.5*((0.7*Days_lockdown + 0.4*Days_Tier2)/para.maxtime)^2;

[Burden, Stringency];

F = weights(1)*Burden + weights(2)*Stringency + exp(weights(3)*(Peak_hospital-Hc));