% basic cost function to evaluate

% weights: (w1,w2) determine how we weigh disease burden and lockdown
% stringency component
% para: ODE parameter set for main simulation
% NB: thresholds defined as para.Imin, para.Imax in para
% Hospital capacity defined as para.Hmax in para

function F = CostFunction(weights, para, Peak_hospital, FinalHospital, Days_lockdown, Days_Tier2)


Burden = FinalHospital/sum(para.N);
Stringency = 0.5*((0.7*Days_lockdown + 0.4*Days_Tier2)/para.maxtime)^2;


F = weights(1)*Burden + weights(2)*Stringency + exp(weights(3)*(Peak_hospital-para.Hmax));