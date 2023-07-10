% function to retrieve epidemiological metrics from simulation output

function [Peak_incidence, Peak_hospital, FinalSize, FinalHospital, Last_lockdown, Days_lockdown, NPhases, nx] = PostProcessor_HH(Classes)

% Find times where social distancing is enforced
nx = size(Classes.SD,1);

% Compute daily new infections, hospitalisations and deaths
death_rates = [0.001; 0.01; 0.1];

DL = 0;

for i = 1:4:nx
    DL = DL + Classes.SD(i+2,1) - Classes.SD(i,1);
end

% Compute metrics
Peak_incidence = round(max(sum(Classes.IS,2) + sum(Classes.IPH,2) + sum(Classes.IH,2)));
Peak_hospital = round(max(sum(Classes.IH,2)));
FinalSize = round(sum(Classes.Cases(end,:)));
FinalHospital = round(sum(Classes.Hosp(end,:)));
FinalDeaths = round(Classes.Cases(end,:)*death_rates);
Last_lockdown = round(Classes.SD(end,1));
Days_lockdown = round(DL);
NPhases = ceil(nx/2);
    
% save metrics to table
save("../mats/PostProcess_metrics.mat","Days_lockdown","Last_lockdown","NPhases","Peak_incidence","Peak_hospital", ...
    "FinalSize","FinalHospital","FinalDeaths")