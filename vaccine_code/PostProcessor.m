% function to retrieve epidemiological metrics from simulation output

function [Classes, Peak_incidence, Peak_hospital, FinalSize, FinalHospital, Last_lockdown, Days_lockdown, Days_Tier2, NLockdowns, nx, ix1, ix2, DL, DT2] = PostProcessor(Classes)

% Compute daily new infections, hospitalisations and deaths
death_rates = [0.001; 0.01; 0.1];

% Find times where social distancing is enforced
nx = size(Classes.SD,1);
ix1 = find(Classes.SD(:,2)==1);
ix2 = find(Classes.SD(:,2)==2);

% append SD for lockdown computation if we end in a restriction
if Classes.SD(end,1) ~= 0
    Classes.SD(end+1:end+2,:) = [Classes.t(end) 0; Classes.t(end) 0];
end

j2 = 1;
j1 = 1;
for i = ix2'
    DL(j2) = Classes.SD(i+2,1) - Classes.SD(i,1);
    j2 = j2 + 1;
end
for i = ix1'
    DT2(j1) = Classes.SD(i+2,1) - Classes.SD(i,1);
    j1 = j1 + 1;
end

% Compute metrics
Peak_incidence = round(max(sum(Classes.IS,2) + sum(Classes.IPH,2) + sum(Classes.IH,2)));
Peak_hospital = round(max(sum(Classes.IH,2)));
FinalSize = round(sum(Classes.Cases(end,:)));
FinalHospital = round(sum(Classes.Hosp(end,:)));
FinalDeaths = round(Classes.Cases(end,:)*death_rates);
Last_lockdown = round(Classes.SD(end,1));
Days_lockdown = round(sum(DL));
Days_Tier2 = round(sum(DT2));
NLockdowns = length(find(Classes.SD(:,2)==2));
    
% save metrics to table
save("./mats/PostProcess_metrics.mat","Days_lockdown","Last_lockdown","NLockdowns","Peak_incidence","Peak_hospital", ...
    "FinalSize","FinalHospital","FinalDeaths")