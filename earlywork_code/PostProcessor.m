% function to retrieve epidemiological metrics from simulation output

function [Peak_incidence, Peak_hospital, FinalSize, FinalHospital, Last_lockdown, Days_lockdown, NPhases, nx] = PostProcessor(Prelim, Classes, para, t_init)

load('mats/Distributions.mat','Distribution_Symptoms_to_Hospital')
load('mats/Probabilities.mat','Sympt_2_hosp')

% Find times where social distancing is enforced
nx = size(Classes.SD,1);

% Compute daily new infections, hospitalisations and deaths
hosp_rates = [0.1; 0.15; 0.3];
death_rates = [0.001; 0.01; 0.1];
    
Classes.DailyInf = [0 0 0; Prelim.Cases(2:end,:) - Prelim.Cases(1:end-1,:); Classes.Cases(2:end,:) - Classes.Cases(1:end-1,:)];
HH = Classes.DailyInf*hosp_rates;
newHH = zeros(para.maxtime+t_init+1,1);
lag = [0:length(Distribution_Symptoms_to_Hospital)-1];
for t = t_init+1:length(newHH)
    newHH(t) = newHH(t) + Distribution_Symptoms_to_Hospital*HH(t-lag);
end
Classes.Hospital = newHH(t_init+1:end);
Classes.Deaths = Classes.DailyInf*death_rates;

DL = 0;

for i = 1:4:nx
    DL = DL + Classes.SD(i+2,1) - Classes.SD(i,1);
end

% Compute metrics
Peak_incidence = round(max(sum(Classes.IS,2)));
Peak_hospital = round(max(Classes.Hospital));
FinalSize = round(sum(Classes.R(end,:)));
FinalHospital = round(sum(Classes.Hospital));
FinalDeaths = round(Classes.Cases(end,:)*death_rates);
Last_lockdown = round(Classes.SD(end,1));
Days_lockdown = round(DL);
NPhases = ceil(nx/2);
    
% save metrics to table
save("./mats/PostProcess_metrics.mat","Days_lockdown","Last_lockdown","NPhases","Peak_incidence","Peak_hospital", ...
    "FinalSize","FinalHospital","FinalDeaths")