clear all

load('./mats/metrics.mat')

format long g
Ss = ['S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'];
mymetrics = [Peak_incidence', Peak_hospital', FinalSize', FinalHospital', Last_lockdown', Days_lockdown', NPhases']