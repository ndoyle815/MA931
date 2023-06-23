clear all

load('../mats/metrics.mat')

format long g
Ss = ['S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6'];
mymetrics = [Peak_incidence', Peak_hospital', FinalSize', FinalHospital', Last_lockdown', Days_lockdown', NPhases']