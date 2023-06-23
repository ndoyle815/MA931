clear all

load('../mats/Distributions.mat')
load('../mats/Probabilities.mat')

f = figure(1);
f.Position = [200 400 900 600];
bar(Distribution_Symptoms_to_Hospital(2:end))
xlabel('Time (days)')
ylabel('Probability')
title('Lag from symptomatic infection to hospitalisation, $D^{S \rightarrow H}(t)$')

saveas(f,strcat('../images/Dist_StoH.png'))