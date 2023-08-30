clear all

load('../mats/Distributions.mat')
load('../mats/Probabilities.mat')

mu = 8;

f = figure(1);
f.Position = [100 400 750 500];
bar(Distribution_Symptoms_to_Hospital(2:end))
hold on
plot([0:30], exppdf([0:30],mu), 'r', 'LineWidth', 2.5)
xlabel('Time (days)')
ylabel('Probability')
title('Lag from symptomatic infection to hospitalisation, $D^{S \rightarrow H}(t)$')

saveas(f,strcat('../images/Dist_StoH.png'))

DHT = [1, Distribution_Hosp_Time] - [Distribution_Hosp_Time, 0]

mu = 10;

f = figure(2);
f.Position = [1000 400 750 500];
bar(DHT)
hold on
plot([0:60], exppdf([0:60],mu), 'r', 'LineWidth', 2.5)
xlabel('Time (days)')
ylabel('Probability')
title('Time spent in hospitalised class, $D^{H \rightarrow R}(t)$')

saveas(f,strcat('../images/Dist_HtoR.png'))