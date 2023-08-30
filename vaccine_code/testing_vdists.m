clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',18)

% model parameters
gamma = 1/7;                     % infectious period
sigma = 1/5.28;                  % incubation period
omega = 1/800;                   % recovered period
red = 4/3;
tau = 0.25;                      % relative infectiousness of asymptomatic
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
N = [100000; 250000; 80000];     % population structure
n = size(N,1);                   % number of age classes
atrisk_prop = N(end)/sum(N);     % at-risk proportion of population
bedsper1000 = 2.5;               % UK hospital beds per 1,000 population
Hmax = bedsper1000*sum(N)/1000;  % capacity
strat = 0;

% transmission matrix
beta = 0.7.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];

% vaccine transmission coefficient
nu_a = [1;1;1];

% generate vaccine parameters
efficacy = 0.9;
vstart_times = [180:60:1080];

% parameters varying sigmoid curve (ie. rollout speed)
tc = 100;
kappa = 0.05;
stagger = tc/2;


% Define time to run model for
t_init = 30;    % preliminary run
maxtime = max(vstart_times);  % main simulation

% define strategy numbers and switching thresholds
%thresholds = [50 150 100 700; 50 150 100 200; 50 150 550 650; 275 350 425 500];
thresholds = [50 150 100 700; 50 150 100 200; 150 350 500 650; 275 350 425 500];
strategies = [1:length(thresholds)];

% Define model parameters as a structure
para0 = struct('beta',beta,'gamma',gamma,'sigma',sigma,'omega',omega,'tau',tau, ...
               'da',da,'N',N,'n',n,'strategy',strat,'init',0,'maxtime',t_init, ...
               'tgap',18,'tdelay',3,'tdiff',7,'hosp_rates',[0.1; 0.15; 0.3], ...
               'epsilon',1/8,'delta',1/10,'rho',0.1,'red',red,'nu_a',nu_a, ...
               'efficacy',efficacy,'kappa',kappa,'tc',tc,'stagger',stagger);

% dummy thresholds to allow infections to build with no intervention
para0.U12 = 20000;
para0.U01 = 20000;
para0.L10 = 20000;
para0.L21 = 20000;
para0.vstart = 2*maxtime;

% parameters for main simulation
para = para0;
para.maxtime = maxtime;

% define functional weights
weights = [0.1:0.05:0.9];
w3 = 2;
Hc = 1500;

ns = length(strategies);
nw = length(weights);
nv = length(vstart_times);

% stores cost function outputs
fs = zeros(nw,nv,ns);

tic
for strat = strategies
    % set switching thresholds
    para.L10 = thresholds(strat,1);
    para.U01 = thresholds(strat,2);
    para.L21 = thresholds(strat,3);
    para.U12 = thresholds(strat,4);

    % run preliminary simulation to get ICs
    [Prelim, ICs] = Get_ICs_HH(para0);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.U12
        para.init = 1;
    else
        para.init = 2;
    end

    for v = 1:nv
        para.vstart = vstart_times(v);

        % Run main simulation
        [Classes] = SEIR_demo_discretised_vacc2(para,ICs);

        % use post-processor to compute metrics of interest
        [~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, ~, ~, ~] = PostProcessor_HH(Classes);
    
        for w = 1:nw
            % evaluate cost function
            fs(w,v,strat) = CostFunction_HH([weights(w), 1-weights(w), w3], para, Hc, Peak_hospital, FinalHospital, Days_lockdown, Days_Tier2);
        end
    end
end
toc

%% Plotting
set(0,'defaultaxesfontsize',18)

f = figure(1);
f.Position = [100 400 1100 500];

cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

subplot(1,3,[1 2])
ax = gca;
ax.Position(2) = ax.Position(2) + 0.03;

for strat = strategies
    h1 = surf(-fs(:,:,strat),'FaceColor','none','FaceColor',cols(strat,:));
    hold on
end
yticks([1:2:nw])
xticks([1:2:nv])
yticklabels(weights(1:2:end))
xticklabels(vstart_times(1:2:end))
ylim([1 nw])
xlim([1 nv])
ylabel('Weight w1')
label_x = xlabel('Vaccine arrival time (days)');
label_x.Position(1) = label_x.Position(1) + 3;
label_x.Position(2) = label_x.Position(2) + 1;
zlabel('Negative cost')

subplot(1,3,3)
ax = gca;
ax.Position(1) = ax.Position(1) + 0.03;
ax.Position(2) = ax.Position(2) + 0.06;
ax.Position(4) = ax.Position(4) - 0.1;

for strat = strategies
    h2 = surf(-fs(:,:,strat),'FaceColor','none','FaceColor',cols(strat,:));
    hold on
end
view([0 90])
yticks([1:2:nw])
xticks([2:3:nv])
xtickangle(0)
yticklabels(weights(1:2:end))
xticklabels(vstart_times(2:3:end))
ylim([1 nw])
xlim([1 nv])
ylabel('Weight w1')
label_x = xlabel('Vaccine arrival time (days)');
label_x.Position(2) = label_x.Position(2) - 0.6;

if efficacy < 0.4
    legend('S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','Fontsize',14,'Location','south','Interpreter','latex')
end

saveas(f,strcat('./vacc_images/optstrat_',num2str(efficacy),'.png'));

%% Computing expectations

%vstart_unifdist = 1/length(vstart_times).*ones(1,length(vstart_times));  % uniform distribution
%vstart_poisdist = poisspdf([0:length(vstart_times)-1],3);  % poisson distribution
%vdists = [vstart_unifdist; vstart_poisdist];
vdists = [discretenormal([0:length(vstart_times)-1],3,1.5); discretenormal([0:length(vstart_times)-1],7.5,1.5); ...
          discretenormal([0:length(vstart_times)-1],12,1.5); discretenormal([0:length(vstart_times)-1],7.5,1); ...
          discretenormal([0:length(vstart_times)-1],7.5,2)];
diststr = {'$(\mu,\sigma) = (3,1.5)$','$(\mu,\sigma) = (7.5,1.5)$','$(\mu,\sigma) = (12,1.5)$', ...
           '$(\mu,\sigma) = (7.5,1)$','$(\mu,\sigma) = (7.5,2)$'};

whichdists = [1 2 3];
%whichdists = [4 2 5];
whichdists = [1];
ndists = length(whichdists);

markers = {'o','^','s','d'};

f = figure(2);
f.Position = [600 800 900 350];

for i = 1:ndists
    wd = whichdists(i);
    vstart_dist = vdists(wd,:);

    for strat = strategies
        Ecosts(:,strat) = sum(fs(:,:,strat).*vstart_dist,2);
        Varcosts(:,strat) = sum((fs(:,:,strat).^2).*vstart_dist,2) - Ecosts(:,strat).^2;
        SD2costs(:,strat) = 2.*(Varcosts(:,strat)).^0.5;
    end
    
    %subplot(2,1,1)
    %for strat = strategies
    %    h = plot(weights,Ecosts(:,strat),"Color",cols(strat,:),'LineWidth',2.5,'LineStyle',linetype{i});
    %    hold on
    %end
    %axis([min(weights) max(weights) 0.02 0.1])
    %xlabel('Weight w1')
    %ylabel('Expected Cost')
    %title(strcat( diststr(whichdists), ', vaccine efficacy =',{' '}, num2str(efficacy) ))
    %title(strcat('vaccine efficacy =',{' '}, num2str(efficacy)))
    %legend('S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','Location','eastoutside','Interpreter','latex')
    %grid on

    %subplot(2,1,2)
    %for strat = strategies
    %    h = plot(weights,SD2costs(:,strat),"Color",cols(strat,:),'LineWidth',2.5,'LineStyle',linetype{i});
    %    hold on
    %end
    %axis([min(weights) max(weights) 0 0.1])
    %xlabel('Weight w1')
    %ylabel('Variance of Cost')
    %title(strcat( diststr(whichdists), ', vaccine efficacy =',{' '}, num2str(efficacy) ))
    %title(strcat('vaccine efficacy =',{' '}, num2str(efficacy)))
    %legend('S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','Location','northeast','Interpreter','latex')
    %grid on

    subplot(ndists,1,i)
    
    for strat = [4 3 1 2]
        e = errorbar(weights,Ecosts(:,strat),SD2costs(:,strat),"Color",cols(strat,:),'Marker',markers{strat}, ...
                     'MarkerSize',10,'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5);
        %hold on
        %s = scatter(weights,Ecosts(:,strat),60,"Color",cols(strat,:),'Marker',markers{strat}, ...
        %            'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5);
        hold on
    end
    axis([min(weights)-0.025 max(weights)+0.025 0.9*min(Ecosts-SD2costs,[],'all') 1.1*max(Ecosts+SD2costs,[],'all')])  
    ylabel('Cost')
    title(diststr{wd})
    if wd == 5
        xlabel('Weight w1')
    end
    if wd == 3
        legend('S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','Location','northeast','Interpreter','latex')
    end
    grid on
    set(gca,'TickLength',[0 0])

end
saveas(f,strcat('./vacc_images/dist',num2str(wd),'_',num2str(efficacy),'.png'))

% expected arrival time
Exp_time = round(sum(vstart_times.*vstart_dist));