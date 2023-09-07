% Script to evaluate control strategies with vaccination, taking the
% arrival date distributions into account and producing Figure 5+6 from the
% report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',18)

% load default parameters
para0 = load('./mats/Parameters.mat');

% vaccination start times
vstart_times = [180:10:1080];
vstarts = [2*max(vstart_times), 360];

% Define time to run model for
t_init = 30;     % preliminary run
maxtime = max(vstart_times);  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 150 100 700; 50 150 100 200; 150 350 500 650; 275 350 425 500];
strategies = [1:length(thresholds)];

% add control thresholds defined by strategy
para = para0;
para.maxtime = maxtime;
para.Hmax = 1500;        % modify hospital capacity

% define vaccine efficacy
para.efficacy = 0.9;

% define functional weights
weights = [0.1:0.01:0.9];
w3 = 2;

% stores cost function outputs
ns = length(strategies);
nw = length(weights);
nv = length(vstart_times);
fs = zeros(nw,nv,ns);

tic
for strat = strategies
    % set switching thresholds
    para.T10 = thresholds(strat,1);
    para.T01 = thresholds(strat,2);
    para.T21 = thresholds(strat,3);
    para.T12 = thresholds(strat,4);

    % run preliminary simulation to get ICs
    [Prelim, Prelim_ICs] = Get_ICs(para0);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.T12
        para.init = 1;
    else
        para.init = 2;
    end

    for v = 1:nv
        para.vstart = vstart_times(v);

        % Run main simulation
        [Classes] = ODEmodel(para, Prelim_ICs);

        % use post-processor to compute metrics of interest
        [Classes, ~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, ~, ~, ~, ~, ~] = PostProcessor(Classes);
    
        for w = 1:nw
            % evaluate cost function
            fs(w,v,strat) = CostFunction([weights(w), 1-weights(w), w3], para, Peak_hospital, FinalHospital, Days_lockdown, Days_Tier2);
        end
    end
end
toc

%% Plotting
set(0,'defaultaxesfontsize',18)

f = figure(1);
f.Position = [100 400 1200 400];

cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

subplot(1,3,[1 2])
ax = gca;
ax.Position(2) = ax.Position(2) + 0.1;
ax.Position(1) = ax.Position(1) - 0.03;
ax.Position(4) = ax.Position(4) - 0.025;

for strat = strategies
    fstrat = fs(:,:,strat);
    h1 = surf(fstrat,'FaceColor',cols(strat,:),'LineStyle','none');
    hold on
    
    %%Extract X,Y and Z data from surface plot
    x=h1.XData;
    y=h1.YData;
    z=h1.ZData;
    %%Create vectors out of surface's XData and YData
    x=x(1,:);
    y=y(:,1);
    %%Divide the lengths by the number of lines needed
    xnumlines = 16; % 10 lines
    ynumlines = 17; % 10 partitions
    xspacing = round(length(x)/xnumlines);
    yspacing = round(length(y)/ynumlines);
    %%Plot the mesh lines 
    % Plotting lines in the X-Z plane
    for i = 1:yspacing:length(y)
        Y1 = y(i)*ones(size(x)); % a constant vector
        Z1 = z(i,:);
        plot3(x,Y1,Z1,'-k','LineWidth',0.5);
    end
    % Plotting lines in the Y-Z plane
    for i = 1:xspacing:length(x)
        X2 = x(i)*ones(size(y)); % a constant vector
        Z2 = z(:,i);
        plot3(X2,y,Z2,'-k','LineWidth',0.5);
    end
end
yticks([1:10:nw])
xticks([1:12:nv])
yticklabels(weights(1:10:end-5))
xticklabels(vstart_times(1:12:end))
ylim([1 nw])
xlim([1 nv])
label_y = ylabel('Weight $w_1$');
label_x = xlabel('Vaccine arrival time (days)');
%label_x.Position(1) = label_x.Position(1) + 0.01;
%label_x.Position(2) = label_x.Position(2) - 1;
%zt = get(gca, 'ZTick');
%set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))
f.CurrentAxes.ZDir = 'Reverse';
zlabel('Cost')

subplot(1,3,3)
ax = gca;
ax.Position(1) = ax.Position(1) + 0.03;
ax.Position(2) = ax.Position(2) + 0.06;
ax.Position(4) = ax.Position(4) - 0.02;

for strat = strategies
    fstrat = fs(:,:,strat);
    h2(strat) = surf(fstrat,'FaceColor','none','FaceColor',cols(strat,:),'LineStyle','none');
    hold on

    %%Extract X,Y and Z data from surface plot
    h2strat = h2(strat);
    x=h2strat.XData;
    y=h2strat.YData;
    z=h2strat.ZData;
    %%Create vectors out of surface's XData and YData
    x=x(1,:);
    y=y(:,1);
    %%Divide the lengths by the number of lines needed
    xnumlines = 16; % 10 lines
    ynumlines = 17; % 10 partitions
    xspacing = round(length(x)/xnumlines);
    yspacing = round(length(y)/ynumlines);
    %%Plot the mesh lines 
    % Plotting lines in the X-Z plane
    for i = 1:yspacing:length(y)
        Y1 = y(i)*ones(size(x)); % a constant vector
        Z1 = z(i,:);
        plot3(x,Y1,Z1,'-k','LineWidth',0.5);
    end
    % Plotting lines in the Y-Z plane
    for i = 1:xspacing:length(x)
        X2 = x(i)*ones(size(y)); % a constant vector
        Z2 = z(:,i);
        plot3(X2,y,Z2,'-k','LineWidth',0.5);
    end
end
view([0 270])
f.CurrentAxes.YDir = 'Reverse';
yticks([1:10:nw])
xticks([1:18:nv])
xtickangle(0)
yticklabels(weights(1:10:end))
xticklabels(vstart_times(1:18:end))
ylim([1 nw])
xlim([1 nv])
ylabel('Weight $w_1$')
label_x = xlabel('Vaccine arrival time (days)');
label_x.Position(2) = label_x.Position(2) - 0.6;
if para.efficacy < 0.4
    legend(h2, 'S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','Location','south','Interpreter','latex')
end

saveas(f,strcat('./vacc_images/optstrat_',num2str(para.efficacy),'.png'));

%% retrieving fs on coarser grid to compute expectations
weights2 = weights([1:5:nw]);
vstart_times2 = vstart_times([1:6:nv]);

for strat = strategies
    newfs(:,:,strat) = fs([1:5:nw],[1:6:nv],strat);
end

%% Computing expectations

% vaccine arrival date distributions
vdists = [discretenormal([0:length(vstart_times2)-1],3,1.5); discretenormal([0:length(vstart_times2)-1],7.5,1.5); ...
          discretenormal([0:length(vstart_times2)-1],12,1.5); discretenormal([0:length(vstart_times2)-1],7.5,1); ...
          discretenormal([0:length(vstart_times2)-1],7.5,2)];
diststr = {'$(\mu,\sigma) = (3,1.5)$','$(\mu,\sigma) = (7.5,1.5)$','$(\mu,\sigma) = (12,1.5)$', ...
           '$(\mu,\sigma) = (7.5,1)$','$(\mu,\sigma) = (7.5,2)$'};

%whichdists = [1 2 3];
%whichdists = [4 2 5];
whichdists = [5];
ndists = length(whichdists);

markers = {'o','^','s','d'};

f = figure(1+whichdists);
if whichdists == 3
    f.Position = [800 500 1070 300];
else
    f.Position = [800 100 750 300];
end


for i = 1:ndists
    wd = whichdists(i);
    vstart_dist = vdists(wd,:);

    for strat = strategies
        % calculate expected costs with uncertainty
        Ecosts(:,strat) = sum(newfs(:,:,strat).*vstart_dist,2);
        Varcosts(:,strat) = sum((newfs(:,:,strat).^2).*vstart_dist,2) - Ecosts(:,strat).^2;
        SD2costs(:,strat) = 2.*(Varcosts(:,strat)).^0.5;
    end

    subplot(ndists,1,i)
    
    for strat = [4 3 1 2]
        e = errorbar(weights2,Ecosts(:,strat),SD2costs(:,strat),"Color",cols(strat,:),'Marker',markers{strat}, ...
                     'MarkerSize',10,'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5);
        %hold on
        %s = scatter(weights,Ecosts(:,strat),60,"Color",cols(strat,:),'Marker',markers{strat}, ...
        %            'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5);
        hold on
    end
    axis([min(weights2)-0.025 max(weights2)+0.025 0.9*min(Ecosts-SD2costs,[],'all') 1.1*max(Ecosts+SD2costs,[],'all')])
    xticks([min(weights2):0.2:max(weights2)])
    ylabel('Cost')
    %title(diststr{wd})
    %if wd == 5
        xlabel('Weight $w_1$')
    %end
    if wd == 3
        legend('S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','FontSize',17.6,'Location','eastoutside','Interpreter','latex')
    end
    grid on
    set(gca,'TickLength',[0 0])

end
saveas(f,strcat('./vacc_images/dist',num2str(wd),'_',num2str(para.efficacy),'.png'))

% expected arrival time
Exp_time = round(sum(vstart_times2.*vstart_dist));
