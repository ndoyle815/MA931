% script to run simulations of control model with and without vaccination,
% produces Figure 1 (cumulative vaccinations) and Figure 3 (simulations) from the report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% load default parameters
para0 = load('./mats/Parameters.mat');

% vaccination start times
vstart_times = [180:60:1080];
vstarts = [2*max(vstart_times), 360];
vstarts = 360;
% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 800;  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 150 100 700; 50 150 100 200; 150 350 500 650; 275 350 425 500];
strategies = [1:length(thresholds)];

% plotting preperation
figure('Position',[200 400 500*length(vstarts) 1000])
stratnames = {'S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)'};
stpos = [50 200 350 700; 50 200 350 500; 150 350 500 650; 175 325 475 625];

tic
for vs = 1:length(vstarts)
    for strat = strategies
        % Preliminary run - no control, 30 day build-up
        para0.init = 0;
        para0.maxtime = t_init;
        para0.vstart = vstarts(vs);
        % Dummy thresholds to allow infections to build with no intervention
        para0.T12 = 20000;
        para0.T01 = 20000;
        para0.T10 = 20000;
        para0.T21 = 20000;    

        % run preliminary simulation
        [Prelim, Prelim_ICs] = Get_ICs(para0);
        
        % add control thresholds defined by strategy
        para = para0;
        para.maxtime = maxtime;
        para.T10 = thresholds(strat,1);
        para.T01 = thresholds(strat,2);
        para.T21 = thresholds(strat,3);
        para.T12 = thresholds(strat,4);

        % starting control state
        if sum(Prelim.IH(end,:)) < para.T12
            para.init = 1;
        else
            para.init = 2;
        end

        % Run model
        [Classes] = ODEmodel(para,Prelim_ICs);

        % Post-Process for epidemic metrics
        [Classes, ~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, nx, ix1, ix2, ~, ~] = PostProcessor(Classes);

        plotidx = length(vstarts)*(strat-1)+vs;  % index for subfigure

        % plotting
        subplot(4,length(vstarts),plotidx)
        yyaxis left
        ax1 = gca;
        ax1.Position(1) = ax1.Position(1) + 0.04;
        ax1.YColor = 'k';
        ax1.FontSize = 16;
        ax1.FontSizeMode = 'manual';

        for i = ix1'
            patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'y', 'Facealpha',0.3, 'EdgeAlpha',0)
            hold on
        end
        for i = ix2'
            patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
            hold on
        end
        plot(Classes.t, sum(Classes.IH,2), 'k', 'LineWidth', 2.5)
        if plotidx == length(vstarts)
            xline(para.vstart,'-','Vaccine Arrival','Color','b','Linewidth',2,'FontSize',14,'Interpreter','latex','LabelOrientation','horizontal')
        else
            xline(para.vstart,'-','Color','b','Linewidth',2)
        end
        
        axis([0 maxtime 0 max([para.Hmax,1200])])
    
        if plotidx > length(vstarts)*3
            xlabel('Time (days)')
        end
        %if mod(plotidx,2) == 1
            label_y = ylabel('$I^H(t)$','Rotation',0);
            label_y.Position(1) = label_y.Position(1) - 0;
        %end
    
        yyaxis right

        plot(Classes.t, para.T01.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
        hold on
        plot(Classes.t, para.T12.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
        hold on
        plot(Classes.t, para.T10.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
        hold on
        plot(Classes.t, para.T21.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
        yticks(stpos(strat,:))
    
        if para.T21 > para.T01
            yticklabels({'$T_{10}$','$T_{01}$','$T_{21}$','$T_{12}$'})
        else
            yticklabels({'$T_{10}$','$T_{21}$','$T_{01}$','$T_{12}$'})
        end
        axis([0 maxtime 0 max([para.Hmax,1200])])
        ax2 = gca;
        ax2.Position(3) = ax2.Position(3) - 0.02;
        ax2.YColor = 'k';
        ax2.TickDir = 'none';
        ax2Y = ax2.YAxis(2,1);
        ax2Y.FontSize = 14;
        ax2.FontSizeMode = 'manual';
    
        title(stratnames{strat})
        grid on
    
        % Evaluate cost function
        %F = CostFunction_HH([0.3 0.7 2], para, Peak_hospital, FinalHospital, Days_lockdown, Days_Tier2);
    end
end
toc

%save figure
saveas(gcf,strcat('./vacc_images/','vacc2_',num2str(para.vstart),'.png'))


% Plotting cumulative vaccinations
f = figure(2);
f.Position = [1250 400 600 400];
plot(Classes.t, Classes.V(:,1)./1000, 'LineWidth', 2.5, 'DisplayName', '0-19')
hold on
plot(Classes.t, Classes.V(:,2)./1000, 'LineWidth', 2.5, 'DisplayName', '20-64')
hold on
plot(Classes.t, Classes.V(:,3)./1000, 'LineWidth', 2.5, 'DisplayName', '65+')
hold on
plot(Classes.t, sum(Classes.V,2)./1000, 'k', 'LineWidth', 2.5, 'DisplayName', 'Total')
%axis([para.vstart-260 maxtime 0 440])

xline(para.vstart,'--','DisplayName','Arrival','Color','k')
xlabel('Time (days)')
ylabel('Population (thousands)')
title('Cumulative Vaccinations')
legend('Interpreter','latex','Location','west')
grid on

saveas(gcf,'./vacc_images/Tvacc.png')
