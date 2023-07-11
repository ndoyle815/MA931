% script to run a single simulation
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% model parameters
gamma = 1/7;                     % infectious period
sigma = 1/5;                     % latency period
omega = 0*1/400;                 % recovered period
tau = 0.25;                      % relative infectiousness of asymptomatic
da = [0.05; 0.2; 0.7];           % probability of symptomatic infection
N = 2.*[50000; 125000; 40000];   % population structure
n = size(N,1);                   % number of age classes
atrisk_prop = N(end)/sum(N);     % at-risk proportion of population
bedsper1000 = 2.5;               % UK hospital beds per 1,000 population
Hmax = bedsper1000*sum(N)/1000;  % capacity
strat = 0;

% transmission matrix
beta = 0.7.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 600;  % main simulation

% Define model parameters as a structure
para0 = struct('beta',beta,'gamma',gamma,'sigma',sigma,'omega',omega,'tau',tau,'da',da,...
               'N',N,'n',n,'strategy',strat,'init',0,'maxtime',t_init,'tgap',18,'tdelay',3,...
               'tdiff',7,'hosp_rates',[0.1; 0.15; 0.3],'epsilon',1/8,'delta',1/10,'rho',0.1);

% dummy thresholds to allow infections to build with no intervention
para0.U12 = 20000;
para0.U01 = 20000;
para0.L10 = 20000;
para0.L21 = 20000;

% add control thresholds defined by strategy
para = para0;
para.init = 1;
para.maxtime = maxtime;

% define strategy numbers and thresholds
strategies = [1:6];
thresholds = [50 200 125 700; 50 150 100 200; 100 300 400 550; 150 250 350 450; 50 300 400 500; 300 600 400 700];

% define functional weights
%weights = [0.5, 0.5, 0.1; 0.8, 0.2, 0.1; 1, 0, 0.1];
weights = [0:0.02:1];
w3 = 0.1;
Hc = Hmax;

ns = length(strategies);
nw = length(weights);

% stores cost function outputs
fs = zeros(ns,nw);
hs = zeros(ns,1);
dls = zeros(ns,1);

tic
for w = 1:nw
    for strat = strategies
        % set switching thresholds
        para.L10 = thresholds(strat,1);
        para.U01 = thresholds(strat,2);
        para.L21 = thresholds(strat,3);
        para.U12 = thresholds(strat,4);

        % evaluate cost function
        [fs(strat,w),hs(strat),dls(strat)] = CostFunction_HH([weights(w), 1-weights(w), w3], para0, para, Hc);
    end
end
toc

fs;
[~,idx] = min(fs)

figure(1)
scatter(weights,idx,30,'b','filled')
axis([0 1 0.5 6.5])
xlabel('w1 (w2 = 1 - w1)')
ylabel('Optimal strategy')
title(strcat('Cost function, Hc = ',num2str(Hc),', ','maxtime =  ',num2str(maxtime)))
grid on

saveas(gcf,strcat('../images/CostFunction_',num2str(Hc),'_',num2str(maxtime),'.png'))


% Hpeaks = [252, 164, 164, 241, 253, 339];
% 
% figure(2)
% clf
% bar(Hpeaks)
% xlabel('Strategy')
% title('Peak daily hospital admissions')
% grid on
% 
% saveas(gcf,'./images/Peakadmissions.png')

% mygreen = [0 0.5 0];
% figure(1)
% scatter([1:6],fs(:,1),'filled','color','b')
% hold on
% scatter([1:6],fs(:,2),'filled','color','r')
% hold on
% scatter([1:6],fs(:,3),'filled')
% axis([1 6 0.02 0.06])
% %set(gca, 'YScale', 'log')
% 
% fscaled = zeros(ns,nw);
% 
% [~,mins] = min(fs);
% [~,maxs] = max(fs);
% 
% for w = 1:nw
%     for strat = strategies
%         fscaled(strat,w) = (fs(strat,w) - fs(mins(w),w))/(fs(maxs(w),w) - fs(mins(w),w));
%     end
% end
% 
% figure(2)
% scatter([1:6],fscaled(:,1),'filled','color','b')
% hold on
% scatter([1:6],fscaled(:,2),'filled','color','r')
% hold on
% scatter([1:6],fscaled(:,3),'filled')
% 
% hscaled = zeros(ns,1);
% dlscaled = zeros(ns,1);
% 
% [~,hmin] = min(hs);
% [~,hmax] = max(hs);
% [~,dlmin] = min(dls);
% [~,dlmax] = max(dls);
% 
% for strat = strategies
%     hscaled(strat) = (hs(strat) - hs(hmin))/(hs(hmax) - hs(hmin));
%     dlscaled(strat) = (dls(strat) - dls(dlmin))/(dls(dlmax) - dls(dlmin));
% end
% 
% markers = {'o' 's' 'd' 'x' '*' '+'};
% 
% figure(3)
% hold on
% for s = 1:ns
%     plot(hscaled(s),dlscaled(s),'LineStyle','none','Marker',markers{s},'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10);
% end
% xlabel('Total Hospitalisations')
% ylabel('Days in Lockdown')
% grid on
% legend(["S1","S2","S3","S4","S5","S6"])
% axis([-0.1 1.1 -0.1 1.1])