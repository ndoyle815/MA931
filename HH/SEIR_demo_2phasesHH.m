% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions

function [Classes] = SEIR_demo_2phasesHH(para,ICs)

%Run ODE using ODE45
SD = [0,para.init];
opts = odeset('RelTol',1e-6,'MaxStep',0.1);
[t, pop] = ode45(@diff_SEIR_model, [0:1:para.maxtime], [ICs.S ICs.E1 ICs.E2 ICs.E3 ICs.IA ICs.IS ICs.IPH ICs.IH ICs.RA ICs.RS ICs.Cases ICs.Hosp], opts, para);

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E1',pop(:,para.n+1:2*para.n),'E2',pop(:,2*para.n+1:3*para.n), ...
            'E3',pop(:,3*para.n+1:4*para.n),'IA',pop(:,4*para.n+1:5*para.n),'IS',pop(:,5*para.n+1:6*para.n), ...
            'IPH',pop(:,6*para.n+1:7*para.n),'IH',pop(:,7*para.n+1:8*para.n),'RA',pop(:,8*para.n+1:9*para.n), 'RS',pop(:,9*para.n+1:10*para.n), ...
            'Cases',pop(:,10*para.n+1:11*para.n),'Hosp',pop(:,11*para.n+1:12*para.n),'SD',SD,'t',t);


%Diff equations

function dPop = diff_SEIR_model(t,pop,para)

S  = pop(1 : para.n);
E1 = pop(para.n+1 : 2*para.n);
E2 = pop(2*para.n+1 : 3*para.n);
E3 = pop(3*para.n+1 : 4*para.n);
IA = pop(4*para.n+1 : 5*para.n);
IS = pop(5*para.n+1 : 6*para.n);
IPH = pop(6*para.n+1 : 7*para.n);
IH = pop(7*para.n+1 : 8*para.n);
RA  = pop(8*para.n+1 : 9*para.n);
RS = pop(9*para.n+1 : 10*para.n);
Cases = pop(10*para.n+1 : 11*para.n);
Hosp = pop(11*para.n+1 : 12*para.n);

%NB: We only record recoveries to determine final epidemic size based off
%symptomatic infections only

% "social distancing" control: 70% decrease in contact rates if total infections
% above a given threshold

if SD(end,2) == 0
    factor = 1;
    if sum(IH) > para.U12 && t - SD(end,1) > para.tgap - para.tdiff
        SD(end+1,:) = [t,2.5];
    elseif sum(IH) > para.U01 && t - SD(end,1) > para.tgap - para.tdiff
        SD(end+1,:) = [t,0.5];
    end
elseif SD(end,2) == 0.5
    factor = 1 - 0.4*SD(end-1,2); % = 1 if previous state 0 or 0.6 if previous state 1
    if t - SD(end,1) > para.tdelay
        SD(end+1,:) = [t,1-SD(end-1,2)]; % = 1 if previous state 0 or 0 if previous state 1
    end
elseif SD(end,2) == 1
    factor = 0.6;
    if sum(IH) > para.U12 && t - SD(end,1) > para.tgap - para.tdiff
        SD(end+1,:) = [t,1.5];
    elseif sum(IH) < para.L10 && t - SD(end,1) > para.tgap
        SD(end+1,:) = [t,0.5];
    end
elseif SD(end,2) == 1.5
    factor = 0.6 - 0.3*(SD(end-1,2)-1); % = 0.6 if previous state 1 or 0.3 if previous state 2
    if t - SD(end,1) > para.tdelay
        SD(end+1,:) = [t,3-SD(end-1,2)]; % = 2 if previous state 1 or 1 if previous state 2
    end
elseif SD(end,2) == 2
    factor = 0.3;
    if sum(IH) < para.L21 && t - SD(end,1) > para.tgap
        SD(end+1,:) = [t,1.5];
    end
elseif SD(end,2) == 2.5
    factor = 1;
    if t - SD(end,1) > para.tdelay
        SD(end+1,:) = [t,2];
    end
end

% ODE equations
dS = -factor.*S.*(para.beta*(para.tau.*IA + IS + IPH + para.rho.*IH))./para.N + para.omega.*RS + para.red*para.omega.*RA;
dE1 = factor.*S.*(para.beta*(para.tau.*IA + IS + IPH + para.rho.*IH))./para.N - 3*para.sigma.*E1;
dE2 = 3*para.sigma.*E1 - 3*para.sigma.*E2;
dE3 = 3*para.sigma.*E2 - 3*para.sigma.*E3;
dIA = 3*para.sigma.*(1-para.da).*E3 - para.gamma.*IA;
dIS = 3*para.sigma.*para.da.*(1-para.hosp_rates).*E3 - para.gamma.*IS;
dIPH = 3*para.sigma.*para.da.*para.hosp_rates.*E3 - para.epsilon.*IPH;
dIH = para.epsilon.*IPH - para.delta.*IH;
dRA = para.gamma.*(IA) - para.red*para.omega.*RA;
dRS  = para.gamma.*(IS) + para.delta.*IH - para.omega.*RS;
dCases = 3*para.sigma.*para.da.*E3;
dHosp = para.epsilon.*IPH;

dPop = [dS; dE1; dE2; dE3; dIA; dIS; dIPH; dIH; dRA; dRS; dCases; dHosp];


end

end