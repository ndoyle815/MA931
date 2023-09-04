% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions
% control and vaccination incorporated at discrete timesteps

function [Classes] = ODEmodel(para,ICs)

%unpack ICs
pop = zeros(1,12*para.n);
pop(1,1:para.n) = ICs.S;
pop(1,1*para.n+1:2*para.n)   = ICs.E1;
pop(1,2*para.n+1:3*para.n)   = ICs.E2;
pop(1,3*para.n+1:4*para.n)   = ICs.E3;
pop(1,4*para.n+1:5*para.n)   = ICs.IA;
pop(1,5*para.n+1:6*para.n)   = ICs.IS;
pop(1,6*para.n+1:7*para.n)   = ICs.IPH;
pop(1,7*para.n+1:8*para.n)   = ICs.IH;
pop(1,8*para.n+1:9*para.n)   = ICs.RA;
pop(1,9*para.n+1:10*para.n)  = ICs.RS;
pop(1,10*para.n+1:11*para.n) = ICs.Cases;
pop(1,11*para.n+1:12*para.n) = ICs.Hosp;
pop(1,12*para.n+1:13*para.n) = ICs.V;

%setup
tn = 0;
SD = [tn, para.init];  % This class will record date and control index at every switch
opts = odeset('RelTol',1e-6);%,'MaxStep',0.1);

% number to be vaccinated
Vmax = para.efficacy.*para.N';

% the main iteration
while tn < para.maxtime
    
    % begin vaccine-associated decrease in transmission
    if tn >= para.vstart
        N_eligible = pop(end,1:para.n) + pop(end,8*para.n+1:9*para.n) + pop(end,9*para.n+1:10*para.n);
        % do vaccine movements
        vacc = [LS(tn-para.vstart, para.kappa, para.tc+2*para.stagger), ...
               LS(tn-para.vstart, para.kappa, para.tc+para.stagger), ...
               LS(tn-para.vstart, para.kappa, para.tc)].*Vmax;
        dvacc = vacc - pop(end,12*para.n+1:13*para.n);
        Vacc_S = dvacc.*(pop(end,1:para.n)./N_eligible);
        Vacc_RA = dvacc.*(pop(end,8*para.n+1:9*para.n)./N_eligible);
        Vacc_RS = dvacc.*(pop(end,9*para.n+1:10*para.n)./N_eligible);
    
        pop(end,1:para.n) = pop(end,1:para.n) - Vacc_S;
        pop(end,8*para.n+1:9*para.n) = pop(end,8*para.n+1:9*para.n) - Vacc_RA;
        pop(end,9*para.n+1:10*para.n) = pop(end,9*para.n+1:10*para.n) - Vacc_RS;
        pop(end,12*para.n+1:13*para.n) = pop(end,12*para.n+1:13*para.n) + Vacc_S + Vacc_RA + Vacc_RS;
    end
    

    % "social distancing" control: 70% decrease in contact rates if hospital 
    % occupancy is above a given threshold (lockdown) or a 40% decrease for softer
    % restrictions and smaller thresholds (Intermediate Control)
    if SD(end,2) == 0
        para.factor = 1;
        if sum(pop(end,7*para.n+1:8*para.n)) > para.T12 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 2.5];
        elseif sum(pop(end,7*para.n+1:8*para.n)) > para.T01 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 0.5];
        end
    elseif SD(end,2) == 0.5
        para.factor = 1 - 0.4*SD(end-1,2); % = 1 if previous state 0 or 0.6 if previous state 1
        if tn - SD(end,1) >= para.tdelay
            SD(end+1,:) = [tn, 1-SD(end-1,2)]; % = 1 if previous state 0 or 0 if previous state 1
        end
    elseif SD(end,2) == 1
        para.factor = 0.6;
        if sum(pop(end,7*para.n+1:8*para.n)) > para.T12 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 1.5];
        elseif sum(pop(end,7*para.n+1:8*para.n)) < para.T10 && tn - SD(end,1) >= para.tgap
            SD(end+1,:) = [tn, 0.5];
        end
    elseif SD(end,2) == 1.5
        para.factor = 0.6 - 0.3*(SD(end-1,2)-1); % = 0.6 if previous state 1 or 0.3 if previous state 2
        if tn - SD(end,1) >= para.tdelay
            SD(end+1,:) = [tn, 3-SD(end-1,2)]; % = 2 if previous state 1 or 1 if previous state 2
        end
    elseif SD(end,2) == 2
        para.factor = 0.3;
        if sum(pop(end,7*para.n+1:8*para.n)) < para.T21 && tn - SD(end,1) >= para.tgap
            SD(end+1,:) = [tn, 1.5];
        end
    elseif SD(end,2) == 2.5
        para.factor = 1;
        if tn - SD(end,1) >= para.tdelay
            SD(end+1,:) = [tn, 2];
        end
    end
    
    inits = [pop(end,1:para.n) pop(end,para.n+1:2*para.n) pop(end,2*para.n+1:3*para.n) ...
            pop(end,3*para.n+1:4*para.n) pop(end,4*para.n+1:5*para.n) pop(end,5*para.n+1:6*para.n) ...
            pop(end,6*para.n+1:7*para.n) pop(end,7*para.n+1:8*para.n) pop(end,8*para.n+1:9*para.n) ...
            pop(end,9*para.n+1:10*para.n) pop(end,10*para.n+1:11*para.n) pop(end,11*para.n+1:12*para.n) ...
            pop(end,12*para.n+1:13*para.n)];
    
    % Run ODE using ODE45
    [t, newpop] = ode45(@diff_SEIR_model, [tn:tn+1], inits, opts, para);
    
    pop = [pop; newpop(end,:)];
    tn = tn + 1;

end

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E1',pop(:,para.n+1:2*para.n),'E2',pop(:,2*para.n+1:3*para.n), ...
            'E3',pop(:,3*para.n+1:4*para.n),'IA',pop(:,4*para.n+1:5*para.n),'IS',pop(:,5*para.n+1:6*para.n), ...
            'IPH',pop(:,6*para.n+1:7*para.n),'IH',pop(:,7*para.n+1:8*para.n),'RA',pop(:,8*para.n+1:9*para.n), 'RS',pop(:,9*para.n+1:10*para.n), ...
            'Cases',pop(:,10*para.n+1:11*para.n),'Hosp',pop(:,11*para.n+1:12*para.n),'V',pop(:,12*para.n+1:13*para.n),'SD',SD,'t',[0:tn]);


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
RA = pop(8*para.n+1 : 9*para.n);
RS = pop(9*para.n+1 : 10*para.n);
Cases = pop(10*para.n+1 : 11*para.n);
Hosp = pop(11*para.n+1 : 12*para.n);
V = pop(12*para.n+1 : 13*para.n);

% ODE equations
dS = -para.factor.*S.*(para.nu_a.*para.beta*(para.tau.*IA + IS + IPH + para.rho.*IH))./para.N + para.omega.*RS + para.red*para.omega.*RA;
dE1 = para.factor.*S.*(para.nu_a.*para.beta*(para.tau.*IA + IS + IPH + para.rho.*IH))./para.N - 3*para.epsilon.*E1;
dE2 = 3*para.epsilon.*E1 - 3*para.epsilon.*E2;
dE3 = 3*para.epsilon.*E2 - 3*para.epsilon.*E3;
dIA = 3*para.epsilon.*(1-para.da).*E3 - para.gamma.*IA;
dIS = 3*para.epsilon.*para.da.*(1-para.hosp_rates).*E3 - para.gamma.*IS;
dIPH = 3*para.epsilon.*para.da.*para.hosp_rates.*E3 - para.zeta.*IPH;
dIH = para.zeta.*IPH - para.delta.*IH;
dRA = para.gamma.*(IA) - para.red*para.omega.*RA;
dRS  = para.gamma.*(IS) + para.delta.*IH - para.omega.*RS;
dCases = 3*para.epsilon.*para.da.*E3;
dHosp = para.zeta.*IPH;
dV = [0; 0; 0];  % S and R individuals already moved in discretised manner

dPop = [dS; dE1; dE2; dE3; dIA; dIS; dIPH; dIH; dRA; dRS; dCases; dHosp; dV];


end

% Logistic sigmoid used to calculate "decrease in transmission due to
% vaccine rollout"
function f = LS(t,A,c)

f = 1/(1+exp(-A*(t-c)));

end

end