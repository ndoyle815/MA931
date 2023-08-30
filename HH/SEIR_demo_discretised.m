% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions

function [Classes] = SEIR_demo_discretised(para,ICs)

%load('mats/Distributions.mat','Distribution_Symptoms_to_Hospital','Distribution_Hosp_Time')

%unpack ICs
pop = zeros(1,12*para.n);
pop(1,1:para.n) = ICs.S;
pop(1,1*para.n+1:2*para.n) = ICs.E1;
pop(1,2*para.n+1:3*para.n) = ICs.E2;
pop(1,3*para.n+1:4*para.n) = ICs.E3;
pop(1,4*para.n+1:5*para.n) = ICs.IA;
pop(1,5*para.n+1:6*para.n) = ICs.IS;
pop(1,6*para.n+1:7*para.n) = ICs.IPH;
pop(1,7*para.n+1:8*para.n) = ICs.IH;
pop(1,8*para.n+1:9*para.n) = ICs.RA;
pop(1,9*para.n+1:10*para.n) = ICs.RS;
pop(1,10*para.n+1:11*para.n) = ICs.Cases;
pop(1,11*para.n+1:12*para.n) = ICs.Hosp;

%setup
tn = 0;
SD = [tn, para.init];
opts = odeset('RelTol',1e-6);%,'MaxStep',0.1);

% stores hospitalisations
%HH = 0;
%StoHlag = [0:length(Distribution_Symptoms_to_Hospital)-1];
%HtoRlag = [0:length(Distribution_Hosp_Time)-1];

% the main iteration
while tn < para.maxtime
    
    % compute new hospitalisations
    %newHH = Distribution_Symptoms_to_Hospital*(pop(end-StoHlag,7*para.n+1:8*para.n)*para.hosp_rates);

    % "social distancing" control: 70% decrease in contact rates if total infections
    % above a given threshold
    if SD(end,2) == 0
        para.factor = 1;
        if sum(pop(end,7*para.n+1:8*para.n)) > para.U12 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 2.5];
        elseif sum(pop(end,7*para.n+1:8*para.n)) > para.U01 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 0.5];
        end
    elseif SD(end,2) == 0.5
        para.factor = 1 - 0.4*SD(end-1,2); % = 1 if previous state 0 or 0.6 if previous state 1
        if tn - SD(end,1) >= para.tdelay
            SD(end+1,:) = [tn, 1-SD(end-1,2)]; % = 1 if previous state 0 or 0 if previous state 1
        end
    elseif SD(end,2) == 1
        para.factor = 0.6;
        if sum(pop(end,7*para.n+1:8*para.n)) > para.U12 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 1.5];
        elseif sum(pop(end,7*para.n+1:8*para.n)) < para.L10 && tn - SD(end,1) >= para.tgap
            SD(end+1,:) = [tn, 0.5];
        end
    elseif SD(end,2) == 1.5
        para.factor = 0.6 - 0.3*(SD(end-1,2)-1); % = 0.6 if previous state 1 or 0.3 if previous state 2
        if tn - SD(end,1) >= para.tdelay
            SD(end+1,:) = [tn, 3-SD(end-1,2)]; % = 2 if previous state 1 or 1 if previous state 2
        end
    elseif SD(end,2) == 2
        para.factor = 0.3;
        if sum(pop(end,7*para.n+1:8*para.n)) < para.L21 && tn - SD(end,1) >= para.tgap
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
            pop(end,9*para.n+1:10*para.n) pop(end,10*para.n+1:11*para.n) pop(end,11*para.n+1:12*para.n)];
    
    % Run ODE using ODE45
    [t, newpop] = ode45(@diff_SEIR_model, [tn:tn+1], inits, opts, para);
    
    pop = [pop; newpop(end,:)];
    tn = tn + 1;

end

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E1',pop(:,para.n+1:2*para.n),'E2',pop(:,2*para.n+1:3*para.n), ...
            'E3',pop(:,3*para.n+1:4*para.n),'IA',pop(:,4*para.n+1:5*para.n),'IS',pop(:,5*para.n+1:6*para.n), ...
            'IPH',pop(:,6*para.n+1:7*para.n),'IH',pop(:,7*para.n+1:8*para.n),'RA',pop(:,8*para.n+1:9*para.n), 'RS',pop(:,9*para.n+1:10*para.n), ...
            'Cases',pop(:,10*para.n+1:11*para.n),'Hosp',pop(:,11*para.n+1:12*para.n),'SD',SD,'t',[0:tn]);


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

% ODE equations
dS = -para.factor.*S.*(para.beta*(para.tau.*IA + IS + IPH + para.rho.*IH))./para.N + para.omega.*RS + para.red*para.omega.*RA;
dE1 = para.factor.*S.*(para.beta*(para.tau.*IA + IS + IPH + para.rho.*IH))./para.N - 3*para.sigma.*E1;
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