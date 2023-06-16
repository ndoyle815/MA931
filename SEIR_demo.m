% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions
% maxtime: maximum time to simulate

function [Classes] = SEIR_demo(para,ICs,maxtime)

%Run ODE using ODE45
SD = [0,1];
opts = odeset('RelTol',1e-5,'MaxStep',0.1);
[t, pop] = ode45(@diff_SEIR_model, [0:1:maxtime], [ICs.S ICs.E1 ICs.E2 ICs.E3 ICs.IS ICs.IA ICs.R ICs.Cases], opts, para);

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E1',pop(:,para.n+1:2*para.n),'E2',pop(:,2*para.n+1:3*para.n), ...
            'E3',pop(:,3*para.n+1:4*para.n),'IS',pop(:,4*para.n+1:5*para.n),'IA',pop(:,5*para.n+1:6*para.n), ...
            'R',pop(:,6*para.n+1:7*para.n),'Cases',pop(:,7*para.n+1:8*para.n),'SD',SD,'t',t);


%Diff equations

function dPop = diff_SEIR_model(t,pop,para)

S  = pop(1 : para.n);
E1 = pop(para.n+1 : 2*para.n);
E2 = pop(2*para.n+1 : 3*para.n);
E3 = pop(3*para.n+1 : 4*para.n);
IS = pop(4*para.n+1 : 5*para.n);
IA = pop(5*para.n+1 : 6*para.n);
R  = pop(6*para.n+1 : 7*para.n);
Cases = pop(7*para.n+1 : 8*para.n);

%NB: We only record recoveries to determine final epidemic size based off
%symptomatic infections only

% "social distancing" control: 70% decrease in contact rates if total infections
% above a given threshold
if SD(end,2) == 0
    if sum(IS) > para.Imax% || I(end) > para.Imax_risk
        factor = 0.3;
        SD(end+1,:) = [t,1];
    else
        factor = 1;
    end
else
    if sum(IS) < para.Imin% && I(end) < para.Imin_risk
        factor = 1;
        SD(end+1,:) = [t,0];
    else
        factor = 0.3;
    end
end

% ODE equations
dS = -factor.*S.*(para.beta*(IS + para.tau.*IA))./para.N;
dE1 = factor.*S.*(para.beta*(IS + para.tau.*IA))./para.N - 3*para.sigma.*E1;
dE2 = 3*para.sigma.*E1 - 3*para.sigma.*E2;
dE3 = 3*para.sigma.*E2 - 3*para.sigma.*E3;
dIS = 3*para.sigma.*para.da.*E3 - para.gamma.*IS;
dIA = 3*para.sigma.*(1-para.da).*E3 - para.gamma.*IA;
dR  = para.gamma.*(IS);
dCases = 3*para.sigma.*para.da.*E3;

dPop = [dS; dE1; dE2; dE3; dIS; dIA; dR; dCases];


end

end