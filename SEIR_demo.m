% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions
% maxtime: maximum time to simulate

function [Classes] = SEIR_demo(para,ICs,maxtime)

%Run ODE using ODE45
SD = [0,1];
opts = odeset('RelTol',1e-5,'MaxStep',0.1);
[t, pop] = ode45(@diff_SEIR_model, [0:0.1:maxtime], [ICs.S ICs.E1 ICs.E2 ICs.E3 ICs.I ICs.R], opts, para);

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E1',pop(:,para.n+1:2*para.n),'E2',pop(:,2*para.n+1:3*para.n), ...
            'E3',pop(:,3*para.n+1:4*para.n),'I',pop(:,4*para.n+1:5*para.n),'R',pop(:,5*para.n+1:6*para.n),'SD',floor(SD),'t',t);


%Diff equations

function dPop = diff_SEIR_model(t,pop,para)

S = pop(1 : para.n);
E1 = pop(para.n+1 : 2*para.n);
E2 = pop(2*para.n+1 : 3*para.n);
E3 = pop(3*para.n+1 : 4*para.n);
I = pop(4*para.n+1 : 5*para.n);
R = pop(5*para.n+1 : 6*para.n);

% "social distancing" control: 70% decrease in contact rates if total infections
% above a given threshold
if SD(end,2) == 0
    if sum(I) > para.Imax || I(end) > para.Imax_risk
        factor = 0.3;
        SD(end+1,:) = [t,1];
    else
        factor = 1;
    end
else
    if sum(I) < para.Imin && I(end) < para.Imin_risk
        factor = 1;
        SD(end+1,:) = [t,0];
    else
        if para.init == 0
            factor = 0.3;
        else
            factor = 1;
        end
    end
end

% ODE equations
dS = -factor.*S.*(para.beta*I)./para.N;
dE1 = factor.*S.*(para.beta*I)./para.N - 3*para.sigma.*E1;
dE2 = 3*para.sigma.*E1 - 3*para.sigma.*E2;
dE3 = 3*para.sigma.*E2 - 3*para.sigma.*E3;
dI = 3*para.sigma.*E3 - para.gamma.*I;
dR = para.gamma.*I;

dPop = [dS; dE1; dE2; dE3; dI; dR];


end

end