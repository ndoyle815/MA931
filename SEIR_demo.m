% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions
% maxtime: maximum time to simulate
% n: number of age categories

function [Classes] = SEIR_demo(para,ICs,maxtime)

%Run ODE using ODE45
SD = [0,1];
opts = odeset('RelTol',1e-5,'MaxStep',0.1);
[t, pop] = ode45(@diff_SEIR_model, [0:0.1:maxtime], [ICs.S ICs.E ICs.I ICs.R], opts, para);

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E',pop(:,para.n+1:2*para.n),'I',pop(:,2*para.n+1:3*para.n), ...
                 'R',pop(:,3*para.n+1:4*para.n),'SD',SD,'t',t);


%Diff equations

function dPop = diff_SEIR_model(t,pop,para)

S = pop(1 : para.n);
E = pop(para.n+1 : 2*para.n);
I = pop(2*para.n+1 : 3*para.n);
R = pop(3*para.n+1 : 4*para.n);

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
        %t
        %sum(I)
        factor = 1;
        SD(end+1,:) = [t,0];
    else
        factor = 0.3;
    end
end

% ODE equations
dS = -factor.*S.*(para.beta*I)./para.N;
dE = factor.*S.*(para.beta*I)./para.N - para.sigma.*E;
dI = para.sigma.*E - para.gamma.*I;
dR = para.gamma.*I;

dPop = [dS; dE; dI; dR];


end

end