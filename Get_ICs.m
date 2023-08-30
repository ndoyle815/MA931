% function to perform preliminary simulation to obtain suitable ICs for
% main simulation

function [Prelim, Prelim_ICs] = Get_ICs(para)

% Define trivial initial conditions as a structure
E0 = 0.0025;  % initial exposed
ICs = struct('S',(1-E0).*para.N, 'E1',E0.*para.N, 'E2',zeros(para.n,1), 'E3',zeros(para.n,1), ...
             'IS',zeros(para.n,1), 'IA',zeros(para.n,1), 'R',zeros(para.n,1), 'Cases',zeros(para.n,1));

% Preliminary run
[Prelim] = SEIR_demo(para,ICs);

% Obtain final conditions to begin main simulation
Prelim_ICs = struct('S',Prelim.S(end,:), 'E1',Prelim.E1(end,:), 'E2',Prelim.E2(end,:), ...
                    'E3',Prelim.E3(end,:), 'IS',Prelim.IS(end,:), 'IA',Prelim.IA(end,:), ...
                    'R',Prelim.R(end,:),'Cases',Prelim.Cases(end,:));