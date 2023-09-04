% function to perform preliminary simulation to obtain suitable ICs for
% main simulation

function [Prelim, Prelim_ICs] = Get_ICs(para)

nages = zeros(para.n,1);
% Define trivial initial conditions as a structure
ICs = struct('S',(1-para.E0).*para.N, 'E1',para.E0.*para.N, 'E2',nages, 'E3',nages, 'IA',nages, ...
             'IS',nages, 'IPH',nages, 'IH',nages, 'RA',nages, ...
             'RS',nages, 'Cases',nages, 'Hosp',nages, 'V',nages);

% Preliminary run
[Prelim] = ODEmodel(para,ICs);

% Obtain final conditions to begin main simulation
Prelim_ICs = struct('S',Prelim.S(end,:), 'E1',Prelim.E1(end,:), 'E2',Prelim.E2(end,:), 'E3',Prelim.E3(end,:), ...
                    'IA',Prelim.IA(end,:), 'IS',Prelim.IS(end,:), 'IPH',Prelim.IPH(end,:), 'IH',Prelim.IH(end,:), ...
                    'RA',Prelim.RA(end,:), 'RS',Prelim.RS(end,:), 'Cases',Prelim.Cases(end,:), ... 
                    'Hosp',Prelim.Hosp(end,:), 'V',Prelim.V(end,:));
