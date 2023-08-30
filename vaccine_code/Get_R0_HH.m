% function which uses the NextGeneration Matrix (NGM) approach to compute
% the basic reproduction number R_0 of the SEIR age model

function R = Get_R0_HH(para)

% matrix elements ordered x = (E_1a, E_2a, E_3a, I_a^A, I_a^S, I_a^PH, I_a^H)

T = zeros(7*para.n);
Sigma = zeros(7*para.n);

% transmission matrix T
T(1,3*para.n+1) = para.tau*para.beta(1,1);
T(1,3*para.n+2) = para.tau*para.beta(1,2);
T(1,3*para.n+3) = para.tau*para.beta(1,3);
T(1,4*para.n+1) = para.beta(1,1);
T(1,4*para.n+2) = para.beta(1,2);
T(1,4*para.n+3) = para.beta(1,3);
T(1,5*para.n+1) = para.beta(1,1);
T(1,5*para.n+2) = para.beta(1,2);
T(1,5*para.n+3) = para.beta(1,3);
T(1,6*para.n+1) = para.rho*para.beta(1,1);
T(1,6*para.n+2) = para.rho*para.beta(1,2);
T(1,6*para.n+3) = para.rho*para.beta(1,3);

T(2,3*para.n+1) = para.tau*para.beta(2,1);
T(2,3*para.n+2) = para.tau*para.beta(2,2);
T(2,3*para.n+3) = para.tau*para.beta(2,3);
T(2,4*para.n+1) = para.beta(2,1);
T(2,4*para.n+2) = para.beta(2,2);
T(2,4*para.n+3) = para.beta(2,3);
T(1,5*para.n+1) = para.beta(2,1);
T(1,5*para.n+2) = para.beta(2,2);
T(1,5*para.n+3) = para.beta(2,3);
T(1,6*para.n+1) = para.rho*para.beta(2,1);
T(1,6*para.n+2) = para.rho*para.beta(2,2);
T(1,6*para.n+3) = para.rho*para.beta(2,3);

T(3,3*para.n+1) = para.tau*para.beta(3,1);
T(3,3*para.n+2) = para.tau*para.beta(3,2);
T(3,3*para.n+3) = para.tau*para.beta(3,3);
T(3,4*para.n+1) = para.beta(3,1);
T(3,4*para.n+2) = para.beta(3,2);
T(3,4*para.n+3) = para.beta(3,3);
T(1,5*para.n+1) = para.beta(3,1);
T(1,5*para.n+2) = para.beta(3,2);
T(1,5*para.n+3) = para.beta(3,3);
T(1,6*para.n+1) = para.rho*para.beta(3,1);
T(1,6*para.n+2) = para.rho*para.beta(3,2);
T(1,6*para.n+3) = para.rho*para.beta(3,3);

% transition matrix Sigma
% E1,2,3
Sigma(1,1) = -para.n*para.sigma;
Sigma(2,2) = -para.n*para.sigma;
Sigma(3,3) = -para.n*para.sigma;
Sigma(1*para.n+1,1) = para.n*para.sigma;
Sigma(1*para.n+2,2) = para.n*para.sigma;
Sigma(1*para.n+3,3) = para.n*para.sigma;
Sigma(1*para.n+1,1*para.n+1) = -para.n*para.sigma;
Sigma(1*para.n+2,1*para.n+2) = -para.n*para.sigma;
Sigma(1*para.n+3,1*para.n+3) = -para.n*para.sigma;
Sigma(2*para.n+1,1*para.n+1) = para.n*para.sigma;
Sigma(2*para.n+2,1*para.n+2) = para.n*para.sigma;
Sigma(2*para.n+3,1*para.n+3) = para.n*para.sigma;
Sigma(2*para.n+1,2*para.n+1) = -para.n*para.sigma;
Sigma(2*para.n+2,2*para.n+2) = -para.n*para.sigma;
Sigma(2*para.n+3,2*para.n+3) = -para.n*para.sigma;

% IA
Sigma(3*para.n+1,2*para.n+1) = (1 - para.da(1))*para.n*para.sigma;
Sigma(3*para.n+2,2*para.n+2) = (1 - para.da(2))*para.n*para.sigma;
Sigma(3*para.n+3,2*para.n+3) = (1 - para.da(3))*para.n*para.sigma;
Sigma(3*para.n+1,3*para.n+1) = -para.gamma;
Sigma(3*para.n+2,3*para.n+2) = -para.gamma;
Sigma(3*para.n+3,3*para.n+3) = -para.gamma;

% IS
Sigma(4*para.n+1,2*para.n+1) = (1 - para.hosp_rates(1))*para.da(1)*para.n*para.sigma;
Sigma(4*para.n+2,2*para.n+2) = (1 - para.hosp_rates(2))*para.da(2)*para.n*para.sigma;
Sigma(4*para.n+3,2*para.n+3) = (1 - para.hosp_rates(3))*para.da(3)*para.n*para.sigma;
Sigma(4*para.n+1,4*para.n+1) = -para.gamma;
Sigma(4*para.n+2,4*para.n+2) = -para.gamma;
Sigma(4*para.n+3,4*para.n+3) = -para.gamma;

% IPH
Sigma(5*para.n+1,2*para.n+1) = para.hosp_rates(1)*para.da(1)*para.n*para.sigma;
Sigma(5*para.n+2,2*para.n+2) = para.hosp_rates(2)*para.da(2)*para.n*para.sigma;
Sigma(5*para.n+3,2*para.n+3) = para.hosp_rates(3)*para.da(3)*para.n*para.sigma;
Sigma(5*para.n+1,5*para.n+1) = -para.epsilon;
Sigma(5*para.n+2,5*para.n+2) = -para.epsilon;
Sigma(5*para.n+3,5*para.n+3) = -para.epsilon;

% IH
Sigma(6*para.n+1,5*para.n+1) = para.epsilon;
Sigma(6*para.n+2,5*para.n+2) = para.epsilon;
Sigma(6*para.n+3,5*para.n+3) = para.epsilon;
Sigma(6*para.n+1,6*para.n+1) = -para.delta;
Sigma(6*para.n+2,6*para.n+2) = -para.delta;
Sigma(6*para.n+3,6*para.n+3) = -para.delta;

% Next generation matrix
K = -T*inv(Sigma);
R = max(eig(K));