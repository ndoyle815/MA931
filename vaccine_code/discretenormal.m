function P = discretenormal(k, mu, sigma)

% If the RV is Y this function returns P(Y=k) for a discrete analogue to
% the normal distribution whose continuous derivative is ~ N(mu, sigma)
% ks can be a float or array
% 
% NB: Y will map to the desired RV for vaccine arrival date T by 
% T = 180 + 60Y days

% k must be sorted if an array

n = length(k);
P = zeros(1,n);

lambda = exp(-(1-2*mu)/(2*sigma^2));
q = exp(-1/sigma^2);

% calculate denominator
jmin = mu-20*sigma;
jmax = mu+20*sigma;

denom = 0;
ki = 1;

for j = jmin:jmax
    Pj = (lambda^j)*q^(j*(j-1)/2);
    denom = denom + Pj;
    
    if j < k(1)
        P(1) = P(1) + Pj;
    elseif j > k(end)
        P(end) = P(end) + Pj;
    else
        P(ki) = P(ki) + Pj;
        ki = ki + 1;
    end

end

P = P/denom;

end