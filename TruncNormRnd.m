function r = TruncNormRnd(mu,sigma,l,u,N)

% draw from a truncated normal distribution

% INPUT:
% mu: mean of the distribution
% sigma: standard deviation
% l: lower limit
% u: upper limit
% N: number of values to be drawn 

% OUTPUT:
% r: value drawn from the truncated normal distribution 

t=truncate(makedist('normal',mu,sigma),l,u);

r=random(t,N,1);

end

