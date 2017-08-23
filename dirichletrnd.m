 function r = drchrnd(a)
% take a sample from a dirichlet distribution
p = length(a);
r = gamrnd(a,1,1,p);
r = r ./ repmat(sum(r,2),1,p);