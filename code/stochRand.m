function n = stochRand(x)

% Stochastic rounding of x to adjacent integer
% Each element of the input x will be round to either floor(x) or ceil(x)
% with probability depending on how close x is to these two neighbouring
% integers

n0 = floor(x);
y = mod(x, 1);
U = rand(size(x));
n = n0 + (U < y);

