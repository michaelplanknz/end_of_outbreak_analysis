function f = discDist(pdfFunc, a, b )

% Calculate a discrete probability mass function from a continuous PDF on
% non-negative values using the Method of Cori et al
% The discrete PMF will be truncated (at the maximum value b) and
% normalised to sum to 1
% If the starting value a is 1 instead of 0, the first element (probability of 1) will
% include all mass between 0 and 1, plus a share of mass between 1 and 2

nn = a:b;

nValues = length(nn);
f = zeros(1, nValues);

for iValue = 1:nValues
    ig1 = @(x)(pdfFunc(x).*(1+nn(iValue)-x ) );
    ig2 = @(x)(pdfFunc(x).*(1-nn(iValue)+x ) );
    f(iValue) = quad(ig1, nn(iValue)-1, nn(iValue)) + quad(ig2, nn(iValue), nn(iValue)+1);
end
f = f/sum(f);




