function [sh, sc] = estimateR(nCasesImp, nCasesLoc, GTD, relInfImp )

% Estimate shape and scale parameters for the gamma posterior distribution
% for R assuming a uniform prior, from time series data for the number of
% imported cases nCasesImp and local cases nCasesLoc.
% GTD is a vector contain the Pprobability mass functio nfor the GT
% distribution on 1, 2, 3, ... (days
% relInfImp is a parameter in [0,1] specifying the assumed trasmission
% portnetial of imported cases relative to local cases

nCasesTot = nCasesLoc + relInfImp*nCasesImp;
if length(GTD) <= length(nCasesTot)-1
    F = ones(1, length(nCasesTot)-1);
    F(1:length(GTD)) = cumsum(GTD);
else
    F = cumsum(GTD(1:length(nCasesTot)-1));
end

sh = 1+sum(nCasesLoc);
sc = 1./sum( F .* nCasesTot(end-1:-1:1) ); 
