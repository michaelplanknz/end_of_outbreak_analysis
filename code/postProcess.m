function [PUE, pNoInf, pNoInfOrCases] = postProcess(GammaRT, PhiRT, RpreInt, par)

% Calculate probability of ultimate extinction for an outbreak starting
% with a single fully infectious seed case under each particle's
% pre-intervention R
PUE1 = calcPUE(RpreInt, par.k);


% Probability of ultimate extinction given information up to time t
PUE = exp(-(1-PUE1).*RpreInt.*GammaRT);

% Probability of no futher infections given information up to time t
pNoInf = exp(-RpreInt.*GammaRT);

% Probability of no futher infections and no future reported cases given information up to time t
pNoInfOrCases = exp(-RpreInt.*GammaRT) .* PhiRT;

