function [PUE, pNoInf, pNoInfOrCases] = postProcess(t, Rt, GammaRT, PhiRT, RpreInt, par)

% Calculate probability of ultimate extinction for an outbreak starting
% with a single fully infectious seed case under each particle's
% pre-intervention R
PUE1 = calcPUE(RpreInt, par.k);


% Probability of ultimate extinction given information up to time t
PUE_i = exp(-(1-PUE1).*RpreInt.*GammaRT);
PUE = mean(PUE_i);

% Probability of no futher infections given information up to time t
pNoInf_i = exp(-RpreInt.*GammaRT);
pNoInf = mean(pNoInf_i);

% Probability of no futher infections and no future reported cases given information up to time t
pNoInfOrCases_i = exp(-RpreInt.*GammaRT) .* PhiRT;
pNoInfOrCases = mean(pNoInfOrCases_i);

