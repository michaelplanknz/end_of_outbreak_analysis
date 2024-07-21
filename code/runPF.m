function [Rt, It, Yt, Zt, Ct, GammaRT, ESS, LL] = runPF(t, nCasesLoc, nInfImp, par)

% Function to run the particle filter on given input data
%
% USAGE: [Rt, It, Yt, Zt, Ct, GammaRT, LL, ESS] = runPF(t, nCasesLoc, nInfImp, par)
%
% INPUTS: t - vector of daily dates defining the model simulation period
%         nCasesLoc - corresponding vector of daily local cases
%         nInfImp - correspodnig time series of daily imported infecitons
%         (seed cases)
%         par - structure of model parameters with fields
%         - par.nParticles - number of particles to use
%         - par.GTD - vector of probability masses for generation time t = 1, 2, .. days
%         - par.RTD - vector of probability masses for infection to reporting time t = 0, 1, 2, .. days
%         - par.R_shape, par.R_scale - shape and scale parameters for the
%         prior distribution for initial Rt
%         - par.k - dispersion parameter for individual transmissibility (inf for Poisson offspring)
%         - par.deltat, par.sigmat - vector (of same length as t) of means and s.d. for the daily change in Rt 
%         - par.obsModel - noise model for observed daily cases data, either
%         "bin" or "negbin" 
%         - par.kObs (if obsModel is "negbin") - dispersion parameter for
%         observed daily cases (inf for Poisson dist)
%         - par.pReport - probability of infecitons being reported as cases
%         - par.resampleLag - fixed lag (number of time steps) for bootstrap particle filter resampling
%
% OUTPUTS: Rt - matrix of reproduction numbers - (i,j) element corresponds to particle i on day j
%          It - matrix of daily infections  - (i,j) element corresponds to particle i on day j
%          Yt - matrix of transmissibility of people infected on a given day - (i,j) element corresponds to particle i on day j
%          Zt - matrix of infections by assigned date of report (independent of whether they actually reported as a case or not)  - (i,j) element corresponds to particle i on day j
%          Ct - matrix of simulated daily cases  - (i,j) element corresponds to particle i on day j
%          ESS - vector of the number of unique particles at each time step
%          LL - particle marginal log likelihood - can be used in PMMH to fit fixed parameters


GTDF = [1, 1-cumsum(par.GTD)];      % Pad the survival function with a leading 1 as the first element applies to t  he previous days infections, which have all of their transmission potential remaining the next day

nSteps = length(t);
iLastData = length(nCasesLoc);

% Initialise variables for renewal equation particle filter
Rt = zeros(par.nParticles, nSteps);
It = zeros(par.nParticles, nSteps);
Yt = zeros(par.nParticles, nSteps);
GammaRT = zeros(par.nParticles, nSteps);
Zt = zeros(par.nParticles, nSteps);
ESS = par.nParticles*ones(nSteps, 1);
lmw = zeros(1, nSteps);


% For the inital value of Rt (needed for first time step), draw from priors 
Rt(:, 1) = gamrnd(par.R_shape, par.R_scale, par.nParticles, 1);

if isfinite(par.k)    
   Yt(:, 1) = gamrnd( par.k * nInfImp(1), 1/par.k );
else
   Yt(:, 1) = nInfImp(1);
end

% Loop through time steps
for iStep = 2:nSteps
   Rt(:, iStep) = max(0, Rt(:, iStep-1) + par.deltat(iStep) + par.sigmat(iStep)*randn(par.nParticles, 1));        

   ind = iStep-1:-1:max(1, iStep-length(par.GTD));
   It(:, iStep) = poissrnd(  Rt(:, iStep) .* sum(par.GTD(1:length(ind)).*Yt(:, ind), 2 ) );       % renewal equation
                         
   if isfinite(par.k)    
       Yt(:, iStep) = gamrnd( par.k * (It(:, iStep) + nInfImp(iStep)), 1/par.k );
   else
       Yt(:, iStep) = It(:, iStep) + nInfImp(iStep);
   end
   GammaRT(:, iStep) = sum(GTDF(1:length(ind)).*Yt(:, ind), 2 );
   ind = iStep:-1:max(1, iStep+1-length(par.RTD));
   Zt(:, iStep) = sum(par.RTD(1:length(ind)) .* (It(:, ind)), 2);             % infections by date of report



   % Particle resampling (only during period for which data is available)
   if iStep <= iLastData
        if par.obsModel == "negbin" & isfinite(par.kObs)   
           weights = nbinpdf(nCasesLoc(iStep), par.kObs, par.kObs./(par.pReport*Zt(:, iStep)+par.kObs));
        elseif par.obsModel == "negbin" & ~isfinite(par.kObs)
           weights = poisspdf(nCasesLoc(iStep), par.pReport*Zt(:, iStep) );
        elseif par.obsModel =="bin"
            weights = binopdf(nCasesLoc(iStep), Zt(:, iStep), par.pReport );
        else
            error(sprintf('Invalid observation model: %s', par.obsModel));
        end
        lmw(iStep) = log(mean(weights));
        
        resampInd = randsample(par.nParticles, par.nParticles, true, weights);

        ESS(iStep) = length(unique(resampInd));

        % Resample particles according to weights
        % NB GammaRT is not resampled because it represents the future
        % force of infection estimated with information known up to a given
        % point in time 
        iResample = max(1, iStep-par.resampleLag);
        Rt(:, iResample:end) = Rt(resampInd, iResample:end);
        It(:, iResample:end) = It(resampInd, iResample:end);
        Yt(:, iResample:end) = Yt(resampInd, iResample:end);
        Zt(:, iResample:end) = Zt(resampInd, iResample:end);
   end

end

LL = sum(lmw);

 % Generate samples from the reported case distribution to construct
 % prediction intervals:
if par.obsModel == "negbin" & isfinite(par.kObs)  
    Ct = nbinrnd(par.kObs, par.kObs./(par.pReport*Zt+par.kObs));     
elseif par.obsModel == "negbin" & ~isfinite(par.kObs)  
    Ct = poissrnd(par.pReport*Zt);
elseif par.obsModel == "bin"
    Ct = binornd(Zt, par.pReport );
end




