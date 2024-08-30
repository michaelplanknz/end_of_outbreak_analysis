function particles = runPF(t, nCasesLoc, nCasesImp, par)

% Function to run the particle filter on given input data
%
% USAGE: particles = runPF(t, nCasesLoc, nInfImp, par)
%
% INPUTS: t - vector of daily dates defining the model simulation period
%         nCasesLoc - corresponding vector of daily local cases
%         nCasesImp - correspodnig time series of daily imported cases
%         par - structure of model parameters with fields
%         - par.nParticles - number of particles to use
%         - par.GTD - vector of probability masses for generation time t = 1, 2, .. days
%         - par.RTD - vector of probability masses for infection to reporting time t = 0, 1, 2, .. days
%         - par.R_shape, par.R_scale - shape and scale parameters for the
%         prior distribution for initial Rt
%         - par.k - dispersion parameter for individual transmissibility (inf for Poisson offspring)
%         - par.deltat, par.sigmat - vectors (of same length as t) of means and s.d. for the daily change in Rt 
%         - par.obsModel - noise model for observed daily cases data, either
%         "bin" or "negbin" 
%         - par.kObs (if obsModel is "negbin") - dispersion parameter for
%         observed daily cases (inf for Poisson dist)
%         - par.pReport - probability of infecitons being reported as cases
%         - par.resampleLag - fixed lag (number of time steps) for bootstrap particle filter resampling
%
% OUTPUTS: particles - a structure with the following fields
%          Rt - matrix of reproduction numbers - (i,j) element corresponds to particle i on day j
%          It - matrix of daily infections  - (i,j) element corresponds to particle i on day j
%          Yt - matrix of transmissibility of people infected on a given day - (i,j) element corresponds to particle i on day j
%          Zt - matrix of infections by assigned date of report (independent of whether they actually reported as a case or not)  - (i,j) element corresponds to particle i on day j
%          Ct - matrix of simulated daily cases  - (i,j) element corresponds to particle i on day j
%          GammaRT - matrix of gamma(t) values calculating in real tme (i.e. only using data availab le up to time t)
%          PhiRT - matrix of phi(t) values representing the probability that there will be no further cases reported on or after day t arising form infections that occurred prior to day t
%          LL - particle marginal log likelihood - can be used in PMMH to fit fixed parameters


GTDF = [1, 1-cumsum(par.GTD)];      % Pad the survival function with a leading 1 as the first element applies to t  he previous days infections, which have all of their transmission potential remaining the next day
RTDC = cumsum(par.RTD);


nSteps = length(t);
iLastData = length(nCasesLoc);



% Estimating imported infection dates and accounting for
% unreported imported infections
if par.filterImportsFlag
    c = conv(nCasesImp, fliplr(par.RTD));
    if c(1:length(par.RTD)-1) > 0
        fprintf('Warning: %.3f imported cases have been lost by the backward convolution process')
    end
    nImpPrior = c(length(par.RTD):end)/par.pReport;
    It_imp = zeros(par.nParticles, nSteps);
    It_imp(:, 1) = poissrnd(nImpPrior(1), par.nParticles, 1);
    figure(200);
    plot(t, nCasesImp, t, par.pReport*nImpPrior)
    drawnow
else
    if sum(nCasesImp(1:par.tInfImp)) > 0
        fprintf('Warning: %i imported cases are being ignored because they occured during the first %i days of the time vector passed to runPF()\n', sum(nCasesImp(1:par.tInfImp)), par.tInfImp)
    end
    nImpDet = [nCasesImp(1+par.tInfImp:end), zeros(1, par.tInfImp)];       % imported infections assumed to occur par.tInfImp days before reported
    nUndetTot = round((1-par.pReport)/par.pReport*sum(nImpDet));                       % number of undetected imported cases, under assumed rpeorting probability
    tUndet = randsample(length(nImpDet), nUndetTot, true, nImpDet );            % sample time (index) of undetected cases by reampling with replacement from detected cases
    nImpUndet = histcounts(tUndet, 1:length(nImpDet)+1);                             % number of undetected cases on each day
    It_imp = repmat(nImpDet + nImpUndet, par.nParticles, 1);                                                 % add undeteced cases to data time series
end
                                             

% Initialise variables for renewal equation particle filter
Rt = zeros(par.nParticles, nSteps);
It = zeros(par.nParticles, nSteps);
Yt = zeros(par.nParticles, nSteps);
GammaRT = zeros(par.nParticles, nSteps);
PhiRT = zeros(par.nParticles, nSteps);
Zt = zeros(par.nParticles, nSteps);

Zt_imp = zeros(par.nParticles, nSteps);
lmw = zeros(1, nSteps);


% For the inital value of Rt (needed for first time step), draw from priors 
Rt(:, 1) = gamrnd(par.R_shape, par.R_scale, par.nParticles, 1);

% Loop through time steps
for iStep = 2:nSteps
    % When you get to intevention start time, save the filtering distribution for Rt in the period leading up to intervention 
    if t(iStep) == par.tRampStart
        RpreInt = mean(Rt(:, iStep-par.preIntWindow:iStep-1), 2);    
    end

    % Update reproduction number Rt
  % Rt(:, iStep) = max(0, Rt(:, iStep-1) + par.deltat(iStep) + par.sigmat(iStep)*randn(par.nParticles, 1));        
   Rt(:, iStep) = Rt(:, iStep-1) .* lognrnd(par.deltat(iStep), par.sigmat(iStep), par.nParticles, 1);

   if par.filterImportsFlag
       It_imp(:, iStep) = poissrnd(nImpPrior(iStep), par.nParticles, 1);        % sample imported infections from the prior
   end

   % Set daily aggregate infectivity (Yt)
   if isfinite(par.k)    
       Yt(:, iStep-1) = gamrnd( par.k * (It(:, iStep-1) + par.relInfImp * ~(t(iStep-1) > par.tMIQ ) *It_imp(:, iStep-1)), 1/par.k );
   else
       Yt(:, iStep-1) = It(:, iStep-1) + par.relInfImp * ~(t(iStep-1) > par.tMIQ ) * It_imp(:, iStep-1);
   end

   % Compute It according to renewal equation
   ind = iStep-1:-1:max(1, iStep-length(par.GTD));
   It(:, iStep) = poissrnd(  Rt(:, iStep) .* sum(par.GTD(1:length(ind)).*Yt(:, ind), 2 ) );     

   % Calculate gamma_t (using only current and previous information - no
   % subsequent resampling)
   GammaRT(:, iStep) = sum(GTDF(1:length(ind)).*Yt(:, ind), 2 );

   % Assign infections a notificatoin date and add them to the vecotr of
   % future case notifications time series (Zt)
   ind = iStep:min(iStep+length(par.RTD)-1, nSteps);
   futureCases = mnrnd( It(:, iStep), par.RTD  );
   Zt(:, ind) = Zt(:, ind) + futureCases(:, 1:length(ind));
   
   % Do the same for imported infections
   futureCasesImp = mnrnd( It_imp(:, iStep), par.RTD  );
   Zt_imp(:, ind) = Zt_imp(:, ind) + futureCasesImp(:, 1:length(ind));

   % Particle resampling (only during period for which data is available)
   if iStep <= iLastData
        if par.obsModel == "negbin" & isfinite(par.kObs)   
           if par.filterImportsFlag     
                weightsImp = nbinpdf(nCasesImp(iStep), par.kObs, par.kObs./(par.pReport*Zt_imp(:, iStep)+par.kObs)); % Additional weights for hte likelihood of imported case data
           else
               weightsImp = 1;
           end
            weights = weightsImp .* nbinpdf(nCasesLoc(iStep), par.kObs, par.kObs./(par.pReport*Zt(:, iStep)+par.kObs));
           PhiRT(:, iStep) = prod( nbinpdf(0, par.kObs, par.kObs./(par.pReport * Zt(:, iStep:end) + par.kObs )) , 2);                             % Calculate phi_t = Pr(no future reported cases from existing local infections)
        elseif par.obsModel == "negbin" & ~isfinite(par.kObs)
           if par.filterImportsFlag
                weightsImp = poisspdf(nCasesImp(iStep), par.pReport*Zt_imp(:, iStep));
           else
               weightsImp = 1;
           end
           weights = weightsImp .* poisspdf(nCasesLoc(iStep), par.pReport*Zt(:, iStep) );
           PhiRT(:, iStep) = poisspdf(0, par.pReport * sum(Zt(:, iStep:end), 2 ) );                          
        elseif par.obsModel =="bin"
           if par.filterImportsFlag
                weightsImp = binopdf(nCasesImp(iStep), Zt_imp(:, iStep), par.pReport);
           else
               weightsImp = 1;
           end            
           weights = weightsImp .* binopdf(nCasesLoc(iStep), Zt(:, iStep), par.pReport );
           if max(weights) == 0
                fprintf('Warning: at time step %i/%i, max Zt = %i and local cases = %i, max Zt_imp = %i and imported cases = %i, using particles with minimal SSE\n', iStep, nSteps, max(Zt(:, iStep)), nCasesLoc(iStep), max(Zt_imp(:, iStep)), nCasesImp(iStep) )
                if par.filterImportsFlag                    
                    distImp = nCasesImp(iStep) - Zt_imp(:, iStep).^2;
                else
                    distImp = 0;
                end
                dist = distImp + nCasesLoc(iStep) - Zt(:, iStep).^2;
                weights = (dist == min(dist));

           end
           PhiRT(:, iStep) = binopdf(0, sum(Zt(:, iStep:end), 2 ), par.pReport );
        else
            error(sprintf('Invalid observation model: %s', par.obsModel));
        end
        lmw(iStep) = log(mean(weights));
        

        try       
            resampInd = randsample(par.nParticles, par.nParticles, true, weights);
        catch
            stop
        end
        

        % Resample particles according to weights
        % NB GammaRT is not resampled because it represents the future
        % force of infection estimated with information known up to a given
        % point in time. Similarly for PhiRT 
        iResample = max(1, iStep-par.resampleLag);
        Rt(:, iResample:end) = Rt(resampInd, iResample:end);
        It(:, iResample:end) = It(resampInd, iResample:end);
        Yt(:, iResample:end) = Yt(resampInd, iResample:end);
        Zt(:, iResample:end) = Zt(resampInd, iResample:end);
        It_imp(:, iResample:end) = It_imp(resampInd, iResample:end);
        Zt_imp(:, iResample:end) = Zt_imp(resampInd, iResample:end);
   end

end


 % Generate samples from the reported case distribution to construct
 % prediction intervals:
if par.obsModel == "negbin" & isfinite(par.kObs)  
    Ct = nbinrnd(par.kObs, par.kObs./(par.pReport*Zt+par.kObs)); 
    Ct_imp = nbinrnd(par.kObs, par.kObs./(par.pReport*Zt_imp+par.kObs));  
elseif par.obsModel == "negbin" & ~isfinite(par.kObs)  
    Ct = poissrnd(par.pReport*Zt);
    Ct_imp = poissrnd(par.pReport*Zt_imp);
elseif par.obsModel == "bin"
    Ct = binornd(Zt, par.pReport );
    Ct_imp = binornd(Zt_imp, par.pReport );
end




% Store outputs in a structure called particles:
particles.Rt = Rt;
particles.It = It;
particles.Yt = Yt;
particles.Zt = Zt;
particles.Ct = Ct;
particles.It_imp = It_imp;
particles.Zt_imp = Zt_imp;
particles.Ct_imp = Ct_imp;
particles.RpreInt = RpreInt;
particles.LL = sum(lmw);

% Calculate the different end-of-outbreak probabilities:
[particles.PUE, particles.pNoInf, particles.pNoInfOrCases] = postProcess(GammaRT, PhiRT, RpreInt, par);

