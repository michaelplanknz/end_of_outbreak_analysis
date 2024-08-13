clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataFolder = "../processed_data/";
resultsFolder = "../results/";

% Outbreak labels (processed data files should be in the form
% label_processed.csv)
outbreakLbl = ["covid_NZ_2020", "ebola_DRC_2018"];

nToPlot = 25;     % Number of particles to plot

% Analyse each outbreak in turn
nOutbreaks = length(outbreakLbl);
iRow = 1;
for iOutbreak = 1:nOutbreaks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % File name with raw data
    fNameData = sprintf('%s_processed.csv', outbreakLbl(iOutbreak));

    processed = readtable(dataFolder+fNameData);
    
    scenarioTab{iOutbreak} = getScenarios(outbreakLbl(iOutbreak));
    nScenarios = height(scenarioTab{iOutbreak});


    for iScenario = 1:nScenarios

        % Get model parameters and vector of times for which simulations will be
        % run
        [t, par] = getPar(outbreakLbl(iOutbreak), scenarioTab{iOutbreak}(iScenario, :) );

        % Copy imported data onto the time array returned by getPar
        nCasesLoc = zeros(size(t));
        nCasesImp = zeros(size(t));
        nCasesLoc(ismember(t, processed.t)) = processed.nCasesLoc;
        nCasesImp(ismember(t, processed.t)) = processed.nCasesImp;
    
    
        % Estimating imported infection dates and accounting for
        % unreported imported infections
        nInfImp = [nCasesImp(1+par.tInfImp:end), zeros(1, par.tInfImp)];       % imported infections assumed to occur par.tInfImp days before reported
        nInfImp(t > par.tMIQ) = 0;                                                  % ignore imported cases with assigned infection date after introduction of MIQ
        nUndetTot = round((1-par.pReport)/par.pReport*sum(nInfImp));                       % number of undetected imported cases, under assumed rpeorting probability
        tUndet = randsample(length(nInfImp), nUndetTot, true, nInfImp );            % sample time (index) of undetected cases by reampling with replacement from detected cases
        nUndet = histcounts(tUndet, 1:length(nInfImp)+1);                             % number of undetected cases on each day
        nInfImp = nInfImp + nUndet;                                                 % add undeteced cases to data time series
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Run particle filter
        [Rt, It, Yt, Zt, Ct, GammaRT, PhiRT, ESS, LL] = runPF(t, nCasesLoc, nInfImp, par);
        
        % Do post-processg in particle filter results to calculate
        % pre-intervention R and different end-of-outbreak probabilities for
        % each particle
        [RpreInt, PUE, pNoInf, pNoInfOrCases] = postProcess(t, Rt, GammaRT, PhiRT, par);
        
        results{iRow}.outbreak = outbreakLbl(iOutbreak);
        results{iRow}.iScenario = iScenario;
        results{iRow}.t = t;
        results{iRow}.par = par;
        results{iRow}.Rt_mean = mean(Rt);
        results{iRow}.It_mean = mean(It);
        results{iRow}.Zt_mean = mean(Zt);
        results{iRow}.Ct_mean = mean(Ct);
        results{iRow}.Rt = Rt(1:nToPlot, :);
        results{iRow}.It = It(1:nToPlot, :);
        results{iRow}.Zt = Zt(1:nToPlot, :);
        results{iRow}.Ct = Ct(1:nToPlot, :);
        [freq, edges] = histcounts(RpreInt);
        results{iRow}.RpreInt_edges = edges;
        results{iRow}.RpreInt_freq = freq;
        results{iRow}.PUE = PUE;
        results{iRow}.pNoInf = pNoInf;
        results{iRow}.pNoInfOrCases = pNoInfOrCases;
        iRow = iRow+1;  
    
    end    
end

save(resultsFolder+"results.mat");

