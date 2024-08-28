clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataFolder = "../processed_data/";
resultsFolder = "../results/";

% Outbreak labels (processed data files should be in the form label_processed.csv)
outbreakLbl = ["covid_NZ_2020", "ebola_DRC_2018"];

nScenarios = 8;

% Some computational settings:
ignoreDays = 30;        % for calculating time at which 95% is reached, ignore this many days at the start of the outbreak
pThresh = 0.95;         % threshold probability (for illustrative purposed)
nToPlot = 25;           % Number of particles to save for plotting
qt = [0.05 0.5 0.95];   % Quantiles of key outputs to save

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
    if ~isConsecutive(processed.t) 
        error('Data file needs to have a field called t consisting of consecutive dates');
    end
    
    for iScenario = 1:nScenarios

        % Get model parameters and vector of times for which simulations will be
        % run
        [t, par] = getPar(outbreakLbl(iOutbreak), iScenario);

        % Copy imported data onto the time array returned by getPar
        nCasesLoc = zeros(size(t));
        nCasesImp = zeros(size(t));
        nCasesLoc(ismember(t, processed.t)) = processed.nCasesLoc;
        nCasesImp(ismember(t, processed.t)) = processed.nCasesImp;
 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Run particle filter
        [Rt, It, Yt, Zt, Ct, It_imp, Zt_imp, GammaRT, PhiRT, RpreInt, ESS, LL] = runPF(t, nCasesLoc, nCasesImp, par);
        
        % Do post-processg in particle filter results to calculate
        % pre-intervention R and different end-of-outbreak probabilities for
        % each particle
        [PUE, pNoInf, pNoInfOrCases] = postProcess(t, Rt, GammaRT, PhiRT, RpreInt, par);
        
        results.outbreak(iRow) = outbreakLbl(iOutbreak);
        results.iScenario(iRow) = iScenario;
        results.t{iRow} = t;
        results.par{iRow} = par;
        results.Rt_quantiles{iRow} = quantile(Rt, qt);
        results.It_quantiles{iRow} = quantile(It, qt);
        results.Zt_quantiles{iRow} = quantile(Zt, qt);
        results.It_imp_quantiles{iRow} = quantile(It_imp, qt);
        results.Zt_imp_quantiles{iRow} = quantile(Zt_imp, qt);
        results.Ct_quantiles{iRow} = quantile(Ct, qt);
        % Just save selection of particles instead of the full set
        results.Rt{iRow} = Rt(1:nToPlot, :);
        results.It{iRow} = It(1:nToPlot, :);
        results.Zt{iRow} = Zt(1:nToPlot, :);
        results.Ct{iRow} = Ct(1:nToPlot, :);
        % Just save mean, SD, and histogram counts for RpreInt
        results.RpreInt_mean{iRow} = mean(RpreInt);
        results.RpreInt_sd{iRow} = std(RpreInt);
        [edges, freq] = histcounts(RpreInt);
        results.RpreInt_edges{iRow} = edges;
        results.RpreInt_freq{iRow} = freq;
        results.PUE{iRow} = PUE;
        results.pNoInf{iRow} = pNoInf;
        results.pNoInfOrCases{iRow} = pNoInfOrCases;
        % Save time at which elimination probabilities first reach 95%:
        results.tPUE95{iRow} = t(find(1:length(t) >= ignoreDays & PUE >= pThresh, 1, "first"));  
        results.tpNoInf95{iRow} = t(find(1:length(t) >= ignoreDays & pNoInf >= pThresh, 1, "first"));
        results.tpNoInfOrCases95{iRow} = t(find(1:length(t) >= ignoreDays & pNoInfOrCases >= pThresh, 1, "first"));
        iRow = iRow+1;  
    end    
end

save(resultsFolder+"results.mat");
