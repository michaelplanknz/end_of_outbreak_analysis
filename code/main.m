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
        particles = runPF(t, nCasesLoc, nCasesImp, par);
        
        % Store results for saving in a structrue:
        results.outbreak(iRow) = outbreakLbl(iOutbreak);
        results.iScenario(iRow) = iScenario;
        results.t{iRow} = t;
        results.par{iRow} = par;
        results.Rt_quantiles{iRow} = quantile(particles.Rt, qt);
        results.It_quantiles{iRow} = quantile(particles.It, qt);
        results.Zt_quantiles{iRow} = quantile(particles.Zt, qt);
        results.Ct_quantiles{iRow} = quantile(particles.Ct, qt);
        results.It_imp_quantiles{iRow} = quantile(particles.It_imp, qt);
        results.Zt_imp_quantiles{iRow} = quantile(particles.Zt_imp, qt);
        results.Ct_imp_quantiles{iRow} = quantile(particles.Ct_imp, qt);
        % Just save mean, SD, and histogram counts for RpreInt
        results.RpreInt_mean{iRow} = mean(particles.RpreInt);
        results.RpreInt_sd{iRow} = std(particles.RpreInt);
        [edges, freq] = histcounts(particles.RpreInt);
        results.RpreInt_edges{iRow} = edges;
        results.RpreInt_freq{iRow} = freq;
        % Save mean and quantiles for end-of-outbreak probabilities
        results.PUE_mean{iRow} = mean(particles.PUE);
        results.pNoInf_mean{iRow} = mean(particles.pNoInf);
        results.pNoInfOrCases_mean{iRow} = mean(particles.pNoInfOrCases);
        results.PUE_quantiles{iRow} = quantile(particles.PUE, qt);
        results.pNoInf_quantiles{iRow} = quantile(particles.pNoInf, qt);
        results.pNoInfOrCases_quantiles{iRow} = quantile(particles.pNoInfOrCases, qt);
        % Save time at which mean elimination probabilities first reach 95%:
        results.tPUE95{iRow} = t(find(1:length(t) >= ignoreDays & results.PUE_mean{iRow} >= pThresh, 1, "first"));  
        results.tpNoInf95{iRow} = t(find(1:length(t) >= ignoreDays & results.pNoInf_mean{iRow} >= pThresh, 1, "first"));
        results.tpNoInfOrCases95{iRow} = t(find(1:length(t) >= ignoreDays & results.pNoInfOrCases_mean{iRow} >= pThresh, 1, "first"));
        iRow = iRow+1;  
    end    
end

save(resultsFolder+"results.mat");
