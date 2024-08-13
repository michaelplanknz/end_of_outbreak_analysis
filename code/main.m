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

kValues = [inf, 1, 0.2];
nk = length(kValues);

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
    


    for ik = 1:nk

        % Get model parameters and vector of times for which simulations will be
        % run
        [t, par] = getPar(outbreakLbl(iOutbreak));

        par.k = kValues(ik);

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
        
        results(iRow).outbreak = outbreakLbl;
        results(iRow).par = par;
        results(iRow).t = t;
        results(iRow).RpreInt = RpreInt;
        results(iRow).PUE = PUE;
        results(iRow).pNoInf = pNoInf;
        results(iRow).pNoInfOrCases = pNoInfOrCases;
        iRow = iRow+1;
    
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nPlot = 25;     % Number of particles to plot
    grey = [0.7 0.7 0.7];
    
    
    figure;
    tiledlayout(2, 2);
    
    nexttile;
    plot(t, Rt(1:nPlot, :), 'Color', grey)
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('reproduction number')
    
    nexttile;
    plot(t, It(1:nPlot, :), 'Color', grey)
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('daily infections')
    
    nexttile;
    plot(t, par.pReport*Zt(1:nPlot, :), 'Color', grey)
    hold on
    plot(t, nCasesLoc )
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('expected daily cases')
    
    nexttile;
    plot(t, Ct(1:nPlot, :), 'Color', grey)
    hold on
    plot(t, nCasesLoc )
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('simulated daily cases')

    figure;
    histogram(RpreInt)
    ylabel('pre-intervention reproduction number')
    

    iMinPlot = 15;
    figure;
    yyaxis left
    bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
    bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
    bh(2).FaceColor = [0 0 1]; 
    ylabel('reported daily cases')
    yyaxis right
    plot(t(iMinPlot:end), PUE(iMinPlot:end), 'b-', t(iMinPlot:end), pNoInf(iMinPlot:end), 'b--', t(iMinPlot:end), pNoInfOrCases(iMinPlot:end), 'b:' )
    xline(par.tRampStart, 'k:');
    ylabel('P(end of outbreak)')
    legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or reported cases", 'Location', 'northwest')
    xlim( [processed.t(1)-1, t(end)  ] )
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

end