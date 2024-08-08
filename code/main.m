clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataFolder = "../processed_data/";
resultsFolder = "../results/";

% Outbreak labels (processed data files should be in the form
% label_processed.csv)
%outbreakLbl = ["covid_NZ_2020", "ebola_DRC_2018"];
outbreakLbl = ["ebola_DRC_2018"];


% Analyse each outbreak in turn
nOutbreaks = length(outbreakLbl);
for iOutbreak = 1:nOutbreaks

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % File name with raw data
    fNameData = sprintf('%s_processed.csv', outbreakLbl(iOutbreak));

    processed = readtable(dataFolder+fNameData);
    
    % Get model parameters and vector of times for which simulations will be
    % run
    [t, par] = getPar(outbreakLbl(iOutbreak));

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
    
    % Calculate pre-intervention R for each particle
    RpreInt = mean( Rt(:, t >= par.tRampStart-par.preIntWindow & t < par.tRampStart), 2);
    
    % Calculate gamma function using full data (hindsight)
    GTDF = [1, 1-cumsum(par.GTD)];      % Pad the survival function with a leading 1 as the first element applies to t  he previous days infections, which have all of their transmission potential remaining the next day

    Gamma = conv2(Yt, GTDF);
    Gamma = [zeros(par.nParticles, 1), Gamma(:, 1:length(t)-1)];          % offset by one so a case on day t contributes to future infection potential starting on day t+1
    
    
    % Calculate probability of ultimate extinction for an outbreak starting
    % with a single fully infectious seed case under each particle's
    % pre-intervention R
    PUE1 = calcPUE(RpreInt, par.k);
    
    
    % Probability of no futher infections given information up to time t
    pEnd = exp(-RpreInt.*GammaRT);
    pEndAvg = mean(pEnd);
    
    % Probability of ultimate extinction
    PUE = exp(-(1-PUE1).*RpreInt.*GammaRT);
    PUEAvg = mean(PUE);

    pNothing = exp(-RpreInt.*GammaRT) .* PhiRT;
    pNothingAvg = mean(pNothing);


    
    
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
    yyaxis left
    bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
    bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
    bh(2).FaceColor = [0 0 1]; 
    ylabel('reported daily cases')
    yyaxis right
    plot(t, PUEAvg, 'b-', t, pEndAvg, 'b--', t, pNothingAvg, 'b:' )
    ylabel('P(end of outbreak)')
    legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or reported cases", 'Location', 'northwest')
    xlim( [processed.t(1)-1, t(end)  ] )
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

end