clear 
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataFolder = "../processed_data/";
resultsFolder = "../results/";

load(resultsFolder+"results.mat");
    

baseScenario = 6;
grey = [0.7 0.7 0.7];

    
for iOutbreak = 1:nOutbreaks
    fNameData = sprintf('%s_processed.csv', outbreakLbl(iOutbreak));
    processed = readtable(dataFolder+fNameData);

    iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.iScenario == baseScenario);

    t = results.t{iRow};
    par = results.par{iRow};
    Rt = results.Rt{iRow};
    It = results.It{iRow};
    Zt = results.Zt{iRow};
    Ct = results.Ct{iRow};
    RpreInt_edges = results.RpreInt_edges{iRow};
    RpreInt_freq = results.RpreInt_freq{iRow};
    PUE = results.PUE{iRow};
    pNoInf = results.pNoInf{iRow};
    pNoInfOrCases = results.pNoInfOrCases{iRow};



    figure;
    x = mean([RpreInt_edges(1:end-1); RpreInt_edges(2:end)] );
    dx = x(2)-x(1);
    y = RpreInt_freq/sum(RpreInt_freq)/dx;
    bar(x, y)
    xlabel('pre-intervention reproduction number')
    ylabel('probability density')
    


    figure;
    tiledlayout(2, 2);
    
    nexttile;
    plot(t, Rt, 'Color', grey)
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('reproduction number')
    grid on
    
    nexttile;
    plot(t, It, 'Color', grey)
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('daily infections')
    grid on
    
    nexttile;
    plot(t, par.pReport * Zt, 'Color', grey)
    hold on
    plot(processed.t, processed.nCasesLoc )
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('expected daily cases')
    grid on
    
    nexttile;
    plot(t, Ct, 'Color', grey)
    hold on
    plot(processed.t, processed.nCasesLoc )
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('simulated daily cases')
    grid on


    iMinPlot = 25;          % don't plot the first part of the p(end of out break) curves which are >0 at the start of the outbreak
    scenarioKey = [6, 7, 8];

    figure;
    tiledlayout(3, 1);
    for iPlot = 1:length(scenarioKey)
        iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.iScenario == scenarioKey(iPlot) );

        t = results.t{iRow};
        par = results.par{iRow};
        PUE = results.PUE{iRow};
        pNoInf = results.pNoInf{iRow};
        pNoInfOrCases = results.pNoInfOrCases{iRow};

        nexttile;
        yyaxis left
        bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
        bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
        bh(2).FaceColor = [0 0 1]; 
        ylabel('reported daily cases')
        yyaxis right
        plot(t(iMinPlot:end), PUE(iMinPlot:end), 'b-', t(iMinPlot:end), pNoInf(iMinPlot:end), 'b--', t(iMinPlot:end), pNoInfOrCases(iMinPlot:end), 'b:' )
        xline(par.tRampStart, 'k:');
        ylabel('P(end of outbreak)')
        if iPlot == length(scenarioKey)
            legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or reported cases", 'Location', 'southeast')
        end
        xlim( [processed.t(1)-1, t(end)  ] )
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        grid on
    end


    scenarioKey = 1:6;
    colOrd = colororder;
    col = repelem( colOrd(1:3, :), 2, 1);
    lt = repmat( ["-"; "--"], 3, 1);

    figure;
    yyaxis left
    bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
    bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
    bh(2).FaceColor = [0 0 1]; 
    ylabel('reported daily cases')
    yyaxis right

    for iPlot = 1:length(scenarioKey)
        iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.iScenario == scenarioKey(iPlot) );
        t = results.t{iRow};
        par = results.par{iRow};
        pNoInf = results.pNoInf{iRow};
        if isequal(par.RTD, 1)
            scLabel(iPlot) = string(sprintf('alpha=%.1f, no delay', par.pReport));
        else
            scLabel(iPlot) = string(sprintf('alpha=%.1f, delay', par.pReport));
        end
        plot(t(iMinPlot:end), pNoInf(iMinPlot:end), 'Color', col(iPlot, :), 'LineStyle', lt(iPlot), 'Marker', 'none' )
        hold on
    end

    xline(par.tRampStart, 'k:');
    ylabel('P(end of outbreak)')
    legend(["data - imported cases", "data - local cases", scLabel], 'Location', 'northwest')
    xlim( [processed.t(1)-1, t(end)  ] )
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    grid on

end
