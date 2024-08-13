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

    t = results{iRow}.t;
    par = results{iRow}.par;
    Rt = results{iRow}.Rt;
    It = results{iRow}.It;
    Zt = results{iRow}.Zt;
    Ct = results{iRow}.Ct;
    RpreInt_edges = results{iRow}.RpreInt_edges;
    RpreInt_freq = results{iRow}.RpreInt_freq;
    PUE = results{iRow}.PUE;
    pNoInf = results{iRow}.pNoInf;
    pNoInfOrCases = results{iRow}.pNoInfOrCases;



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
    
    nexttile;
    plot(t, It, 'Color', grey)
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('daily infections')
    
    nexttile;
    plot(t, par.pReport * Zt, 'Color', grey)
    hold on
    plot(t, results{iRow}.nCasesLoc )
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('expected daily cases')
    
    nexttile;
    plot(t, Ct, 'Color', grey)
    hold on
    plot(processed.t, processed.nCasesLoc )
    xline(par.tRampStart, 'k:');
    xline(par.tRampEnd, 'k:');
    ylabel('simulated daily cases')



    iMinPlot = 25;          % don't plot the first part of the p(end of out break) curves which are >0 at the start of the outbreak
    scenarioKey = [6, 7, 8];

    figure;
    tiledlayout(3, 1);
    for iTile = 1:length(scenarioKey)
            iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.iScenario == scenarioKey(iTile) );

            t = results{iRow}.t;
            par = results{iRow}.par;
            PUE = results{iRow}.PUE;
            pNoInf = results{iRow}.pNoInf;
            pNoInfOrCases = results{iRow}.pNoInfOrCases;
    
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
            legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or reported cases", 'Location', 'northwest')
            xlim( [processed.t(1)-1, t(end)  ] )
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'k';
    end


    scenarioKey = 1:6;
    colOrd = colororder;
    col = repelem( colOrd(1:3, :), [2, 1]);
    lt = repmat( ["-"; "--"], 3, 1);

    figure;
    yyaxis left
    bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
    bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
    bh(2).FaceColor = [0 0 1]; 
    ylabel('reported daily cases')
    yyaxis right

    for iCurve = 1:length(scenarioKey)
        t = results{iRow}.t;
        pNoInf = results{iRow}.pNoInf;
        plot(t(iMinPlot:end), pNoInf(iMinPlot:end), 'LineColor', col(iCurve, :), 'LineType', lt(iCurve) )
        hold on
    end

    xline(par.tRampStart, 'k:');
    ylabel('P(end of outbreak)')
    legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or reported cases", 'Location', 'northwest')
    xlim( [processed.t(1)-1, t(end)  ] )
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

end
