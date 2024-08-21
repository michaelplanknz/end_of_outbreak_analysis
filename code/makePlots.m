clear 
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFlag = true;

dataFolder = "../processed_data/";
resultsFolder = "../results/";
figuresFolder = "../figures/";

load(resultsFolder+"results.mat");


baseScenario = 1;
greyCol = [0.7 0.7 0.7];
lightGreen = [0.8 1 0.7];
darkGreen = [0.2 0.6 0.1];

iFig = 1;    
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
    Rt_quantiles = results.Rt_quantiles{iRow};
    It_quantiles = results.It_quantiles{iRow};
    Zt_quantiles = results.Zt_quantiles{iRow};
    Ct_quantiles = results.Ct_quantiles{iRow};
%     RpreInt_sh = results.RpreInt_sh{iRow};
%     RpreInt_sc = results.RpreInt_sc{iRow};
    PUE = results.PUE{iRow};
    pNoInf = results.pNoInf{iRow};
    pNoInfOrCases = results.pNoInfOrCases{iRow};

%     h = figure(iFig);
%     x = 0:0.01:4;
%     y = gamcdf(x, RpreInt_sh, RpreInt_sc);
%     xlabel('pre-intervention R')
%     ylabel('probability density')
%     grid on
%     if saveFlag 
%         saveas(h, figuresFolder+sprintf('fig%i.png', iFig));
%     end
%     iFig = iFig+1;


%     h = figure(iFig);
%     h.Position = [ 273   236   778   612];
%     tiledlayout(2, 2);
% 
%     nexttile;
%     plot(t, Rt, 'Color', greyCol)
%     xline(par.tRampStart, 'k:');
%     ylabel('reproduction number')
%     xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
%     grid on
%     title('(a)')
%     
%     nexttile;
%     plot(t, It, 'Color', greyCol)
%     xline(par.tRampStart, 'k:');
%     ylabel('daily local infections')
%     xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
%     grid on
%     title('(b)')
%     
%     nexttile;
%     plot(t, par.pReport * Zt, 'Color', greyCol)
%     hold on
%     plot(processed.t, processed.nCasesLoc, 'bo' )
%     xline(par.tRampStart, 'k:');
%     ylabel('expected daily local cases')
%     xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
%     grid on
%     title('(c)')
%     
%     nexttile;
%     plot(t, Ct, 'Color', greyCol)
%     hold on
%     plot(processed.t, processed.nCasesLoc, 'bo' )
%     xline(par.tRampStart, 'k:');
%     ylabel('simulated daily local cases')
%     xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
%     grid on
%     title('(d)')
%     if saveFlag 
%         saveas(h, figuresFolder+sprintf('fig%i.png', iFig));
%     end
%     iFig = iFig+1;

    % Alternative verison with median and 90% CI
    h = figure(iFig);
    h.Position = [ 273   236   778   612];
    tiledlayout(2, 2);

    nexttile;
    fill( [t, fliplr(t)], [Rt_quantiles(1, :), fliplr(Rt_quantiles(3, :))], lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, Rt_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    xline(par.tRampStart, 'k:');
    ylabel('reproduction number')
    xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
    grid on
    title('(a)')
   
    nexttile;
    fill( [t, fliplr(t)], [It_quantiles(1, :), fliplr(It_quantiles(3, :))], lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, It_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    xline(par.tRampStart, 'k:');
    ylabel('daily local infections')
    xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
    grid on
    title('(b)')
    
    nexttile;
    fill( [t, fliplr(t)], par.pReport * [Zt_quantiles(1, :), fliplr(Zt_quantiles(3, :))], lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, par.pReport * Zt_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    plot(processed.t, processed.nCasesLoc, 'bo')
    xline(par.tRampStart, 'k:');
    ylabel('expected daily local cases')
    xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
    grid on
    title('(c)')
    
    nexttile;
    fill( [t, fliplr(t)], [Ct_quantiles(1, :), fliplr(Ct_quantiles(3, :))], lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, Ct_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    plot(processed.t, processed.nCasesLoc, 'bo' )
    xline(par.tRampStart, 'k:');
    ylabel('simulated daily local cases')
    xlim([t(1), t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
    grid on
    title('(d)')
    if saveFlag 
        saveas(h, figuresFolder+sprintf('fig%i.png', iFig));
    end
    iFig = iFig+1;






    iMinPlot = 30;          % don't plot the first part of the p(end of out break) curves which are >0 at the start of the outbreak
    scenarioKey = [1, 7, 8];

    h = figure(iFig);
    h.Position = [           437          46        1049         918];
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
        ylabel('daily case notifications')
        if iOutbreak == 1
            ylim([0 100])
        else
            ylim([0 10])
        end
        yyaxis right
        plot(t(iMinPlot:end), PUE(iMinPlot:end), 'b-', t(iMinPlot:end), pNoInf(iMinPlot:end), 'b--', t(iMinPlot:end), pNoInfOrCases(iMinPlot:end), 'b:' )
        xline(par.tRampStart, 'k:');
        ylabel('P(end of outbreak)')
        if iPlot == length(scenarioKey)
            legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or notifications", 'Location', 'southeast')
        end
        xlim( [processed.t(1)-1, t(end)  ] )
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        grid on
        if iPlot == 1
           title('(a)')
        elseif iPlot == 2
           title('(b)')
        else
            title('(c)')
        end

    end
    if saveFlag 
        saveas(h, figuresFolder+sprintf('fig%i.png', iFig));
    end
    iFig = iFig+1;

    scenarioKey = 1:6;
    colOrd = colororder;
    col = repelem( colOrd(1:3, :), 2, 1);
    lt = repmat( ["-"; "--"], 3, 1);

    h = figure(iFig);
    h.Position = [          242         308        1277         587];
    yyaxis left
    bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
    bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
    bh(2).FaceColor = [0 0 1]; 
    ylabel('daily case notifications')
    if iOutbreak == 1
       ylim([0 100])
    else
        ylim([0 10])
    end
    yyaxis right

    for iPlot = 1:length(scenarioKey)
        iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.iScenario == scenarioKey(iPlot) );
        t = results.t{iRow};
        par = results.par{iRow};
        pNoInf = results.pNoInf{iRow};
        RTmean = sum( (0:length(par.RTD)-1).*par.RTD );
        scLabel(iPlot) = "\alpha="+string(sprintf('%.1f, ', par.pReport))+"t_n="+string(sprintf('%.1f days', RTmean));
        plot(t(iMinPlot:end), pNoInf(iMinPlot:end), 'Color', col(iPlot, :), 'LineStyle', lt(iPlot), 'Marker', 'none' )
        hold on
    end

    xline(par.tRampStart, 'k:');
    ylabel('P(end of outbreak)')
    legend(["data - imported cases", "data - local cases", scLabel], 'Location', 'southeast')
    xlim( [processed.t(1)-1, t(end)  ] )
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    grid on
    if saveFlag 
        saveas(h, figuresFolder+sprintf('fig%i.png', iFig));
    end
    iFig = iFig+1;
end
