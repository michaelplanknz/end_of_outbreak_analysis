clear 
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFlag = false;

dataFolder = "../processed_data/";
resultsFolder = "../results/";
figuresFolder = "../figures/";

load(resultsFolder+"results.mat");


baseScenario = 1;
greyCol = [0.7 0.7 0.7];
lightGreen = [0.8 1 0.7];
darkGreen = [0.2 0.6 0.1];
lightRed = [1 0.8 0.7];
darkRed = [0.6 0.2 0.1];
lightBlue = [0.7 0.8 1];
darkBlue = [0.1 0.2 0.6];


if saveFlag
    fOut = resultsFolder+"table.tex";
else
    fOut = 'temp.tex';
end
fid = fopen(fOut, 'w');
fprintf(fid, '\\begin{tabular}{lll} \n');
fprintf(fid, '\\hline\n');








% Main figures

iMinPlot = 30;          % don't plot the first part of the p(end of out break) curves which are >0 at the start of the outbreak
iFig = 1;    
for iOutbreak = 1:2
    fNameData = sprintf('%s_processed.csv', outbreakLbl(iOutbreak));
    processed = readtable(dataFolder+fNameData);

    iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.sensitivityFlag == sensitivityFlag(iOutbreak) & results.iScenario == baseScenario);

    t = results.t{iRow};
    par = results.par{iRow};
    Rt_quantiles = results.Rt_quantiles{iRow};
    It_quantiles = results.It_quantiles{iRow};
    Zt_quantiles = results.Zt_quantiles{iRow};
    Ct_quantiles = results.Ct_quantiles{iRow};
    RpreInt_mean = results.RpreInt_mean{iRow};
    RpreInt_sd = results.RpreInt_sd{iRow};
    RpreInt_bins = results.RpreInt_bins{iRow};
    RpreInt_freq = results.RpreInt_freq{iRow};
    PUE = results.PUE_mean{iRow};
    pNoInf = results.pNoInf_mean{iRow};
    pNoInfOrCases = results.pNoInfOrCases_mean{iRow};
    PUE_quantiles = results.PUE_quantiles{iRow};
    pNoInf_quantiles = results.pNoInf_quantiles{iRow};
    pNoInfOrCases_quantiles = results.pNoInfOrCases_quantiles{iRow};


    h = figure(iFig);
    h.Position = [   190         412        1128         315];
    tiledlayout(1, 3, "TileSpacing", "Compact");

    nexttile;
    [x, y] = getFillArgs(t, Rt_quantiles(1, :), Rt_quantiles(3, :) );
    fill( x, y, lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, Rt_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    xline(par.tRampStart, 'k:');
    ylabel('reproduction number')
    xlim([processed.t(1)-1, t(end) ]);
    grid on
    title('(a)')
   
    nexttile;
    [x, y] = getFillArgs(t, It_quantiles(1, :), It_quantiles(3, :) );
    fill( x, y, lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, It_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    xline(par.tRampStart, 'k:');
    ylabel('daily new local infections')
    xlim([processed.t(1)-1, t(end) ]);
    grid on
    title('(b)')

    nexttile;
    yyaxis left
    [x, y] = getFillArgs(t, Ct_quantiles(1, :), Ct_quantiles(3, :) );
    fill( x, y, lightGreen, 'LineStyle', 'none'  )
    hold on
    plot(t, Ct_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
    plot(processed.t, processed.nCasesLoc, 'k.' )
    if outbreakLbl(iOutbreak) == "covid_NZ_2020"
        ylim([0 100])
    else
        ylim([0 10])
    end
    ylabel('daily local case notifications')
    yyaxis right
    [x, y] = getFillArgs(t(iMinPlot:end), pNoInf_quantiles(1, iMinPlot:end), pNoInf_quantiles(3, iMinPlot:end) );
    fill(x, y, lightBlue, 'LineStyle', 'none' )
    plot(t(iMinPlot:end), pNoInf(iMinPlot:end), 'b-')    
    ylabel('end of outbreak probability (P_0)')
    xline(par.tRampStart, 'k:');
    %xlim([processed.t(1)-1, processed.t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
    xlim([processed.t(1)-1, t(end) ]);
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    grid on
    title('(c)')

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
    if outbreakLbl(iOutbreak) == "covid_NZ_2020"
       ylim([0 100])
    else
        ylim([0 10])
    end
    yyaxis right

    scLabel = strings(1, length(scenarioKey));
    for iPlot = 1:length(scenarioKey)
        iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.sensitivityFlag == sensitivityFlag(iOutbreak) & results.iScenario == scenarioKey(iPlot) );
        t = results.t{iRow};
        par = results.par{iRow};
        pNoInf = results.pNoInf_mean{iRow};
        RTmean = sum( (0:length(par.RTD)-1).*par.RTD );
        scLabel(iPlot) = "\alpha="+string(sprintf('%.1f, ', par.pReport))+"t_n="+string(sprintf('%.1f days', RTmean));
        plot(t(iMinPlot:end), pNoInf(iMinPlot:end), 'Color', col(iPlot, :), 'LineStyle', lt(iPlot), 'Marker', 'none' )
        hold on
    end

    xline(par.tRampStart, 'k:');
    ylabel('end of outbreak probability (P_0)')
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






    scenarioKey = [1, 7, 8];

    h = figure(iFig);
    h.Position = [           437          46        1049         918];
    tiledlayout(3, 1, "TileSpacing", "Compact");
    for iPlot = 1:length(scenarioKey)
        iRow = find(results.outbreak == outbreakLbl(iOutbreak) & results.sensitivityFlag == sensitivityFlag(iOutbreak) & results.iScenario == scenarioKey(iPlot) );

        t = results.t{iRow};
        par = results.par{iRow};
        PUE = results.PUE_mean{iRow};
        pNoInf = results.pNoInf_mean{iRow};
        pNoInfOrCases = results.pNoInfOrCases_mean{iRow};

        nexttile;
        yyaxis left
        bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
        bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
        bh(2).FaceColor = [0 0 1]; 
        ylabel('daily case notifications')
        if outbreakLbl(iOutbreak) == "covid_NZ_2020"
            ylim([0 100])
        else
            ylim([0 10])
        end
        yyaxis right
        plot(t(iMinPlot:end), PUE(iMinPlot:end), 'b-', t(iMinPlot:end), pNoInf(iMinPlot:end), 'b--', t(iMinPlot:end), pNoInfOrCases(iMinPlot:end), 'b:' )
        xline(par.tRampStart, 'k:');
        ylabel('end of outbreak probability')
        if iPlot == length(scenarioKey)
            legend('data - imported cases', 'data - local cases',  'PUE', 'P_0', "P_{00}", 'Location', 'southeast')
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




    if ~sensitivityFlag(iOutbreak)
    
        scenarioKey = 1:6;
        rowKey = find(results.outbreak == outbreakLbl(iOutbreak) & ismember(results.iScenario, scenarioKey) );
        if outbreakLbl(iOutbreak) == "covid_NZ_2020"
            fprintf(fid, '\\multicolumn{3}{l}{ \\bf Covid-19} \\\\ \n');
        elseif outbreakLbl(iOutbreak) == "ebola_DRC_2018"
            fprintf(fid, '\\multicolumn{3}{l}{ \\bf Ebola} \\\\ \n');
        end
        fprintf(fid, ' &  $t_n=%.1f$ days & $t_n=%.1f$ days \\\\ \n', results.par{rowKey(1)}.RTmean, results.par{rowKey(2)}.RTmean);
        fprintf(fid, '$\\alpha=%.1f$ & %s & %s \\\\ \n ', results.par{rowKey(1)}.pReport, string(results.tpNoInf95{rowKey(1)}) , string(results.tpNoInf95{rowKey(2)}) );
        fprintf(fid, '$\\alpha=%.1f$ & %s & %s \\\\ \n ', results.par{rowKey(3)}.pReport, string(results.tpNoInf95{rowKey(3)}) , string(results.tpNoInf95{rowKey(4)}) );
        fprintf(fid, '$\\alpha=%.1f$ & %s & %s \\\\ \n ', results.par{rowKey(5)}.pReport, string(results.tpNoInf95{rowKey(5)}) , string(results.tpNoInf95{rowKey(6)}) );
        fprintf(fid, '\\hline\n');
    end
end

fprintf(fid, '\\end{tabular} \n');
fclose('all');

type(fOut);

if ~saveFlag
    delete(fOut);
end











% Additional figures
ttls = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"];

for iOutbreak = 1:2
    fNameData = sprintf('%s_processed.csv', outbreakLbl(iOutbreak));
    processed = readtable(dataFolder+fNameData);
    
    scenarioSet{1} = 1:6;
    scenarioSet{2} = [1 7 8];
    
    for iSet = 1:2
        h = figure(iFig);
        h.Position = [     255          41        1128        179*length(scenarioSet{iSet}) ];
        tiledlayout(length(scenarioSet{iSet}), 3, "TileSpacing", "Compact");
    
        for ii = 1:length(scenarioSet{iSet})
            iScenario = scenarioSet{iSet}(ii);
            iRow = iScenario + nScenarios*(iOutbreak-1);
    
            t = results.t{iRow};
            par = results.par{iRow};
            Rt_quantiles = results.Rt_quantiles{iRow};
            It_quantiles = results.It_quantiles{iRow};
            Zt_quantiles = results.Zt_quantiles{iRow};
            It_imp_quantiles = results.It_imp_quantiles{iRow};
            Zt_imp_quantiles = results.Zt_imp_quantiles{iRow};
            Ct_quantiles = results.Ct_quantiles{iRow};
            pNoInf_quantiles = results.pNoInf_quantiles{iRow};
            PUE = results.PUE_mean{iRow};
            pNoInf = results.pNoInf_mean{iRow};
            pNoInfOrCases = results.pNoInfOrCases_mean{iRow};
    
    
            
    
        
            nexttile;
            [x, y] = getFillArgs(t, Rt_quantiles(1, :), Rt_quantiles(3, :) );
            fill( x, y, lightGreen, 'LineStyle', 'none'  )
            hold on
            plot(t, Rt_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
            xline(par.tRampStart, 'k:');
            ylabel('reproduction number')
            xlim([processed.t(1)-1, t(end) ]);
            grid on
           
            nexttile;
            [x, y] = getFillArgs(t, It_quantiles(1, :), It_quantiles(3, :) );
            fill( x, y, lightGreen, 'LineStyle', 'none'  )
            hold on
            plot(t, It_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
            xline(par.tRampStart, 'k:');
            ylabel('daily new local infections')
            xlim([processed.t(1)-1, t(end) ]);
            grid on
            title(ttls(ii));
        
            nexttile;
            yyaxis left
            [x, y] = getFillArgs(t, Ct_quantiles(1, :), Ct_quantiles(3, :) );
            fill( x, y, lightGreen, 'LineStyle', 'none'  )
            hold on
            plot(t, Ct_quantiles(2, :), 'Color', darkGreen, 'LineStyle', '-')
            plot(processed.t, processed.nCasesLoc, 'k.' )
            if outbreakLbl(iOutbreak) == "covid_NZ_2020"
                ylim([0 100])
            else
                ylim([0 10])
            end
            ylabel('daily local case notifications')
            yyaxis right
            [x, y] = getFillArgs(t(iMinPlot:end), pNoInf_quantiles(1, iMinPlot:end), pNoInf_quantiles(3, iMinPlot:end) );
            fill(x, y, lightBlue, 'LineStyle', 'none' )
            plot(t(iMinPlot:end), pNoInf(iMinPlot:end), 'b-')    
            ylabel('end of outbreak probability (P_0)')
            xline(par.tRampStart, 'k:');
            %xlim([processed.t(1)-1, processed.t(find(processed.nCasesLoc > 0, 1, 'last')+7) ]);
            xlim([processed.t(1)-1, t(end) ]);
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'b';
            grid on
    
        end
        if saveFlag 
          saveas(h, figuresFolder+sprintf('fig%i.png', iFig));
        end
        iFig = iFig+1;
    end
end


