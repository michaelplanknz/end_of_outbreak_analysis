

    h = figure;
histogram(RpreInt);
    xlabel('pre-intervention reproduction number')
    ylabel('probability density')

% x = 0:0.01:4;
% y = gampdf(x, RpreInt_sh, RpreInt_sc);
% plot(x, y)
% xlabel('pre-intervention reproduction number')
% ylabel('probability density')

grey = [0.7 0.7 0.7];
      h = figure;
    h.Position = [ 273   236   778   612];
    tiledlayout(2, 2);

    nexttile;
    plot(t, Rt(1:nToPlot, :), 'Color', grey)
    xline(par.tRampStart, 'k:');
%    xline(par.tRampEnd, 'k:');
    ylabel('reproduction number')
    grid on
    
    nexttile;
    plot(t, It(1:nToPlot, :), 'Color', grey)
    xline(par.tRampStart, 'k:');
%    xline(par.tRampEnd, 'k:');
    ylabel('daily infections')
    grid on
    
    nexttile;
    plot(t, par.pReport * Zt(1:nToPlot, :), 'Color', grey)
    hold on
    plot(processed.t, processed.nCasesLoc )
    xline(par.tRampStart, 'k:');
  %  xline(par.tRampEnd, 'k:');
    ylabel('expected daily cases')
    grid on
    
    nexttile;
    plot(t, Ct(1:nToPlot, :), 'Color', grey)
    hold on
    plot(processed.t, processed.nCasesLoc )
    xline(par.tRampStart, 'k:');
  %  xline(par.tRampEnd, 'k:');
    ylabel('simulated daily cases')
    grid on


iMinPlot = 30;
    h = figure;
          yyaxis left
        bh = bar(processed.t, [processed.nCasesImp, processed.nCasesLoc]', 'stacked' );
        bh(1).FaceColor =  [0.67578 0.84375 0.89844];  
        bh(2).FaceColor = [0 0 1]; 
        ylabel('reported daily cases')
        yyaxis right
        plot(t(iMinPlot:end), PUE(iMinPlot:end), 'b-', t(iMinPlot:end), pNoInf(iMinPlot:end), 'b--', t(iMinPlot:end), pNoInfOrCases(iMinPlot:end), 'b:' )
        xline(par.tRampStart, 'k:');
        ylabel('P(end of outbreak)')

            legend('data - imported cases', 'data - local cases',  'ultimate extinction', 'no future transmission', "no future transmission or reported cases", 'Location', 'southeast')

        xlim( [processed.t(1)-1, t(end)  ] )
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        grid on