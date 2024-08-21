function PUE = calcPUE(R, k)

% Calculate probability of ultimate extinction for an outbreak starting
% with a single fully infectious seed case with reproduction numbers in the array R and dispersion parameter k (set k=inf for Poisson)

% Maximum number of PUE solutions to calculate
nMax = 1000;

Runique = unique(R);
nUnique = length(Runique);


% If the number of unique values is less than the maximum allowed number of
% solutions, just calculate PUE for each unique value of R 
if nUnique <= nMax
    Rarr = Runique;
else    % otherwise use nMax equally spaced values for R and then later interpolate 
    Rarr = linspace(min(R), max(R), nMax);
end
nPoints = length(Rarr);
PUEarr = zeros(1, nPoints);

x0 = 0.1;
opts = optimset('display', 'off');
if isfinite(k)      % NegBin PUE
    for iPoint = 1:nPoints
        % pNB = k/(k+Rarr(iPoint));
        % myF = @(x)( (pNB/(1-(1-pNB)*x))^k - x );
        myF = @(x)( (k/(k+Rarr(iPoint)*(1-x)))^k - x);
        PUEarr(iPoint) = fsolve(myF, x0, opts);
    end
else            % Poission PUE
    for iPoint = 1:nPoints
      myF = @(x)(exp(Rarr(iPoint).*(x-1)) - x);
      PUEarr(iPoint) = fsolve(myF, x0, opts);
  end
end

if nPoints > 1
    PUE = interp1(Rarr, PUEarr, R);
else
    PUE = PUEarr;
end
