% Calculate gamma function using full data (hindsight)
GTDF = [1, 1-cumsum(par.GTD)];      % Pad the survival function with a leading 1 as the first element applies to t  he previous days infections, which have all of their transmission potential remaining the next day

Gamma = conv2(Yt, GTDF);
Gamma = [zeros(par.nParticles, 1), Gamma(:, 1:length(t)-1)];          % offset by one so a case on day t contributes to future infection potential starting on day t+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Probability of no further infections after time t given the benefit of
    % hindsight
    pEndHind = exp(-RpreInt.*Gamma);
    pEndHindAvg = mean(pEndHind);

    % Probability of ultimate extinction
    PUEHind = exp(-(1-PUE1).*RpreInt.*Gamma);
    PUEHindAvg = mean(PUEHind);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code snippets from different methods for calculation of RpreInt



