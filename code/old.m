    % Probability of no further infections after time t given the benefit of
    % hindsight
    pEndHind = exp(-RpreInt.*Gamma);
    pEndHindAvg = mean(pEndHind);

    % Probability of ultimate extinction
    PUEHind = exp(-(1-PUE1).*RpreInt.*Gamma);
    PUEHindAvg = mean(PUEHind);
    
    