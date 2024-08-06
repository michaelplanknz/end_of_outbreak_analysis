function [t, par] = getPar(outbreakLbl)


par.nParticles = 1e5;


if outbreakLbl == "covid_NZ_2020"
    
    date0 = datetime(2020, 2, 26);              % date of 1st case (simulation may start earlier than this because date of infection may be earlier)
    date1 = datetime(2020, 6, 30);           % End date for simulation (anything up to 10 August 2020 which was last day of zero reported cases)
    
    par.tMIQ = datetime(2020, 4, 10);               % Imported cases after this date will be ignored
    par.tInfImp = 7;                          % 3 - number of days imported cases are assumed to have been infectious in community before case notifiction (no assumptions about infectious period ending on notification date or starting after arrival date)
    
    t = date0-par.tInfImp:date1;
    
    par.resampleLag = 30;       % fixed lag resampling
    
    par.sigmaR = 0.05;          % S.D. in random walk step for R
    par.k = inf;                % overdispersion parameter for offspring distribution (set to inf for a Poission distribution)

    par.pReport = 0.5;          % Reporting probability
    par.obsModel = "negbin";    % negative binomial daily observation modle (with dispersion parameter kObs, set kObs = infty for Poisson)
    par.kObs = inf;             % overdispersion parameter for daily observed cases (set to inf for a Poisson distribution)
    
    [par.R_shape, par.R_scale] = gamShapeScale(1.8, 0.7 );                    % shape-scale parameters for prior for initial R
     
    % Generation time distribution parameters
    GTmax = 15;                     % maximum time for GT distribution 
    [GTD_shape, GTD_scale] = gamShapeScale(5.05, 1.94);
    pdfFnGTD = @(x)(gampdf(x, GTD_shape, GTD_scale ));

    % Infection to report time distribution parameters
    reportDelayFlag = 1;
    RTmax = 30;
    [RTD_shape, RTD_scale] = gamShapeScale(11.2, 4.72);                   % mean and variance of incubation plus onset to isolation plus isolation to reporting
    pdfFnRTD = @(x)(gampdf(x, RTD_shape, RTD_scale));

    
    par.preIntWindow = 14;              % days before ramp start to use as a window for estimating per-intervention Rt
    
    par.tRampStart = datetime(2020, 3, 22); % start of intervention-related ramp down
    par.tRampEnd = datetime(2020, 3, 27);   % end of intervention-related ramp down
    rampDrop = 1.45;                    % total expected drop in R during ramp down
    rampDropSD = 0.2;



elseif outbreakLbl == "ebola_DRC_2018"

    date0 = datetime(2018, 4, 5);              % date of 1st case (simulation may start earlier than this because date of infection may be earlier)
    date1 = datetime(2018, 8, 14);           % End date for simulation (anything up to 10 August 2020 which was last day of zero reported cases)
    
    par.tMIQ = NaT;               % Imported cases after this date will be ignored
    par.tInfImp = 0;                          % 3 - number of days imported cases are assumed to have been infectious in community before case notifiction (no assumptions about infectious period ending on notification date or starting after arrival date)
    
    t = date0-par.tInfImp:date1;
    
    par.resampleLag = 30;     % fixed lag resampling
    
    par.sigmaR = 0.02;          % S.D. in random walk step for R
    par.k = inf;                % overdispersion parameter for offspring distribution (set to inf for a Poission distribution)

    par.pReport = 1;          % Reporting probability
    par.obsModel = "bin";     %  binomial daily observation model 
    par.kObs = inf;             % overdispersion parameter for daily observed cases (set to inf for a Poisson distribution, ignored if par.obsModel is "bin")

    
    [par.R_shape, par.R_scale] = gamShapeScale(2.5, 1 );                    % shape-scale parameters for prior for initial R
    
    % Generation time distribution parameters
    GTmax = 70;                     % maximum time for GT distribution 
    [GTD_shape, GTD_scale] = gamShapeScale(15.3, 9.3);
    pdfFnGTD = @(x)(gampdf(x, GTD_shape, GTD_scale ));

    % Infection to report time distribution parameters
    reportDelayFlag = 0;
    RTmax = 30;
    [RTD_shape, RTD_scale] = gamShapeScale(11.2, 4.72);                   % mean and variance of incubation plus onset to isolation plus isolation to reporting
    pdfFnRTD = @(x)(gampdf(x, RTD_shape, RTD_scale));    
    
    par.preIntWindow = 14;              % days before ramp start to use as a window for estimating per-intervention Rt
    
    par.tRampStart = datetime(2018, 5, 8); % start of intervention-related ramp down
    par.tRampEnd = datetime(2018, 5, 9);   % end of intervention-related ramp down
    rampDrop = 2.0;                    % total expected drop in R during ramp down
    rampDropSD = 0.4;
else
    error(sprintf('Invalid outbreak label: %s', outbreakLbl));
end


% Consturct GT and reporting time distributions
par.GTD = discDist( pdfFnGTD, 1, GTmax ); % Probability mass of discretised GT distribution on integers 1, 2, ...
if reportDelayFlag == 1
    par.RTD = discDist( pdfFnRTD, 0, RTmax);
else
    par.RTD = 1;
end


% Construct vector of prior daily mean and s.d. change in Rt
nDays = days(par.tRampEnd-par.tRampStart);
dropPerDayMean = rampDrop / nDays;
dropPerDaySD = rampDropSD / sqrt(nDays);
rampCV = 0.2;
par.deltat = zeros(size(t));
par.deltat(t >= par.tRampStart & t < par.tRampEnd) = -dropPerDayMean;
par.sigmat = par.sigmaR * ones(size(t));
par.sigmat(t >= par.tRampStart & t < par.tRampEnd) = sqrt(par.sigmaR^2 +  dropPerDaySD^2 );


