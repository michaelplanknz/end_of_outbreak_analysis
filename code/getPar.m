function [t, par] = getPar(outbreakLbl, iScenario)


par.nParticles = 1e5;


k_scenarios = [inf; inf; inf; inf; inf; inf; 1; 0.2];

if outbreakLbl == "covid_NZ_2020"
    % Set arrays of parameters varied from one scenairo to the next
    pReport_scenarios = [0.4; 0.4; 0.7; 0.7; 1; 1; 0.4; 0.4];
    RTmean_scenarios = [7.7; 11.2; 7.7; 11.2; 7.7; 11.2; 7.7; 7.7];
    RTsd_scenarios   = [3.2; 4.7; 3.2; 4.7; 3.2; 4.7; 3.2; 3.2];

    
    date0 = datetime(2020, 2, 26);              % date of 1st case (simulation may start earlier than this because date of infection may be earlier)
    date1 = datetime(2020, 6, 30);           % End date for simulation (anything up to 10 August 2020 which was last day of zero reported cases)
    
    par.resampleLag = 30;       % fixed lag resampling
 
    par.tRampStart = datetime(2020, 3, 23);
    par.rampDur = 7;                    % duration of intervention-related change in Rt
    
    par.sigmaR = 0.05;                % S.D. in random walk step for log Rt
    par.sigmaR_control = 0.2;          % S.D. in random walk step for log Rt during rapid intervention-associated change
    par.deltaR_control = -0.1;         % Mean of random walk step for log Rt during rapid intervention-associated change

    par.k = k_scenarios(iScenario);                % overdispersion parameter for offspring distribution (set to inf for a Poission distribution)

    par.pReport = pReport_scenarios(iScenario);          % Reporting probability
    par.obsModel = "negbin";    % negative binomial daily observation modle (with dispersion parameter kObs, set kObs = infty for Poisson)
    par.kObs = inf;             % overdispersion parameter for daily observed cases (set to inf for a Poisson distribution)
    
    [par.R_shape, par.R_scale] = gamShapeScale(2, 1 );                    % shape-scale parameters for prior for initial R
     
    % Generation time distribution parameters
    GTmax = 15;                     % maximum time for GT distribution 
    par.GTmean = 5.05;
    par.GTsd = 1.94;
    [GTD_shape, GTD_scale] = gamShapeScale(par.GTmean, par.GTsd);
    pdfFnGTD = @(x)(gampdf(x, GTD_shape, GTD_scale ));

    % Infection to report time distribution parameters
    RTmax = 25;
    par.RTmean = RTmean_scenarios(iScenario);
    par.RTsd = RTsd_scenarios(iScenario);
    [RTD_shape, RTD_scale] = gamShapeScale(par.RTmean, par.RTsd);       
    pdfFnRTD = @(x)(gampdf(x, RTD_shape, RTD_scale));

    par.filterImportsFlag = false;
    par.tInfImp = 5;                          % number of days imported cases are assumed to have been infectious in community before case notifiction (no assumptions about infectious period ending on notification date or starting after arrival date)

    par.relInfImp = 0.5;
    par.tMIQ = datetime(2020, 4, 10);               % Imported cases after this date will be ignored

    par.preIntWindow = 14;              % days before ramp start to use as a window for estimating per-intervention Rt
    




elseif outbreakLbl == "ebola_DRC_2018"
    pReport_scenarios = [0.8; 0.8; 0.9; 0.9; 1; 1; 0.8; 0.8];
    RTmean_scenarios = [6.2; 11.2; 6.2; 11.2; 6.2; 11.2; 6.2; 6.2 ];
    RTsd_scenarios =   [1.6; 4.3; 1.6; 4.3; 1.6; 4.3; 1.6; 1.6];

    date0 = datetime(2018, 4, 5);              % date of 1st case (simulation may start earlier than this because date of infection may be earlier)
    date1 = datetime(2018, 8, 31);           % End date for simulation 
    
    par.resampleLag = 50;     % fixed lag resampling

    par.tRampStart = datetime(2018, 5, 8); % start of intervention-related change in Rt
    par.rampDur = 7;                    % duration of intervention-related change in Rt

    par.sigmaR = 0.05;                % S.D. in random walk step for log Rt
    par.sigmaR_control = 0.2;          % S.D. in random walk step for log Rt during rapid intervention-associated change
    par.deltaR_control = -0.1;         % Mean of random walk step for log Rt during rapid intervention-associated change

    par.obsModel = "negbin";     %  binomial daily observation model 
    par.kObs = inf;             % overdispersion parameter for daily observed cases (set to inf for a Poisson distribution, ignored if par.obsModel is "bin")
    
    [par.R_shape, par.R_scale] = gamShapeScale(2.5, 1 );                    % shape-scale parameters for prior for initial R
    
    % Generation time distribution parameters
    GTmax = 50;                     % maximum time for GT distribution 
    par.GTmean = 15.3;
    par.GTsd = 9.3;
    [GTD_shape, GTD_scale] = gamShapeScale(par.GTmean, par.GTsd);
    pdfFnGTD = @(x)(gampdf(x, GTD_shape, GTD_scale ));

    % Infection to report time distribution parameters
    RTmax = 30;
    par.RTmean = RTmean_scenarios(iScenario);
    par.RTsd = RTsd_scenarios(iScenario);
    [RTD_shape, RTD_scale] = gamShapeScale(par.RTmean, par.RTsd);  
    pdfFnRTD = @(x)(gampdf(x, RTD_shape, RTD_scale));  

    par.filterImportsFlag = false;
    par.tInfImp = round(RTmean_scenarios(iScenario));                          % number of days imported cases are assumed to have been infectious in community before case notifiction (no assumptions about infectious period ending on notification date or starting after arrival date)
    par.relInfImp = 1;
    par.tMIQ = NaT;                                      % Imported cases after this date will be ignored
    
    par.preIntWindow = 14;              % days before ramp start to use as a window for estimating per-intervention Rt
else
    error(sprintf('Invalid outbreak label: %s', outbreakLbl));
end


par.k = k_scenarios(iScenario);                % overdispersion parameter for offspring distribution (set to inf for a Poission distribution)
par.pReport = pReport_scenarios(iScenario);          % Reporting probability


% Time vector
t = date0-11:date1;                         % Pad with a few days at the start to allow for setting the infection date for imported cases to be earlier than the report date


% Consturct GT and reporting time distributions
par.GTD = discDist( pdfFnGTD, 1, GTmax ); % Probability mass of discretised GT distribution on integers 1, 2, ...
if RTmean_scenarios(iScenario) > 0
    par.RTD = discDist( pdfFnRTD, 0, RTmax);
else
    par.RTD = 1;
end


% Construct vector of prior daily mean and s.d. change in Rt
par.deltat = zeros(size(t));
par.deltat(t >= par.tRampStart & t < par.tRampStart+par.rampDur) = par.deltaR_control;

par.sigmat = par.sigmaR * ones(size(t));
par.sigmat(t >= par.tRampStart & t < par.tRampStart + par.rampDur) = par.sigmaR_control;


