clear
close all

rawDataFolder = "../raw_data/";
processedDataFolder = "../processed_data/";

fNameRaw = 'covid_cases_NZ_march-may2020.csv';
fNameProcessed = 'covid_NZ_2020_processed.csv';

processed.t = (datetime(2020, 2, 26):datetime(2020, 6, 30))';      

% Read raw data
opts = detectImportOptions(rawDataFolder+fNameRaw);
opts = setvartype(opts, {'CaseStatus', 'Sex', 'AgeGroup', 'District', 'OverseasTravel', 'InfectionStatus'}, 'categorical');
raw = readtable(rawDataFolder+fNameRaw, opts);
raw.OverseasTravel(isundefined(raw.OverseasTravel)) = "No";

% Extract case time series
casesLoc = groupsummary(raw(raw.OverseasTravel == "No", :), 'ReportDate', "sum", "NumberOfCasesReported");
casesImp = groupsummary(raw(raw.OverseasTravel == "Yes", :), 'ReportDate', "sum", "NumberOfCasesReported");


processed.nCasesLoc = zeros(size(processed.t));
processed.nCasesLoc(ismember(processed.t, casesLoc.ReportDate)) = casesLoc.sum_NumberOfCasesReported';
processed.nCasesImp = zeros(size(processed.t));
processed.nCasesImp(ismember(processed.t, casesImp.ReportDate)) = casesImp.sum_NumberOfCasesReported';

processed = struct2table(processed);

% Save processed data
writetable(processed, processedDataFolder+fNameProcessed);



