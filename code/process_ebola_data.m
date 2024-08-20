clear 
close all

rawDataFolder = "../raw_data/";
processedDataFolder = "../processed_data/";

fNameRaw = 'linelist_20180620.xlsx';
fNameProcessed = "ebola_DRC_2018_processed.csv";

addpath(rawDataFolder);

% Read raw data
opts = detectImportOptions(rawDataFolder+fNameRaw,'ReadVariableNames',true);
raw = readtable(rawDataFolder+fNameRaw, opts);

onsetDate = raw.DateDeD_butDesSignesEtSympt_mes_;

% Inclusion criteria: confirmed or probable case and has valid symptom onset date
inclFlag = ~isnat(onsetDate) & (contains(raw.Classification_pid_miologique_, "1-Conf") | contains(raw.Classification_pid_miologique_, "2-Probable") );

% Set up date vector
dateStart = min(onsetDate(inclFlag));
dateEnd = datetime(2018, 7, 10);
processed.t = (dateStart:dateEnd)';      

edges = [processed.t; processed.t(end)+1];
nCases = histcounts(onsetDate(inclFlag), edges)';

% Assume first reported case is imported and the rest are local
processed.nCasesLoc = [0; nCases(2:end)];
processed.nCasesImp = [nCases(1); zeros(length(processed.t)-1, 1) ];

processed = struct2table(processed);

% Save processed data
writetable(processed, processedDataFolder+fNameProcessed);



