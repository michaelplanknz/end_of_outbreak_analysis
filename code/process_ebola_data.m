clear 
close all

rawDataFolder = "../raw_data/";
processedDataFolder = "../processed_data/";

fNameProcessed = "ebola_DRC_2018_processed.csv";

addpath(rawDataFolder);

% Read raw data
[initial_date, incidence_data] = getEbolaData();

% Set up date vector
iFirstCase = find(incidence_data > 0, 1, 'first');
dateStart = datetime(initial_date) + iFirstCase-1;
dateEnd = datetime(initial_date) + length(incidence_data)-1;

processed.t = (dateStart:dateEnd)';      

% Assume first reported case is imported and the rest are local
processed.nCasesLoc = [0; incidence_data(iFirstCase+1:end)];
processed.nCasesImp = [incidence_data(iFirstCase); zeros(length(processed.t)-1, 1) ];

processed = struct2table(processed);

% Save processed data
writetable(processed, processedDataFolder+fNameProcessed);



