function scenarioTab = getScenarios(outbreakLbl)



if outbreakLbl == "covid_NZ_2020"
    scenarioTab.pReport = [1; 1; 0.7; 0.7; 0.4; 0.4; 0.4; 0.4];
    scenarioTab.RT_mean = [0; 11.2; 0; 11.2; 0; 11.2; 11.2; 11.2];
elseif outbreakLbl == "ebola_DRC_2018"
    scenarioTab.pReport = [1; 1; 0.9; 0.9; 0.8; 0.8; 0.8; 0.8];
    scenarioTab.RT_mean = [0; 6.2; 0; 6.2; 0; 6.2; 6.2; 6.2];
else
    error(sprintf('Invalid outbreak label: %s', outbreakLbl));
end

scenarioTab.k = [inf; inf; inf; inf; inf; inf; 1; 0.2];
scenarioTab = struct2table(scenarioTab);

