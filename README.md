# Estimating end-of-outbreak probability

Data and code to reproduce the results in the article ``Robust estimation of end-of-outbreak probabilities in the presence of case underascertainment and reporting lags''. 

A preprint of the article is available at XXX. Results in this version were generated using the version of the repository tagged v1.0.

# How to use this repository

Processed data (case time series) for the Covid-19 and Ebola virus disease (EVD) case studies are stored in the folder `processed_data`. 

Matlab code is in the folder `code`. 

The script `main.m` reads in the data and runs the models for both case studies and for the set of scenarios shown in the article. `main.m` saves the model output in the file `results.mat` in the `results` folder. 

The script `makePlots.m` reads in `results.mat` and plots the graphs that appear in the article. Set the variable `saveFlag` to `true` to save the graphs as .png files in the `figures` folder, and to write a LaTeX table with the end-of-outbreak declaration dates in `table.tex`. If `saveFlag` is `false`, the grpahs will be plotted but no results will be saved to file. 




