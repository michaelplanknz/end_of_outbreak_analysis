# Estimation of end-of-outbreak probabilities in the presence of case underascertainment and reporting lags

Data and code to reproduce the results in the article "Estimation of end-of-outbreak probabilities in the presence of case underascertainment and reporting lags" by Plank et al. 

A preprint of the article is available at 
[https://doi.org/10.48550/arXiv.2409.16531](
https://doi.org/10.48550/arXiv.2409.16531). 

Results in version 1 [v1] of the preprint were generated using the version of the repository tagged v1.0. 

Results in version 2 [v2] of the preprint and the final published article were generated using the version of the repository tagged v2.0.

# Abstract

Towards the end of an infectious disease outbreak, when a period has elapsed without new case notifications, a key question for public health policy makers is whether the outbreak can be declared over. This requires the benefits of a declaration (e.g., relaxation of outbreak control measures) to be balanced against the risk of a resurgence in cases. To support this decision making, mathematical methods have been developed to quantify the end-of-outbreak probability. Here, we propose a new approach to this problem that accounts for a range of features of real-world outbreaks, specifically: (i) incomplete case ascertainment; (ii) reporting delays; (iii) individual heterogeneity in transmissibility; and (iv) whether cases were imported or infected locally. We showcase our approach using two case studies: Covid-19 in New Zealand in 2020, and Ebola virus disease in the Democratic Republic of the Congo in 2018. In these examples, we found that the date when the estimated probability of no future infections reached 95% was relatively consistent across a range of modelling assumptions. This suggests that our modelling framework can generate robust quantitative estimates that can be used by policy advisors, alongside other sources of evidence, to inform end-of-outbreak declarations.

# How to use this repository

Processed data (case time series) for the Covid-19 and Ebola virus disease (EVD) case studies are stored in the folder `processed_data`. 

Matlab code is in the folder `code`. 

The script `main.m` reads in the data and runs the models for both case studies and for the set of scenarios shown in the article. `main.m` saves the model output in the file `results.mat` in the `results` folder. 

The script `makePlots.m` reads in `results.mat` and plots the graphs that appear in the article. Set the variable `saveFlag` to `true` to save the graphs as .png files in the `figures` folder, and to write a LaTeX table with the end-of-outbreak declaration dates in `table.tex`. If `saveFlag` is `false`, the grpahs will be plotted but no results will be saved to file. 




