# Transferrin-calculations-with-R

Here is DOI for all 3 scripts:
https://doi.org/10.5281/zenodo.999668

#### File #1
Standard_curve_Resazurin.R 

R-script includes mathematical transformation, linear regression building and visualization of the standard curve in resazurin-based time-kill assay. Input is measurements of fluorescence in serially diluted broth culture obtained from FLUOstar Omega plate reader (BMG LABTECH GmbH, Germany).


##### File #2
Resazurin_time-kill_data_analysis.R

R-script to analyse data from Resazurin-based time-kill assay. Raw data of fluorescence measurements obtained from FLUOstar Omega plate reader (BMG LABTECH GmbH, Germany) in txt format.


##### File #3
Antibiotic_dynamics.R

R-script for data analysis in 20-days passage experiment. The passage was started from 1/4 MIC of antibiotics concentration and every 5 days the concentration was increased up to the maximum concentration which bacteria are able to tolerate.

##### File #4
KP-nPDR3_Cipro.R

R-script for data analysis in 24-hours mutants selection experiment. Generates barplot and performs descriptive statistics and ANOVA statistic analysis. The experiment was done to select for mutants resistant to ciprofloxacin and meropenem after 24 hours of the antibiotics exposure (in subinhibitory concentrations = 1/3 MIC) with or without 750 µL/mL of transferrin. K. pneumonia KP-nPDR#3 and KP-nPDR#4 strains. One example is presented.
