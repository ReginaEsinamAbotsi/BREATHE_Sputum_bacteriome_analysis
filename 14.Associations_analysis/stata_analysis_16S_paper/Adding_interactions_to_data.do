*AUTHOR: Regina Esinam Abotsi, Department of Molecular and Cell Biology, University of Cape Town, South Africa.

*DATE: December 20th 2021

*PURPOSE: To prepare data for LMME on 16S copies and Shannon diversity index Breathe 16S paper

log using log_Dec20th2021
clear

* Creating interaction terms 

import delimited /Users/reginaabotsi/Documents/Files/REGINA/MY_PROJECTS/MY_R_projects/Microbiome/Diversity_plots_Breathe_Sputum250620
21/Differential_abundancetesting02July_files/Nicolas/Analysis_compared_with_Yao/Other_analysis_for16spaper/Associations_analysis/stata_analysis_16S_paper/meta_edited24Nov2021No_NAsStata.csv, encoding(utf8) 




*------------------------------------------------------------------------------------------------------------
*------------------------Preparing data for all analyses----------------------------------------------------
*-------------------------------------------------------------------------------------------------------------

*Generating the interaction terms

gen visit12m = (visit=="Week 48")
*[creates a new variable coded 0 if not 12m visit and 1 if 12m visit]
tab visit12m

gen visit18m=(visit=="Week 72")
*[creates a new variable coded 0 if not 18m visit and 1 if 18m visit]
tab visit18m

gen azm = (trial_arm=="AZM") 
*[creates a new variable coded 0 if placebo and 1 if azm]
tab azm

gen visit12m_azm = visit12m*azm 
*[ creates a new variable coded 1 if azm and 12m visit and 0 otherwise]
tab visit12m_azm

gen visit18m_azm=visit18m*azm
*[ creates a new variable coded 1 if azm and 18m visit and 0 otherwise]
tab visit18m_azm


*Generate new variables to allow the right variable to be used as reference

*Need to derive vlsupenrl so supressed will be used as reference, (however vlsupenrl must be used since it has values for all observations)
tab vlenrlsup_baseline
gen vlsup_reg = "Unsuppresed" if vlenrlsup_baseline==">=1000 copies"
tab vlsup_reg, missing


replace vlsup_reg = "Suppresed" if vlenrlsup_baseline=="Supressed"
tab vlsup_reg, missing

*Need to derive normalwaz_reg so not underweight will be used as reference, 
tab normalwaz_baseline
gen normalwaz_reg = "underweight" if normalwaz_baseline=="<-2"
tab normalwaz_reg, missing


replace normalwaz_reg = "Not underweight" if normalwaz_baseline==">=-2"
tab normalwaz_reg, missing

*Need to derive normalhaz_reg so not stunted will be used as reference, 
tab normalhaz_baseline
gen normalhaz_reg = "stunted" if normalhaz_baseline=="<-2"
tab normalhaz_reg, missing


replace normalhaz_reg = "Not stunted" if normalhaz_baseline==">=-2"
tab normalhaz_reg, missing

export delimited using "/Users/reginaabotsi/Documents/Files/REGINA/MY_PROJECTS/MY_R_projects/Microbiome/Diversity_plots_Breathe_Sputu
m25062021/Differential_abundancetesting02July_files/Nicolas/Analysis_compared_with_Yao/Other_analysis_for16spaper/Associations_analys
is/stata_analysis_16S_paper/data_for_regression20thDec2021.csv", replace

