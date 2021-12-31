*AUTHOR: Regina Esinam Abotsi, Department of Molecular and Cell Biology, University of Cape Town, South Africa.

*DATE: March 22nd 2021
*Update: 17th November 2021

*PURPOSE: To conduct GLMM on bacterial carriage using viral load at baseline for all observation to retain maximum number of observations in the model
*Prerequisites: Generate or prepare data using the "Preparing_data_for_modeling17112021.do" script to generate "final_data_for_modeling17112021.csv" data to be used for this next step


log using log_17november2021modeling

clear
*import the data
import delimited "~/final_data_for_modeling17112021.csv", encoding(utf8) 

*------------------------------------------------------------------------------------------------------------
*------------------------NASOPHARYNGEAL SWABS----------------------------------------------------
*-------------------------------------------------------------------------------------------------------------

*Creating subset of data where NP were collected and cultured

*creating NP only data
tab np_available_cultured, missing 
*51 samples unavailable
drop if np_available_cultured=="no"
tab np_available_cultured, missing
*total of 874 observations remain



*NEW model _viral load at ENROLMENT
xi: melogit strep_np i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_np i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

xi: melogit staph_np i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_np i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

xi: melogit haem_np i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_np i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

xi: melogit morax_np i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_np i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

*The number of observations included in the above model is 857 out of 874 which is awesome



*Any bacterial carriage
tab anybact_int

*NEW model _viral load at ENROLMENT
xi: melogit i.anybact_int i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_np i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or





*------------------------------------------------------------------------------------------------------------
*------------------------SPUTUM----------------------------------------------------
*-------------------------------------------------------------------------------------------------------------
clear
*import the data again to obtain observations for which sputum is available but NP is not hence will have been excluded in by the previous analysis
import delimited ~/final_data_for_modeling17112021.csv, encoding(utf8) 

*Creating subset of data where SPUTUM were collected and cultured

*creating SPUTUM only data
tab sp_available_cultured, missing 
*81 samples unavailable
drop if sp_available_cultured=="no"
tab sp_available_cultured, missing
*total of 844 observations remain


*NEW model _viral load at ENROLMENT
xi: melogit strep_sp i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_sp i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

xi: melogit staph_sp i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_sp i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

xi: melogit haem_sp i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_sp i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

xi: melogit morax_sp i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_sp i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or

*The number of observations included in the above model is 827 out of 844 which is awesome



*Any bacterial carriage
tab anybact_int_spt 

gen anybact_int_sptb  = (anybact_int_spt =="yes")
tab anybact_int_sptb

*NEW model _viral load at ENROLMENT
xi: melogit i.anybact_int_sptb i.visit12m i.visit18m i.visit12m_azm i.visit18m_azm i.season1_sp i.normalwaz_reg i.agegpenrl_low i.sex i.site i.vlsupenrl || studyno:, intpoints(30) technique(nr) or


log close

