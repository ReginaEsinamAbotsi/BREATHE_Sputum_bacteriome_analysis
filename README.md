# BREATHE_Sputum_bacteriome_analysis
Details all analysis related to the BREATHE Sputum bacteriome analysis

This github repository houses R notebooks of codes used in the manuscript Abotsi RE. The impact of long-term azithromycin treatment on the sputum bacteriome of African children with HIV-associated chronic lung disease: an analysis from the double-blind, placebo-controlled BREATHE trial Microbiome .... https://doi....


The following README gives an overview of the overall structure of the repository, and important notes on how to run the scripts.

# Directory tree for repository
BREATHE_Sputum_bacteriome_published_2021



      1.Microbiome_published_material -- contains manuscript (pdf) and supplementary material (pdf, xlsx) published in BMC Microbiome

      2.Preprocessing_QC_generating_phyloseq -- contains scripts aused for sequences processing to generating first phyloseq object using TADA - Targeted Amplicon Diversity Analysis using DADA2, implemented in Nextflow found here https://github.com/h3abionet/TADA. The output from the pipeline is also included.

      3.Decontamination_steps -- houses R scripts used to check extraction efficiency,reproducibility between and within run, in silico decontamination etc 

      4.Generating_final_phyloseq_object-- houses R scripts used to generate final phyloseq object made up of the final ASVs, samples and metadata retained after decontamination step

      5.16S_copies_Alpha_diversity-- houses R scripts used to determine the 16S copies and alpha diversity results and figures

     
      6.Within-patient_change_in_beta_diversity-- houses R scripts used to generate the within-patient change in beta diversity results and figures

      7.Between_group_beta_diversity-- houses R scripts used to generate the between-group beta diversity results and figures

      8.PCoA-- houses R scripts used to generate the PCoA plots to visualise beta diversity
      
      9.PERMANOVA-- houses R scripts used to conduct PERMANOVA tests 

      10.Rel.abun.plots.microeco-- houses R scripts used to generate the barplots of relative abundance and group mean barplots
      
      11a.DA_analysis_AZM_only_Baseline_Week48 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the AZM arm only from baseline and 48 week visit using 10 methods
      
      11b.DA_analysis_AZM_only_Baseline_Week72 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the AZM arm only from baseline and 72 week visit using 10 methods
      
      11c.DA_analysis_AZM_only_Week48_Week72 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the AZM arm only from 48 and 72 week visit using 10 methods
      
      11d.DA_Analysis_AZM_Placebo_Baseline -- houses R scripts used to conduct ddifferentially abundant taxa analysis on samples from the AZM and Placebo arm at baseline using 10 methods
      
      11e.DA_Analysis_AZM_Placebo_Week48 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the AZM and Placebo arm at week 48 using 10 methods
      
      11f.DA_Analysis_AZM_Placebo_Week72 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the AZM and Placebo arm at week 72 using 10 methods
     
     11g.DA_analysis_Placebo_only_Baseline_Week48 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the Placebo arm only from baseline and 48 week visit using 10 methods
     
     11h.DA_analysis_Placebo_only_Baseline_Week72 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the Placebo arm only from baseline and 72 week visit using 10 methods
     
     11i.DA_analysis_Placebo_only_Week48_Week72 -- houses R scripts used to conduct differentially abundant taxa analysis on samples from the Placebo arm only from 48 and 72 week visit using 10 methods
     
    11j. DA_Analysis_Supplementary_Tables -- houses supplementary tables from the differentially analysis test on samples from the AZM arm only (Baseline and 48 weeks, and 48 and 72 weeks) and from AZM and Placebo arms at 48 weeks using 10 methods

      12.ANCOM2_Analysis_AZM_only_allvisits-- houses R scripts used to determine the differentially abundant taxa in AZM arm at all visits using Ancom2 

      13.Simper_analysis-- houses R scripts used to conduct SIMPER analysis

      14.Associations-- houses R scripts used to determine associations between clinical and socio-demographic factors and alpha diversity

      
          |-- AZM_Placebo -- contains scripts used to run differentially abundant taxa analysis on AZM and Placebo samples 

                |-- Baseline -- contains scripts used to run differentially abundant taxa analysis on AZM and Placebo samples at baseline

                |-- 12m -- contains scripts used to run differentially abundant taxa analysis on AZM and Placebo samples at 48 weeks

                |-- 18m -- contains scripts used to run differentially abundant taxa analysis on AZM and Placebo samples at 72 weeks

         |-- AZM_only -- contains scripts used to run differentially abundant taxa analysis on only AZM samples 

               |-- AZM012 -- contains scripts used to run differentially abundant taxa analysis on AZM baseline and 48 weeks 

               |-- AZM1218 -- contains scripts used to run differentially abundant taxa analysis on AZM 48 and 72 weeks samples

               |-- AZM018 -- contains scripts used to run differentially abundant taxa analysis on AZM baseline and 72 weeks samples

        |-- Placebo_only -- contains scripts used to run differentially abundant taxa analysis on only Placebo samples 

               |-- Placebo012 -- contains scripts used to run differentially abundant taxa analysis on Placebo baseline and 48 weeks 

               |-- Placebo1218 -- contains scripts used to run differentially abundant taxa analysis on Placebo 48 and 72 weeks samples

               |-- Placebo018 -- contains scripts used to run differentially abundant taxa analysis on Placebo baseline and 72 weeks samples
               
       15.Final_DA_heatmaps-- houses R scripts used to determine the differentially abundant taxa using the other 9 methods and generating heatmap from the results 
       
       Data dictionary -- houses excel files of variables and their interpretation 
      
Figures_made_inppt-- houses differentially abundant taxa analysis results summary in figures on a Microsoft Powerpoint file. Also contains figure on number of samples collected and analysed


