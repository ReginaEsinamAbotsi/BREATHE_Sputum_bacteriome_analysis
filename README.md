# BREATHE_Sputum_bacteriome_analysis
Details all analysis related to the BREATHE Sputum bacteriome analysis

This github repository houses R notebooks of codes used in the manuscript Abotsi RE. The impact of long-term azithromycin treatment on the sputum bacteriome of African children with HIV-associated chronic lung disease: an analysis from the double-blind, placebo-controlled BREATHE trial Microbiome .... https://doi....


The following README gives an overview of the overall structure of the repository, and important notes on how to run the scripts.

# Directory tree for repository
BREATHE_Sputum_bacteriome_published_2021



      1.Microbiome_published_material -- contains manuscript (pdf) and supplementary material (pdf, xlsx) published in BMC Microbiome

      2.Preprocessing_QC_generating_phyloseq -- contains scripts used for sequences processing to generating first phyloseq object using TADA - Targeted Amplicon Diversity Analysis using DADA2, implemented in Nextflow found here https://github.com/h3abionet/TADA

      3.Decontamination_steps -- houses R scripts used to check extraction efficiency,reproducibility between and within run, in silico decontamination etc 

      4.Generating_final_phyloseq_object-- houses R scripts used to generate final phyloseq object made up of the final ASVs, samples and metadata retained after decontamination step

      5.16S_copies_Alpha_diversity-- houses R scripts used to determine the 16S copies and alpha diversity results and figures

      6.Rel.abun.plots.microeco-- houses R scripts used to generate the barplots of relative abundance and group mean barplots

      7.Within-patient_change_in_beta_diversity-- houses R scripts used to generate the within-patient change in beta diversity results and figures

      8.Between_group_beta_diversity-- houses R scripts used to generate the between-group beta diversity results and figures

      9.PCoA-- houses R scripts used to generate the PCoA plots to visualise beta diversity

      10.ANcom_results_script-- houses R scripts used to determine the differentially abundant taxa using Ancom2 

      11.Simper_analysis-- houses R scripts used to conduct SIMPER analysis

      12.Associations-- houses R scripts used to determine associations between clinical and socio-demographic factors and alpha diversity

      13.Final_DA_heatmaps-- houses R scripts used to determine the differentially abundant taxa using the other 9 methods and generating heatmap from the results 

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
