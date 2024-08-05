# Analysis code for "Prenatal exposures and cell type proportions are main drivers of FKBP5 DNA methylation in maltreated and non-maltreated children"

Author: Vera N. Karlbauer  
Contact: vera_karlbauer@psych.mpg.de  
Date created: 25-07-2024  

[![DOI](https://zenodo.org/badge/838393393.svg)](https://zenodo.org/badge/latestdoi/838393393)

## Contains: 
- 01_functions: core functions for statistical analysis  
- 02_analysis: scripts for statistical analysis, including figure & table generation 

### Scripts included in 02_analysis:
- 01_main_table_1_sample_demographics.Rmd: Sample demographics and distribution of key variables for Table 1.  
- 02_main_figure_1_track_plot_sequenced_regions.Rmd: Creates track plot of all sequenced TBS amplicons on FKBP5 (Figure 1) and shows overlap with GREs and CTCF binding sites. 
- 03_appendix_b_table_b2_mean_dnam_per_cpg.Rmd: Creates table with mean DNAm for each CpG at T0 timepoint (Supplementary Figure B2).  
- 04_main_variance_partitioning.Rmd: Variance partitioning analysis to determine drivers of FKBP5 methylation, creates Figure 2.  
- 05_main_diffmeth_maltreatment.Rmd: Tests for effects of childhood maltreatment on FKBP5 DNAm using linear models, creates Supplementary Table B3.  
- 06_main_diffmeth_diagnosis.Rmd: Tests for effects of psychiatric diagnosis on FKBP5 DNAm using linear models, creates Supplementary Table B4.   
- 07_main_diffmeth_genotype.Rmd: Tests for effects of rs1360780 genotype on FKBP5 DNAm using linear models, tests for maltreatment-genotype interaction on FKBP5 DNAm, creates Supplementary Table B5 and B6, creates Figure 3.   
- 08_main_diffmeth_gc_exposure.Rmd: Tests for effects of epigenetic glucocortioid exposure score on FKBP5 DNAm using linear models, tests for glucocortioid-genotype interaction on FKBP5 DNAm, creates Supplementary Table B7 and B8, creates Figure 4.   
- 09_main_gc_associations.Rmd: Analyzes associations of prenatal glucocorticoid exposure with different maternal and childhood phenotypes, creates Table 2. 
- 10_main_correlation_structure.Rmd: Examines correlations between exposure scores, estimated cell type proportions, and biodata, creates Figure 5. 
- 11_appendix_a_exploration_exposure_scores.Rmd: Examines interactions between maltreatment and epigenetic prenatal alcohol exposure score on epigenetic glucocorticoid exposure score, creates Supplementary Figure A1. 

## Folders referenced in the scripts, but not pushed to the repository:  
- 03_results: output from 02_analysis in table form saved as .txt files  
- 04_figures: figures generated from 02_analysis saved as .png files  
- 05_tables: tables in final publication format   
