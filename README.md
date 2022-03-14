# useful-bioinformatics-functions  

This repo contains useful functions for doing typical bioinformatics tasks that I have developed during my PhD, and generalised here to make them useful to others.  

# 1. helper-funs-to-query-api-opentargets-gh.R
  * I created these helper functions as wrappers for the https://genetics.opentargets.org/ API. 
  * There are a large multiude of query designs possible with this structure, but this particular use-case queries using a variantID.
  * It returns a dataframe containing some useful information related to it, namely the nearestGene.id which I originally designed this query to grab. 
  
# 2. plot-manhattan-from-MatrixEQTL-output-file.R  
  * This is a function that allows the user to plug in a raw [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) output file, and to output a Manhattan plot.  

# 3. apply_covariate_correction_to_gene_sample_matrix.py  
  * The function regress_covs_from_gene() takes in an expression matrix and a covariate matrix, and performs linear regression gene-wise to correct out the effects of covariates on the expression data. The outputs are an expression residuals matrix, and a table containing a Shapiro normality test P-value for each gene, enabling QC-ing of the residual data. An optional plot is also output to visualise the distribution of the Shapiro P-values.  
  