# useful-bioinformatics-functions
This repo contains useful functions for bioinformatics-related tasks that I have developed during my PhD. 

# 1. helper-funs-to-query-api-opentargets-gh.R
  * I created these helper functions as wrappers for the https://genetics.opentargets.org/ API. 
  * There are a large multiude of query designs possible with this structure, but this particular use-case queries using a variantID.
  * It returns a dataframe containing some useful information related to it, namely the nearestGene.id which I originally designed this query to grab. 
  
# 2. plot-manhattan-from-MatrixEQTL-output-file.R  
  * This is a function that allows the user to plug in a raw [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) output file, and to output a Manhattan plot. 
