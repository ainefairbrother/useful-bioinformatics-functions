# useful-bioinformatics-functions
This repo contains useful functions that I have developed during my PhD. 

# 1. helper-funs-to-query-api-opentargets-gh.R
  * I created these helper functions as wrappers for the https://genetics.opentargets.org/ API. 
  * There are a large multiude of query designs possible with this structure, but this particular use-case queries using a variantID
  * and returns some useful information related to it, namely the nearestGene.id which I designed this query to fetch. 
