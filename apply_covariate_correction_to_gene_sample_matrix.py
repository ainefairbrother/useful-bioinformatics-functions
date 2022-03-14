###################################################################
# Title: Functions to regress covariates from an expression matrix
# Author: Aine Fairbrother-Browne
# Date: 03/22
#
# Usage:
#
# Inputs:
# The correct_gene_by_sample_matrix() function takes in a gene x sample
# expression matrix, mat_file, and a covariate x sample covariate file,
# covs_file. It regresses the covariates from the expression matrix
# using a linear model, applying it gene-wise. 
#
# Outputs:
# 1. A matrix of the residual values
# 2. A table of Shapiro normality test P-values 
# 3. A distribution plot to show the spread of Shapiro P-values
###################################################################

def regress_covs_from_gene(y):

  """
  Function to regress out covariates from each gene - applies over columns (genes)
  Inserts residuals into predefined residual dataframe with correct row (sample id) and col (gene) labels
  """
  gene_name = y.name
  
  print(gene_name)

  # predictors
  X = covs

  # response var
  y = np.array(y)

  # get list of locations where gene has nan values
  nan_locs = np.where(np.isnan(y))[0]

  # masking NAs
  finiteYmask = np.isfinite(y)
  
  # cleaning NAs from predictors and response vars
  Yclean = y[finiteYmask]
  Xclean = X[finiteYmask]
  
  # fit the model 
  lm.fit(Xclean, Yclean)

  # make predictions from covs - i.e. expected TPMs
  predictions = lm.predict(Xclean)
  
  # reinsert removed nan values one-by-one to make y and predictions the same length
  for i in range(0, len(nan_locs)):
    predictions = np.insert(predictions, nan_locs[i], np.nan)

  # residuals = observations - predictions, normalising so they're normally distributed
  residuals = y - predictions

  # adding residuals to residual df
  residual_df[gene_name] = residuals
  
  # cleaning residuals of NaNs for shapiro fn. 
  shapiro_resids = residuals[np.logical_not(np.isnan(residuals))]

  # get shapiro pvals and add to shapiro_df
  shapiro_pval = shapiro(shapiro_resids)[1]
  shapiro_df.loc[gene_name, 'shapiro_pval'] = shapiro_pval
  
def correct_gene_by_sample_matrix(mat_file="", covs_file="", plot_shapiro_dist=True):

  """
  :param mat: file path of a transcriptomic gene (rows) x sample (cols) matrix, where the values are counts, TPMs, FPKMs, RPKMs, etc (csv)
  :param covs: file path of a covariate (rows) by sample (cols) matrix containing the covariates that you want to correct for. These must be numeric, for instance sex should be encoded as 0|1. (csv)
  :param: plot_shapiro_dist
  :return: residuals df
  
  """
  
  # import libs
  import pandas as pd
  import numpy as np
  from sklearn.linear_model import LinearRegression
  from scipy.stats import shapiro
  import sys
  import os
  import re
  
  # import tpm file
  print("1. Importing mat_file.")
  mat = pd.read_csv(open(mat_file), 
    encoding="utf-8", 
    engine='python', 
    index_col=0, 
    header=0)

  print("mat dimensions=",mat.shape)

  print("2. Importing covs_file.")
  # import covariate/metadata file
  covs = pd.read_csv(open(covs_file), 
    encoding="utf-8", 
    engine='python', 
    index_col=0, 
    header=0)

  print("covs dimensions=",covs.shape)

  print("3. Filtering mat and covs for column names in common.")

  # getting columns (sample IDs) in common between covs and mat
  intersection_cols = covs.columns & mat.columns
  intersection_cols
  print(len(intersection_cols), "column names in common between mat and covs.")

  covs = covs[intersection_cols]
  mat = mat[intersection_cols]

  print("After filtering, the dimensions of covs =", covs.shape, "and dimensions of mat =", mat.shape)

  # fill encoded covs with mean - fills with most common value for the variable in question
  # allows the mlr fn. to run, and should make very little difference to the correction 
  covs = covs.T
  covs = covs.fillna(covs.mean())

  # make samples = rows, genes = cols
  mat = mat.T

  # remove genes (cols) that are all 0
  mat = mat.loc[:, ~mat.eq(0).all()]
  # remove genes (cols) that are all NA
  mat = mat.dropna(axis=1)

  print("4. Checking for duplicate columns")
  if(mat.loc[:,mat.columns.duplicated()].shape[1]>0):
    print(mat.loc[:,mat.columns.duplicated()].shape[1], "duplicate columns detected, these will be removed.")
    if(mat.loc[:,mat.columns.duplicated()].shape[1]==0):
      print("No duplicate columns detected.")

  # remove dup cols in mat
  mat = mat.loc[:,~mat.columns.duplicated()]

  # remove index name
  mat = mat.rename_axis(None, axis=1)
  covs = covs.rename_axis(None, axis=1)

  # performing mlr across genes 

  # defining empty array to hold residual tpm values
  residual_df = pd.DataFrame(np.zeros(shape=(len(mat.index.values), len(mat.columns.values))),
    index=mat.index,
    columns=mat.columns)
  
  # defining empty array to hold shapiro pvals
  shapiro_df = pd.DataFrame(np.zeros(shape=(len(mat.columns.values), 1)),
   index=mat.columns,
   columns=['shapiro_pval'])

  # initiating the lm
  print("5. Initialising the linear regression model")
  lm = LinearRegression()

  # applying regress_covs_from_gene fn. across genes
  # for 160 samples (rows) and 15K genes (cols), this step takes 41 seconds
  print("6. Applying linear regression across", mat.shape[1], "columns (genes)")
  
  mat.apply(regress_covs_from_gene, axis=0);
  
  print("Regression done, residuals generated.")

  # find % genes with normally distributed residuals
  sig_norm_gene_count = (int(shapiro_df[shapiro_df <= 0.05].count()) / len(shapiro_df.index.values)) * 100
  print(round(sig_norm_gene_count, 2), "% of gene residuals have a significantly normal distribution")

  # generate outfile names
  print("7. Writing output")
  residual_out_file = mat_file.replace(".csv", "_residuals.csv") 
  shapiro_out_file = mat_file.replace(".csv", "_shapiro_pvals.csv") 

  # export output files
  residual_df.to_csv(residual_out_file, index=True, header=True)
  shapiro_df.to_csv(shapiro_out_file, index=True, header=True)

  if plot_shapiro_dist==True:
    
    print("8. Generating and writing out a plot to show the distribution of shapiro p-values. The vertical red line indicates P=0.05.")
    
    import matplotlib
    matplotlib.use('Agg')
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # get plot outfile name
    plot_out_file = mat_file.replace(".csv", "_shapiro_pvals.png")
    
    # plotting a density plot for the shapiro(residual) p-values

    plt.figure()
    sns.distplot(shapiro_df['shapiro_pval'], hist=True)
    plt.xlabel('Shapiro pvals') 
    plt.axvline(x=0.05, ymax=10, c='red', alpha=0.5)
    plt.show()
    plt.savefig(plot_out_file)





