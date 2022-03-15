# setup file to generate a conda environment to run apply_covariate_correction_to_gene_sample_matrix.py  

# initiate conda environment containing python v3.8
conda create -n py38 python=3.8
# activate the env
conda activate py38
# install packages - sys, re, os should be installed already
conda install -c conda-forge regex
conda install -c conda-forge matplotlib
conda install -c anaconda seaborn
conda install -c anaconda pandas
conda install -c anaconda numpy
conda install -c anaconda scikit-learn
conda install -c anaconda scipy
conda install -c conda-forge fire
echo "Environment py38 ready to run apply_covariate_correction_to_gene_sample_matrix.py"


