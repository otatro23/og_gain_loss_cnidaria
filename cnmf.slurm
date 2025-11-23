#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=cnmf
#SBATCH --output=cnmf.log
#SBATCH --mail-user=Olivia.Tatro@unh.edu   
#SBATCH --mail-type=END,FAIL  
##SBATCH --open-mode=append

# command line inputs:
# $1 = output dir (genus/cnmf)
# $2 = input counts dir (genus/<counts>)
# $3 = genus for messages

# --------------------------- prep input files ------------------------------
#awk '{print $1"\t"$1}' $2/features.tsv > $2features.2col.tsv
#awk '{print $1"\t"$2"\tGene Expression"}' $2/features.2col.tsv > $2/features.tsv
#gzip $2/features.tsv
#gzip $2/barcodes.tsv
#gzip $2/matrix.mtx

# ---------------------------- activate environment -----------------------------
module purge
module load anaconda/colsa
#conda create --name cnmf_geps --clone template
conda activate cnmf_geps
#conda install python=3.12 
#pip install cnmf
#pip install harmonypy
#pip install scikit-misc

# -------------------------------- Running CNMF ----------------------------------
# the following steps follow the tutorial at: https://github.com/dylkot/cNMF/blob/main/Stepwise_Guide.md

# Step 1: normalize the input matrix and prepare the run parameters
echo "Step 1: normalize the input matrix and prepare the run parameters"
cnmf prepare --output-dir $1 -c $2/matrix.mtx.gz -k 5 10 15 20 25 30 35 40 45 --n-iter 100 --seed 14 --numgenes 2000 --tpm $2/matrix.mtx.gz
echo "Step 1 completed."

# Step 2: factorize the matrix
echo "Step 2: factorize the matrix"
cnmf factorize --output-dir $1
echo "Step 2 completed."

# Step 3: combine the individual spectra results files for each K into a merged file
echo "Step 3: combine the individual spectra results files for each K into a merged file"
cnmf combine --output-dir $1
echo "Step 3 completed."

# Step 4: select an optimal K by considering the trade-off between stability and error
echo "Step 4: select an optimal K by considering the trade-off between stability and error"
cnmf k_selection_plot --output-dir $1
echo "Step 4 completed."

# Step 5: obtain consensus estimates for the programs and their usages at the desired value of K
# need to change the --components to the optimal k value as determined by the previous step. Documentation on the github says to choose the k value where there is a local maximum for stability.
#echo "Step 5: obtain consensus estimates for the programs and their usages at the desired value of K"
#cnmf consensus --output-dir $1 --components 7 --local-density-threshold 2.0 --show-clustering
#echo "Step 5 completed."
# they recomended running this a second time with a lower --local-density-threshold based on the histogram. We want to choose a threshold that eliminates outliers. 
#cnmf consensus --output-dir $1 --components 7 --local-density-threshold 0.1 --show-clustering
#echo "Step 5 pt2 completed."


# ---------------------------------- gep mapping ------------------------------------
