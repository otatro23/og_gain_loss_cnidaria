#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=diamond
#SBATCH --output=diamond.log
#SBATCH --mail-user=Olivia.Tatro@unh.edu   
#SBATCH --mail-type=END,FAIL  

# this script will blast the protein fasta from the sc data to the orthofinder input fasta
# requires: data directory with genus.fa and genus.orthofinder.fa
# all other prep to get this information needs to be done within each genus directory

module purge
module load linuxbrew/colsa
conda activate diamond-2.1.9

#diamond --help

# clytia has a separate script in its directory because I used diamond blastx instead of blastp
genera="aurelia acropora turritopsis"
for genus in $genera ; do
    echo "running blast for: $genus"
    dir="${genus}/data/"
    diamond makedb --in $dir/$genus.orthofinder.fa  -d $dir/${genus}_db
    diamond blastp --query $dir/$genus.fa --db $dir/${genus}_db --out $dir/${genus}.diamond.tsv --outfmt 6 --evalue 1e-5 --threads 24
done


