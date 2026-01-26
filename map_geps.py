#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import time

# ------------------------------- argparse -------------------------------
parser = argparse.ArgumentParser(prog = "map_geps.py", description = "Reads GEP spectra from cNMF from scRNAseq from multiple taxa and compares the GEPs to eachother.")
parser.add_argument("input_dir", help = "directory with tpm matrices for each taxon. File names must start with \"taxon_\".")
parser.add_argument("ogmap", help = "mapping file with columns: \"og\" and \"sc_ID\", which are the gene IDs from each taxon. The genes must be named in the format: \"taxon-geneID\", e.g. acropora-evm.TU.Chr4.452")
parser.add_argument("--output", "-o", default = "map_geps_output", help = "output directory name")
args = parser.parse_args()

# ------------------------------------ read in ogmap -------------------------------------------
ogmap = {} # gene:og
try:
    with open(args.ogmap, "r") as ogmap_file:
        for line in ogmap_file:
            if line.startswith("og"):
                continue
            cols = line.rstrip().split("\t")
            ogmap[cols[1]] = cols[0]
except IOError:
    print(f"Error reading {args.ogmap}.")

# ------------------------------------- read in gene spectra files  -------------------------------------

# dictionary of taxon:gep x og dataframe
gep_spectra = {}

for file in os.scandir(args.input_dir):
    taxon = file.name.split("_")[0]
    print("Taxon:", taxon)

    spectra = pd.read_table(file.path)

    genes = list(spectra.columns)
    genes.pop(0)

    ogs = {}

    # count how many genes are missing from the ogmap
    missed_genes = 0
    total_genes = 0

    for gene in genes:
        total_genes += 1
        try:
            og = ogmap[taxon + "-" + gene]
            ogs.setdefault(og, [])
            ogs[og].append(gene)
        except KeyError as key:
            #print(f"Gene: {key} is not in ogmap. Skipping.")
            missed_genes += 1
    print(f"missed genes: {missed_genes}    total genes: {total_genes}")

    og_spectra = {} # {Orthogroup:gene expression}

    # adds the expression of each gene in an orthogroup to the orthogroup column in spectra
    for og, genes in ogs.items():
        # set the og entry in og_spectra to the first gene
        og_spectra[og] = list(spectra[genes[0]])
        genes.pop(0)
        # iterate through remaining genes and add their expression values to the list at og_spectra[og]
        if len(genes) > 0:
            for gene in genes:
                exp_list = list(spectra[gene])
                for index, exp in enumerate(exp_list):
                    og_spectra[og][index] += exp

    gep_spectra[taxon] = pd.DataFrame(og_spectra)

# ------------------ calculate orthogroup proportion in each gep ---------------
# takes a df and divides each value by the row total
def convert_tpm_prop(df):
    row_sums = df.sum(axis = 1)
    df_prop = df.div(row_sums, axis = 0)
    return df_prop

# {taxon : proportion df}
# the dataframe is the result of adding all of the tpms in each gep and then dividing each value by the total. Basically the proportion of expression of each gene/OG of all gene expression in the gep
gep_prop = {}

for taxon, df in gep_spectra.items():
    gep_prop[taxon] = convert_tpm_prop(df)


# ---------------------------- compare geps --------------------------------
hdf = gep_prop["hydra"]
ndf = gep_prop["nematostella"]

# each gep is a series with names being orthogroups and values beign the proportion of total gep expression it represents
def compare_geps(gep1, gep2):
    ogs = gep1.index
    prop_similar = 0 # the amount of similarity there is in the gene expression between 2 geps
    for og in ogs:
        exp = [gep1[og], gep2.get(og, 0)]
        prop_similar += min(exp)
    return(prop_similar)

# prop_df1 and 2 are the og expression proportions
# this function compares each from 1 species to each gep from the other. 
def compare_taxa(taxon1, prop_df1, taxon2, prop_df2):
    similarities = {} # the keys will be the geps from taxon1, the values are dictionaries with each gep from taxon2 as the keys and the similarities as values
    for gep1 in range(len(prop_df1)):
        gep1_name = taxon1 + "_" + str(gep1)
        similarities[gep1_name] = {}

        for gep2 in range(len(prop_df2)):
            gep2_name = taxon2 + "_" + str(gep2)
            similarity = compare_geps(prop_df1.iloc[gep1], prop_df2.iloc[gep2])
            similarities[gep1_name][gep2_name] = similarity

    # convert to similarity df
    similarities_df = pd.DataFrame(similarities)
    return(similarities_df)

# compare all geps between taxa
os.mkdir(args.output)

for taxon1, prop_df1 in gep_prop.items():
    for taxon2, prop_df2 in gep_prop.items():
        if taxon1 <= taxon2: # ensures each comparison only happens once
            continue
        print(f"{taxon1} vs {taxon2} comparision")
        similarities_df = compare_taxa(taxon1, prop_df1, taxon2, prop_df2)
        #print(similarities_df.max(axis = 1))
        similarities_df.to_csv(f"{args.output}/{taxon1}_{taxon2}.tsv", sep = "\t")
