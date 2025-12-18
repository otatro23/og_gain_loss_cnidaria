#!/usr/bin/env python3


### IN PROGRESS ###

import os
import pandas as pd

# ------------------------------------ read in ogmap -------------------------------------------
ogmap = {} # gene:og
try:
    with open("ogmap.tsv", "r") as ogmap_file:
        for line in ogmap_file:
            if line.startswith("og"):
                continue
            cols = line.rstrip().split("\t")
            ogmap[cols[1]] = cols[0]
except IOError:
    print("Error reading ogmap.tsv")

# ------------------------------------- read in gene spectra files  -------------------------------------

# dictionary of taxon:gep x og dataframe
gep_spectra = {}

for file in os.scandir("z_score_gene_spectra"):
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
            og = "orthogroup:" + ogmap[taxon + "-" + gene]
            ogs.setdefault(og, [])
            ogs[og].append(gene)
        except KeyError as key:
            #print(f"Gene: {key} is not in ogmap. Skipping.")
            missed_genes += 1
    print(f"missed genes: {missed_genes}    total genes: {total_genes}")

    og_spectra = {}

    # adds the expression of each gene in an orthogroup to the orthogroup column in spectra
    for og, genes in ogs.items():
        #spectra[og] = 0
        for gene in genes:
            #spectra[og] += spectra[gene]

    new_df = {}

    for col in spectra.columns:
        if col.startswith("orthogroup:"):
            # removes the "orthogroup:" from the OG columns
            #spectra[col.split(":")[1]] = spectra[col]
            new_df[col.split(":")[1]] = spectra[col]

        elif col.startswith("Unnamed"):
            #renames gep column
            #spectra["gep"] = spectra[col]
            new_df["gep"] = spectra[col]
        else:
            #spectra.pop(col)
            pass

    print(new_df)

    gep_spectra[taxon] = pd.DataFrame(new_df)


# input: taxon name and its gep-spectra df
# output: dictionary of geps and sets of upregulated orthogroups
def get_ogs(taxon, df):
    high_ogs = {} # keys: gep     values: set of OGs upregulated

    for index in range(len(df)):
        row = df.iloc[index]
        gep = taxon + "_" + str(int(row["gep"]))
        og_exp = row.drop("gep")
        og_exp = og_exp[og_exp > 0] # filter for expression zscore over 0 (upregulated OGs)
        high_ogs[gep] = set(og_exp.index) # get the ogs that got by the filtering step

    return high_ogs

ogs = {} # key: taxon     {value: set(upregulated ogs)}

for taxon, df in gep_spectra.items():
    ogs[taxon] = get_ogs(taxon, df)

def jaccard(set1, set2):
    return len(set1.intersection(set2)) / len(set1.union(set2))

# input dictionaries of gep:set(og) for 2 taxa
# output: matrix of jaccard indices of OG sets. Rows are taxon1; cols are taxon2.
def compare(geps_ogs1, geps_ogs2):
    jaccards = pd.DataFrame(index = geps_ogs1.keys(), columns = geps_ogs2.keys())
    print(jaccards)
    for gep1, og_set1 in geps_ogs1.items():
        for gep2, og_set2 in geps_ogs2.items():
            jaccards.loc[gep1, gep2] = jaccard(og_set1, og_set2)
    print(jaccards)

compare(ogs["aurelia"], ogs["nematostella"])



