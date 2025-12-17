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

# goal:
# I want to get a dictionary with keys acropora_1... for each gep. The values will be a dictionary with OG:expression z score, but for now, it will simply be gene:expression z score

geps = {}

for file in os.scandir("test"):
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

    for og, genes in ogs.items():
        spectra[og] = 0
        for gene in genes:
            spectra[og] += spectra[gene]

    print(spectra)

    continue

    for gep_index in range(len(spectra)):

        for gene, value in spectra.iloc[gep_index].items():
            #print(f"Index: {index}, Value: {value}")
            if gene.startswith("Unnamed"):
                gep_name = taxon + "_" + str(int(value))
                print(gep_name)
                geps[gep_name] = {}
            else:
                try:
                    og = ogmap[taxon + "-" + gene]
                    print(og)
                    geps[gep_name].setdefault(og, 0)
                    geps[gep_name][og] += value
                except KeyError as key:
                    print(f"Key {key} not found in ogmap.")


