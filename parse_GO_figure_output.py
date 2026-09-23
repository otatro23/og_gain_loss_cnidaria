#!/usr/bin/env python3

# This script takes the deepfri output directory, deepfri filter value, and the go_figure_data.tsv from the deepfri R script.
# 1. read in the go_figure_data.tsv and create a pandas dataframe with two columns: category and GO term.
# 2. for each combination of category, event, and node, count the number of centroids that were annotated with a GO term in the category
# 3. store results from 2 in a dataframe, and write it as a tsv for stats and plotting in R

import argparse
import csv
import os

parser = argparse.ArgumentParser()
parser.add_argument("deepfri", help = "path to directory of deepfri output files")
parser.add_argument("--filter", "-f", type = int, default = 0.5,  help = "threshold to filter deepfri scores by")
parser.add_argument("go_figure_data", help = "go_figure_data.tsv from the R deepfri analysis script. cut for the description, category, and go members columns")
parser.add_argument("og_summary", help = "og_summary.tsv from the R deepfri analysis script")
args = parser.parse_args()

# structure: go: set of categories
go_cats = {}

with open(args.go_figure_data, "r") as go_data:
    for line in go_data:
        if line.startswith("description"):
            continue

        entries = line.rstrip().split("\t")
        gos = entries[2][1:-1].split(", ") # remove brackets and split to get a list of GO terms
        for go in gos:
            go_clean = go[1:-1] # remove the '' from the Go term
            go_cats.setdefault(go_clean, set())
            go_cats[go_clean].add(entries[1])


# structure: (category, node): [gain_count, loss_count]
counts = {}

categories = ["Transcription Regulator Activity", "Structural Molecule Activity", "Signal Transduction Activity", "Motor & Mechanical Activity", "Molecular Adaptors & Modulators", "Catalytic Activity", "Binding (Small Molecule, Ion, Cofactor)", "Transporter activity", "Specialized Structure", "Plasma Membrane & Cell Periphery", "Nucleus & Chromatin", "Membrane−Bound Organelles", "Extracellular Region & Matrix", "Endomembrane System & Vesicles", "Cytoskeleton & Cell Junctions", "Cytoplasm & Cytosol", "Metabolism & Biosynthesis", "Immune & Defense Response", "Gene Expression & Regulation", "Development & Morphogenesis", "Cellular Transport & Localization", "Cell Cycle & Division", "Cell Communication & Signaling", "Behavior, Movement, & Environmental Response"]

for category in categories:
    for node in range(98, 194):
        counts[(category, str(node))] = [0,0]

print(counts.keys())

# structure: og: set(categories)
og_cats = {}
for file in os.scandir(args.deepfri):
    print(file)
    with open(file, "r") as deepfri:
        for line in deepfri:
            if line.startswith("#"): # skip header
                continue

            if line.startswith("Protein"): # skip header
                continue

            entry = line.split(",")

            if float(entry[2]) < args.filter: # low score
                continue

            query = entry[0] # sequence ID
            og = query.split("_")[0]
            go = entry[1] # GO term

            og_cats.setdefault(og, set())
            if go in go_cats:
                for category in go_cats[go]:
                    og_cats[og].add(category)

with open(args.og_summary, "r") as og_sum:
    for line in og_sum:
        if line.startswith("og"): # skip header
            continue

        entries = line.rstrip().split("\t")

        og = entries[0]
        node = entries[1]
        event = 0 if entries[2] == "gain" else 1

        if node == 129 and event == 0:
            print("here")

        if og in og_cats:
            categories = og_cats[og]
            for category in categories:
                counts[(category, node)][event] += 1


with open("function_gl.tsv", "w") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["category", "node", "gain", "loss"])
    for ((category, node) , (gain, loss)) in counts.items():
        writer.writerow([category, node, gain, loss])
