#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import os

# 2 positional arguments for the results and summary files from compare.py
# 1 optional argument for the output directory path
parser = argparse.ArgumentParser()
parser.add_argument("results", help = "directory with tabular data from compare.py, contains query, deepfri GO term, number of abstractions between DeepFRI GO term and eggNOG GO term, jaccard index")
parser.add_argument("--output_dir", "-o", default = "figures", help = "path to output directory for figures")
args = parser.parse_args()

# comparisons dictionary: key = filter value    values = ([abstractions], [jaccard indices])
comps = {}

# iterate through all of the results files
for tsv in os.scandir(args.results):
    try:
        abs = []
        jaccards = []

        # extract the filter value from the first line and the abstractions and jaccard indices from subsequent lines
        with open(tsv, "r") as file:
            for line in file:
                if line.startswith("#"):
                    filter = line.rstrip().split(" ")[7]
                    continue
                else:
                    values = line.rstrip().split("\t")
                    abs.append(int(values[2]))
                    jaccards.append(float(values[3]))

        # add to dictionary
        comps[filter] = (abs, jaccards)

    except IOError:
        print("Could not open {0}.".format(os.path))

# ------------------------------- plot figures -----------------------------------
try:
    os.mkdir(args.output_dir)
except FileExistsError:
    print("{0} directory already exists.".format(args.output_dir))

# barplot for how many entries there are based on the filter value
    # function: bar()
    # x = filters
    # height = length of abstractions or jaccards list
# box plots for abstractions and jaccard grouped by filter
    # function: boxplot()
    # x = list of lists of ints/floats
    # labels = filters
# histogram of deepfri 0.5 cutoff abstractions
    # function: bar()
    # x = list of abstraction values
    # height = number of annotations in each

# gathering correct objects for all plots:

# sort filters here and access the dictionary values by the filter.
filters = sorted(comps.keys()) # filter thresholds (as strings)
bar_heights = [] # number of annotations for each filter threshold
annotated_bar_heights = [] # number of annotations that had a match somewhere in the GO hierarchy to eggnog
abs_dists = [] # list of lists of abstractions
jaccard_dists = [] # list of lists of jaccard indices
dpi = 500

for filter in filters:
    bar_heights.append(len(comps[filter][0]))

    annotated = []
    for abs in comps[filter][0]:
        if abs != -1:
            annotated.append(abs)
    annotated_bar_heights.append(len(annotated))

    abs_dists.append(annotated)
    jaccard_dists.append(comps[filter][1])

# plot barplot
fig, axes = plt.subplots()
axes.bar(x = filters, height = bar_heights, color = "steelblue")
axes.bar(x = filters, height = annotated_bar_heights, color = "seagreen")
axes.set_xlabel("Filter threshold")
axes.set_ylabel("Total number of DeepFRI annotations")
fig.savefig(args.output_dir + "/barplot.png", format = "png", dpi = dpi)

# plot abstractions boxplot
fig, axes = plt.subplots()
axes.boxplot(x = abs_dists, labels = filters)
axes.set_xlabel("Filter threshold")
axes.set_ylabel("Abstractions between DeepFRI and eggNOG GO terms")
fig.savefig(args.output_dir + "/abs_boxplot.png", format = "png", dpi = dpi)

# plot jacccard boxplot
fig, axes = plt.subplots()
axes.boxplot(x = jaccard_dists, labels = filters)
axes.set_xlabel("Filter threshold")
axes.set_ylabel("Jaccard index")
fig.savefig(args.output_dir + "/jaccard_boxplot.png", format = "png", dpi = dpi)

# plot histogram
if "0.5" in filters:
    abs = {}

    for abs_value in comps["0.5"][0]:
        abs.setdefault(abs_value, 0)
        abs[abs_value] += 1

    labels = sorted(abs.keys())
    heights = []
    for label in labels:
        heights.append(abs[label])

    fig, axes = plt.subplots()
    axes.bar(x = labels, height = heights, color = "steelblue")
    axes.set_xlabel("Abstractions between DeepFRI and eggNOG GO terms")
    axes.set_ylabel("Number of DeepFRI GO terms")
    fig.savefig(args.output_dir + "/barplot0.5.png", format = "png", dpi = dpi)
