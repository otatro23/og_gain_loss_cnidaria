#!/usr/bin/env python3

# command: ./compare.py go_basic.obo eggnog.tsv /mnt/home/plachetzki/omt1027/surf/97_taxa/orthofinder_run/deepfri/deepfri_output/test

# ---------------------------------------- import -------------------------------------------
import argparse
import os

# --------------------------------------- arguments -----------------------------------------
# user must put in the go.obo file and the output of both eggnog and deepfri
parser = argparse.ArgumentParser()
parser.add_argument("go", type = str, help = "path to go.obo file")
parser.add_argument("eggnog", type = str, help = "path to eggnog output")
parser.add_argument("deepfri", type = str, help = "path to deepfri output directory")
parser.add_argument("--deepfri_filter", "-d", type = float, default = 0.5, help = "score to filter deepfri data by")
parser.add_argument("--out", "-o", type = str, default = "compare.out.tsv", help = "name of output file")
parser.add_argument("--abstractions", "-a", type = int, default = 2, help = "number of \"is_a\" relationship levels to consider")
args = parser.parse_args()

# -------------------------------------- parse go.obo ----------------------------------------

# class that contains go terms and instance attributes for their "is_a" relationships, whether they are obsolete or not, and their replacement if they have one
class Go:
    def __init__(self, id):
        self.id = id
        self.is_obsolete = False
        self.replacement = None
        self.parents = []
        self.alt_ids = []

    def add_parent(self, go):
        self.parents.append(go)

    def __str__(self):
        string = self.id

        if len(self.alt_id) > 0:
            for id in self.alt_ids:
                string += "/" + id

        if self.is_obsolete:
            string += " is obsolete."
            if self.replacement is not None:
                string += " Replaced by: " + self.replacement
        elif self.parents != []:
            string += " has parent: " + str(self.parents)
        return(string)


# dictionary of go objects
gos = {}

# read in the go terms
with open(args.go, "r") as go_obo:
    for line in go_obo:

        # create a new Go object with go id. the len() part is to deal with a few lines that had id: but not a GO term. all GO terms will be the same length
        if line.startswith("id:") and len(line.rstrip()) == 14:
            go_term = line.rstrip().split(" ")[1]
            gos[go_term] = Go(go_term)

        # adds go terms with "is_a" relationship to the .is_a attribute of the Go instance
        if line.startswith("is_a"):
            is_a = line.split(" ")[1]
            gos[go_term].add_parent(is_a)

        # changes the Go instance attribute of is_obsolete to True if it is
        if line.startswith("is_obsolete"):
            gos[go_term].is_obsolete = True

        # adds in replacement go if it is obsolete and it has a replacement
        if line.startswith("replaced_by"):
            rep = line.rstrip().split(" ")[1]
            gos[go_term].replacement = rep

        if line.startswith("alt_id"):
            gos[go_term].alt_ids.append(line.rstrip().split(" ")[1])

# create new go object for alternate ids. If one program has the normal id and one has the alt id, both will end up pointing to the normal id
alt_gos = {}
for go in gos.values():
    for alt in go.alt_ids:
        alt_gos[alt] = go
# adding the new gos back into the gos dictionary
for query, go in alt_gos.items():
    gos[query] = go


# format: key = query sequence ID. Value = list(set(eggnog annotations), set(deepfri annotations))
annotations = {}

# read in eggnog data by sequence. add to annotations dictionary
with open(args.eggnog, "r") as eggnog:
    for eggnog_line in eggnog:
        entry = eggnog_line.split("\t")
        query = entry[0].split(":")[1]
        eggnog_gos = entry[6].split(",")
        annotations[query] = [set(eggnog_gos), set()]

# read in deepfri data
for file in os.scandir(args.deepfri):
    #print(file)
    with open(file, "r") as deepfri:
        for line in deepfri:
            if line.startswith("#"): # skip header
                continue

            if line.startswith("Protein"): # skip header
                continue

            entry = line.split(",")

            if float(entry[2]) < args.deepfri_filter: # low score
                continue

            query = entry[0] # sequence ID
            go = entry[1] # GO term

            annotations.setdefault(query, [set(),set()]) # default to a list of two empty lists. If the query isn't there, it wasn't present in the eggnog data
            annotations[query][1].add(go)


# ---------------------------------------------- compare annotations -------------------------------------------
# returns the jaccard index of two sets and the length of overlap in a tuple
def jaccard(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return (intersection / union, intersection)


# returns the all of the gos up to the number of jumps up on the go hierarchy
def get_gos(go_str, abs):

    go = gos[go_str]

    # base case
    if abs == 0:
        return [go.id]

    if len(go.parents) == 0: 
        return [go.id]

    parents = go.parents

    go_tree = []

    # call get_gos on each parent go with one less abstraction because going up to the parent level is an abstraction itself
    for parent in parents:
        go_tree.extend(get_gos(parent, abs - 1))

    # add the current go
    go_tree.append(go.id)

    return(go_tree)

# calls get_gos on each GO term in a collection and returns a concatenated set with all the results
def get_gos_from_collection(go_col, abs):
    go_list = []
    for go_term in go_col:
        go_list.extend(get_gos(go_term, abs))

    return set(go_list)


jaccard_scores = []

# calculate jaccard index and write to output file
with open(args.out, "w") as output:
    for query, go_sets in annotations.items():
        # skip if it doesn't have annotations from both eggnog and deepfri
        if len(go_sets[0]) == 0 or len(go_sets[1]) == 0:
            continue

        eggnog_tree = get_gos_from_collection(go_sets[0], args.abstractions)
        deepfri_tree = get_gos_from_collection(go_sets[1], args.abstractions)

        jaccard_result = jaccard(eggnog_tree, deepfri_tree)

        jaccard_scores.append(jaccard_result[0])

        output.write("\t".join((query, str(len(go_sets[0])), str(len(go_sets[1])), str(jaccard_result[1]), str(jaccard_result[0]))) + "\n")


print("Mean Jaccard index:", sum(jaccard_scores) / len(jaccard_scores))
