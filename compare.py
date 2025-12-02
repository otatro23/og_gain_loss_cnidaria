#!/usr/bin/env python3

# command: ./compare.py go_basic.obo eggnog.tsv /mnt/home/plachetzki/omt1027/surf/97_taxa/orthofinder_run/deepfri/test

# need to deal with replacements in a similar way to the alt ids
# need to append the average jaccard index to a string

# ---------------------------------------- import -------------------------------------------
import argparse
import os
import csv

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

print(args.deepfri_filter)

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

        if len(self.alt_ids) > 0:
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

            annotations.setdefault(query, [set(),set()]) # default to a list of two empty sets. If the query isn't there, it wasn't present in the eggnog data
            annotations[query][1].add(go)


# ---------------------------------------------- compare annotations -------------------------------------------
# returns the jaccard index of two sets and the length of overlap in a tuple
def jaccard(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return (intersection / union, intersection)


# returns all of the gos up to the number of jumps up on the go hierarchy. If abs = -1, the function goes all the way up the tree
def get_gos(go_str, abs = -1):
    go = gos[go_str]

    # base case
    if abs == 0:
        return [go.id]
    if len(go.parents) == 0: 
        return [go.id]

    parents = go.parents
    go_tree = []

    # call get_gos on each parent go with one less abstraction unless -1 specified
    new_abs = abs if abs == -1 else abs - 1
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
with open("old_output", "w") as output:
    for query, go_sets in annotations.items():
        # skip if it doesn't have annotations from both eggnog and deepfri
        if len(go_sets[0]) == 0 or len(go_sets[1]) == 0:
            continue

        eggnog_tree = get_gos_from_collection(go_sets[0], args.abstractions)
        deepfri_tree = get_gos_from_collection(go_sets[1], args.abstractions)

        jaccard_result = jaccard(eggnog_tree, deepfri_tree)

        jaccard_scores.append(jaccard_result[0])

        output.write("\t".join((query, str(len(go_sets[0])), str(len(go_sets[1])), str(jaccard_result[1]), str(jaccard_result[0]))) + "\n")

# ---------------------------- compare annotations by individual deepfri GO ---------------------------

# give it 2 GOs and determine the number of abstractions from GO 1 to get GO 2 or None if they are not related
def get_abs(go_str_1, go_str_2, abs = 0):

    go_1 = gos[go_str_1]

    go_2 = gos[go_str_2]

    # return abs if the strings are ==
    if go_1.id == go_2.id:
        return abs

    # base case: no more parents
    if len(go_1.parents) == 0:
        return None

    parent_abs = []
    # call get_abs on each parent
    for parent in go_1.parents:
        this_abs = get_abs(parent, go_str_2, abs + 1)
        if not this_abs is None: # won't be able to compare None to int, so we just want to add ints to the list
            parent_abs.append(this_abs)

    # if there is nothing other than None in the list, then there were no matches. If there is an int and a none, then it will return the min of the ints. This handles the cases where one GO term has multiple parents and some might match. 
    if len(parent_abs) == 0:
        return None

    return min(parent_abs)

# calls get_abs in both directions and returns minimum
def get_min_abs(go_str_1, go_str_2):
    results = []

    forwards = get_abs(go_str_1, go_str_2)
    backwards = get_abs(go_str_2, go_str_1)

    if not forwards is None:
        results.append(forwards)

    if not backwards is None:
        results.append(backwards)

    if len(results) == 0:
        return None
    else:
        return min(results)

# testing get_min_abs()
#print("Same GO:", get_min_abs("GO:0000001", "GO:0000001"))
#print("1 step:", get_min_abs("GO:0000001", "GO:0048308"))
#print("1 step backwards:", get_min_abs("GO:0048308", "GO:0000001"))
#print("2 steps:", get_min_abs("GO:0000001", "GO:0006996"))
#print("Not in same ontology:", get_min_abs("GO:0000001", "GO:0000006"))

# goal: create a tsv that has query, deepfri GO term, # abs to any eggnog GO (or None), highest jaccard index
# store all jaccards for summary stats
jaccard_scores = []
# store all abs_scores for summary stats
abs_scores = []

# calculate jaccard index and write to output file
with open(args.out, "w") as output:
    output.write("# DeepFRI annotations filtered by Score >= {0}\n".format(args.deepfri_filter))
    writer = csv.writer(output, delimiter = "\t")
    for query, go_sets in annotations.items():
        # skip if it doesn't have annotations from both eggnog and deepfri
        if len(go_sets[0]) == 0 or len(go_sets[1]) == 0:
            continue

        eggnog_trees = []
        for eggnog_go in go_sets[0]:
            eggnog_trees.append(get_gos(eggnog_go))

        #print(query, go_sets)

        # iterate through all deepfri GOs 
        for deepfri_go in go_sets[1]:
            abstractions = []

            # find the jaccard between the deepfri go and each eggnog go and store the max in max_jaccard
            jaccards = []
            deepfri_tree = get_gos(deepfri_go)
            for eggnog_tree in eggnog_trees:
                jaccards.append(jaccard(set(eggnog_tree), set(deepfri_tree))[0])
            max_jaccard = max(jaccards)
            #print(max_jaccard)

            jaccard_scores.append(max_jaccard)

            # find the lowest number of abstractions between the deepfri go and all eggnog gos
            for eggnog_go in go_sets[0]:
                abs_result = get_min_abs(deepfri_go, eggnog_go)
                if not abs_result is None:
                    abstractions.append(abs_result)

            min_abstractions = -1 if len(abstractions) == 0 else min(abstractions)
            #print(min_abstractions)

            if not min_abstractions is None:
                abs_scores.append(min_abstractions)

            writer.writerow([query, deepfri_go, min_abstractions, max_jaccard])

summary_exists = os.path.exists("compare_summary.txt")
print(summary_exists)

with open("compare_summary.txt", "a") as summary:
    writer = csv.writer(summary, delimiter = "\t")

    if not summary_exists: # write header if the file doesnt already exist
        writer.writerow(["mean_num_abs", "mean_jaccard", "filter"])

    mean_abs = sum(abs_scores)/len(abs_scores)
    mean_jaccard = sum(jaccard_scores)/len(jaccard_scores)
    writer.writerow([mean_abs, mean_jaccard, args.deepfri_filter])
