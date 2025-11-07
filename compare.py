#!/usr/bin/env python3

# ------------------- import --------------------
import argparse

# ---------------- arguments ---------------------
# user must put in the go.obo file and the output of both eggnog and deepfri
parser = argparse.ArgumentParser()
parser.add_argument("go", help = "path to go.obo file")
parser.add_argument("eggnog", help = "path to eggnog output")
parser.add_argument("deepfri", help = "path to deepfri output file")
args = parser.parse_args()

# ------------------ parse go.obo ----------------

# class that contains go terms and instance attributes for their "is_a" relationships, whether they are obsolete or not, and their replacement if they have one
class Go:
    def __init__(self, id):
        self.id = id
        self.is_obsolete = False
        self.replacement = None
        self.is_a = []

    def add_is_a(self, go):
        self.is_a.append(go)

    def __str__(self):
        string = self.id
        if self.is_obsolete:
            string += " is obsolete."
            if self.replacement is not None:
                string += " Replaced by: " + self.replacement
        elif self.is_a != []:
            string += " is: " + str(self.is_a)
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
            gos[go_term].add_is_a(is_a)

        # changes the Go instance attribute of is_obsolete to True if it is
        if line.startswith("is_obsolete"):
            gos[go_term].is_obsolete = True

        # adds in replacement go if it is obsolete and it has a replacement
        if line.startswith("replaced_by"):
            rep = line.rstrip().split(" ")[1]
            gos[go_term].replacement = rep

for go in gos.values():
    print(go)

print("number of gos:", len(gos))
