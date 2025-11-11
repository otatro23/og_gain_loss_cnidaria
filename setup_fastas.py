#!/usr/bin/env python3

# goal: ensure that the genes in the fasta file that we will blast match up with the gene IDs in the single cell counts data

# -------------------------------------- acropora ----------------------------------------
# parse through the fasta file and change headers to the feature ids that match the features.tsv

old_fasta = "/mnt/home/plachetzki/omt1027/surf/97_taxa/single_cell/genera/acropora/data/acropora.fa"
map = "/mnt/home/plachetzki/omt1027/surf/97_taxa/single_cell/geps_cnid/acropora/SAMap/feature_map.tsv"
new_fasta = "/mnt/home/plachetzki/omt1027/surf/97_taxa/single_cell/geps_cnid/acropora/SAMap/acropora.fasta"

# parse through mapping file to create a dictionary with key: current ID and value: features.tsv ID
genes = {}
try:
    with open(map, "r") as mapping:
        for line in mapping:
            cols = line.rstrip().split(" ")
            genes[cols[0]] = cols[1]
except IOError:
    print("Could not open {0}".format(map))

# parse through old fasta and write a new fasta with the header IDs switched
try:
    with open(old_fasta, "r") as old:

        try:
            with open(new_fasta, "w") as new:
                for line in old:
                    if line.startswith(">"):
                        nr_id = line.split(" ")[0][1::]
                        new.write(">" + genes.get(nr_id, nr_id) + "\n")
                    else:
                        new.write(line)
        except IOError:
            print("Could not open {0}".format(new_fasta))

except IOError:
    print("Could not open {0}".format(old))
