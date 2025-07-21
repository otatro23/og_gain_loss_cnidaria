# og_gain_loss_cnidaria
This project aims to categorize gains and losses of function throughout the cnidarian lineage based on gene orthogroups.

Data: Proteomes of 72 Cnidarians, Mnemiopsis, Amphimedon, Homo, and Drosophila 

## Workflow
1. Orthofinder was run on the .fasta files to assign genes from the proteomes into orthogroups
2. Usearch was run to find centroids: groups of sequences to summarize the diversity of each orthogroup.
3. DeepFRI was run on each centroid sequence to annotate it with GO terms.
4. EggNOG was also run on each centroid to annotate it with GO terms, COGs, and KEGG Pathways.
5. The concept of Dollo Parsimony, that genes do not evolve twice due to that being statistically unlikely, was used to infer orthogroup gains and losses at each node in Cnidarian phylogeny.
6. Data from EggNOG, DeepFRI, and orthogroup gains and losses were combined to yield a better understanding of gain and loss of function at various nodes in Cnidaria. 
