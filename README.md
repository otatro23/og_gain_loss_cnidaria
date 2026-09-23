# og_gain_loss_cnidaria
In this project, we integrate gene gain and loss inference, functional annotation of uncharacterized gene sequences, and comparative single-cell RNA sequencing analysis to investigate the evolution of novel organismal traits and cell types in Cnidaria. Specifically, we focus on the evolution of the medusa life history stage in Medusozoa, parasitism in Endocnidozoa, symbiosis in Anthozoa, and the evolution and diversification of cnidocytes. 

### Gene gain and loss inference
We gather proteomes from 72 cnidarians and 25 outgroup species and assign genes into orthogroups. Following the conceptual framework of dollo parsimony, we determine the gain and loss history of the orthogroups based on our phylogenetic tree and the presence or absence of genes in each orthogroup in each extant species. 

### Functional annotation
We annotate genes gained and lost at each node with GO terms using DeepFRI, cluster annotations of gained and lost genes using GO-Figure, and group the GO clusters into broad functional categories with Microsoft Copilot. We also use EggNOG to annotate sequences with COGs. For both analyses, we determined 

### scRNA-seq analysis
We gather scRNA-seq datasets from 5 cnidarian species and integrate them into one large dataset with SAMap. We then analyze the collective expression patterns of genes gained at cnidarian focal nodes. We identify and further annotate orthogroups that were gained at many focal nodes that had enriched expression in cnidocytes, reconstructing the evolution and diversification of cnidocytes.
