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

## Citations
1. Edgar, R. C. (2010). Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26(19), 2460–2461. https://doi.org/10.1093/bioinformatics/btq461
2. Emms, D. M., Lui, Y., Belcher, L., Holmes, J., & Kelly, Steven. (2025). OrthoFinder: Scalable phylogenetic orthology inference for comparative genomics.
Gligorijević, V., Renfrew, P. D., Kosciolek, T., Leman, J. K., Berenberg, D., Vatanen, T., Chandler, C., 3. Taylor, B. C., Fisk, I. M., Vlamakis, H., Xavier, R. J., Knight, R., Cho, K., & Bonneau, R. (2021). Structure-based protein function prediction using graph convolutional networks. Nature Communications, 12(1), 3168. https://doi.org/10.1038/s41467-021-23303-9
4. Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research, 47(D1), D309–D314. https://doi.org/10.1093/nar/gky1085
