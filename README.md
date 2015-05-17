# grinn
Graph database and R package for omics integration

Description
=========
grinn is a bioinformatics platform contains an internal graph database (Neo4j), and the R package for omics studies.
The package allows reconstruction of different network types e.g. metabolite-protein-gene, metabolite-protein, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
grinn incorporates data from several databases including KEGG, SMPDB, HMDB, REACTOME and CheBI.
It uses the correlation-based network analysis to estimate relationships among different omics levels independently from domain knowledge, and with the internal graph database it allows fast integration of domain knowledge to annotate unknown metabolites.
