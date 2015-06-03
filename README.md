# grinn
Graph database and R package for omics integration

Version: 2.1 (03 June 2015)

Description
=========
grinn is a bioinformatics platform contains an internal graph database (Neo4j), and the R package for omics studies.
grinn incorporates data from several databases including KEGG, SMPDB, HMDB, REACTOME, CheBI, UniProt and ENSEMBL.
The package allows reconstruction of different network types e.g. metabolite-protein-gene, metabolite-protein, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
grinn applies different correlation-based network analyses to estimate relationships among different omics levels independently from domain knowledge, and with the internal graph database it provides rapid integration of domain knowledge i.e. to aid annotation of unknown metabolites.

Installation
=========
* Install grinn R package for grinn database query, correlation analysis and network construction.

```
install.packages("devtools")
devtools::install_github("kwanjeeraw/grinn")
library(grinn)
```

* Require [Neo4j-community](http://neo4j.com/download/) >= 2.2.0 for the grinn internal database (local version)
please send us an email for the grinn database files.
** Download and then unzip Neo4j server.
** Extract and move the grinn database files to the Neo4j server directory
** Start the Neo4j server, for windows: Double-click on %NEO4J_HOME%\bin\Neo4j.bat, for linux: ./bin/neo4j start 
for more details see [here](http://neo4j.com/docs/stable/server-installation.html)

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinn/)

Updates
=========
see [features](NEWS.md)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinn/blob/master/LICENSE)