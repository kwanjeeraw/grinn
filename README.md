# grinn
Graph database and R package for omics integration

Version: 2.2 (18 June 2015)

Description
=========
grinn is a bioinformatics platform contains an internal graph database (Neo4j), and the R package for omics studies.
grinn incorporates data from several databases including KEGG, SMPDB, HMDB, REACTOME, CheBI, UniProt and ENSEMBL.
The package allows reconstruction of different network types e.g. metabolite-protein-gene, metabolite-protein, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
grinn applies different correlation-based network analyses to estimate relationships among different omics levels independently from domain knowledge, and with the internal graph database it provides rapid integration of domain knowledge i.e. to aid annotation of unknown metabolites.

Installation
=========
* To use grinn internal database (server version, access human database only):
  1. Require R software and dependent R packages including RCurl, WGCNA, stringi, igraph, jsonlite, preprocessCore, impute, GO.db.
  2. Download grinn R package (server version): for windows click [here](http://kwanjeeraw.github.io/grinn/extra/server/grinn_2.2.zip), for linux click [here](http://kwanjeeraw.github.io/grinn/extra/server/grinn_2.2.tgz), and install using the following commands.
```
#Install dependent R packages, if not exist
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute"))
install.packages(c("RCurl", "jsonlite", "igraph", "WGCNA", "stringi")) 

#Install grinn package server version
install.packages("[path to grinn_2.2]", repos = NULL)
library(grinn)
```
* <b>Alternatively,</b> grinn internal databases are available for local use: Arabidopsis database, Mouse database and Human database. Please send us an email for the database files.
  1. Require Neo4j-community >= 2.2.0 for the grinn internal database (local version), please send us an email for the grinn database files, currently available: Arabidopsis database, Mouse database and Human database.

    - Download and then unzip [Neo4j server](http://neo4j.com/download/)

    - Extract and move the grinn database files to the Neo4j server directory

    - Start the Neo4j server, for windows: Double-click on %NEO4J_HOME%\bin\Neo4j.bat, for linux: ./bin/neo4j start 
for more details see [here](http://neo4j.com/docs/stable/server-installation.html)  
  2. Install grinn R package using the following commands.
```
#Install devtools R package, if not exist
install.packages("devtools")

#Install grinn package
devtools::install_github("kwanjeeraw/grinn")
library(grinn)
```

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinn/)

Updates
=========
#### version 2.2 (18/06/15)
* Arabidopsis database V.1 incorporating data from KEGG, SMPDB, REACTOME, CheBI, UniProt and ENSEMBL.
* Add node attributes to correlation networks

see all [features](NEWS.md)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinn/blob/master/LICENSE)
