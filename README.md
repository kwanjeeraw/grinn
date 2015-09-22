# grinn
Graph database and R package for omics integration

Version: 2.3 (22 September 2015)

Description
=========
grinn is a bioinformatics platform contains an internal graph database (Neo4j), and the R package for omics studies.
grinn incorporates data from several databases including KEGG, SMPDB, HMDB, REACTOME, CheBI, UniProt and ENSEMBL.
The package allows reconstruction of different network types e.g. metabolite-protein-gene, metabolite-protein, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
grinn applies different correlation-based network analyses to estimate relationships among different omics levels independently from domain knowledge, and with the internal graph database it provides rapid integration of domain knowledge i.e. to aid annotation of unknown metabolites.

Installation
=========
  1. Require [R software](https://www.r-project.org/)
  2. Install grinn R package and dependent R packages using the following commands
```
#Install dependent R packages, if not exist
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute"))
install.packages(c("RCurl", "jsonlite", "igraph", "WGCNA", "stringi")) 

#Install devtools R package, if not exist
install.packages("devtools")

#Install grinn package
devtools::install_github("kwanjeeraw/grinn")
library(grinn)
```
Grinn databases
=========
Grinn human database is provided by default and can be accessed directly after package installation. 
<b>Alternatively</b> grinn databases are available for <b>local use</b>: Human database, Arabidopsis database, Mouse database Saccharomyces cerevisiae database and Escherichia coli k-12 database.
<b>Local database installation</b>
  1. Require Neo4j-community >= 2.2.0 for the grinn internal database, please send us an email for the database files

    - Download and then unzip [Neo4j server](http://neo4j.com/download/)

    - Extract and move the grinn database files to the Neo4j server directory

    - Start the Neo4j server, for windows: Double-click on %NEO4J_HOME%\bin\Neo4j.bat, for linux: ./bin/neo4j start 
for more details see [here](http://neo4j.com/docs/stable/server-installation.html)  
  2. Switch between databases
```
#Change grinn internal database by providing the database url, default location is "http://grinn.genomecenter.ucdavis.edu:7474/db/data/cypher"
setGrinnDb("http://localhost:7474/db/data/cypher")

#Check current grinn internal database
getGrinnDb()
```

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinn/)

Updates
=========
#### version 2.3 (22/09/15)
* Saccharomyces cerevisiae database V.1 incorporating data from KEGG, SMPDB, REACTOME, CheBI, UniProt and ENSEMBL.
* Escherichia coli k-12 database V.1 incorporating data from KEGG, SMPDB, CheBI and UniProt.
* Bug fixed, return number of edges
* Include functions to get and set grinn database location

see all [features](NEWS.md)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinn/blob/master/LICENSE)
