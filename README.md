# grinn
a Graph database and R package for omic data integration

Version: 2.5 (12 January 2016)

Description
=========
Grinn is a bioinformatics platform contains an internal graph database (Neo4j), and the R package for -omic studies.
Grinn databases incorporate data from several databases including KEGG, SMPDB, HMDB, REACTOME, CheBI, UniProt and ENSEMBL.
The R package allows reconstruction of different network types e.g. metabolite-protein-gene, metabolite-protein, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
Grinn applies different correlation-based network analyses to estimate relationships among different omics levels independently from domain knowledge, and with the internal graph database it provides rapid integration of domain knowledge i.e. to aid annotation of unknown metabolites.

Installation
=========
  1. Require [R software](https://www.r-project.org/)
  2. Install grinn R package using the following commands
```
#Install devtools R package, if not exist
install.packages("devtools")

#Install grinn package
library(devtools)
devtools::install_github("kwanjeeraw/grinn")
library(grinn)
```
Grinn databases
=========
Grinn internal database is a part of the Grinn software to compute the networks. Grinn databases are available for Human, Arabidopsis, Mouse, Saccharomyces cerevisiae and Escherichia coli k-12. The human database is provided by default and can be accessed directly after package installation. 

Alternatively the Grinn databases can be installed <b>locally</b>. The Grinn databases are available at https://sourceforge.net/projects/grinn/files/grinnDatabases/. 

<b>Local database installation</b>
  1. Require Neo4j-community >= 2.2.0 for the Grinn internal database

    - Download and then unzip [Neo4j server](http://neo4j.com/download/)

    - Extract and move the grinn database files to the Neo4j server directory

    - Start the Neo4j server, 
    
    for windows: Double-click on %NEO4J_HOME%\bin\Neo4j.bat 
    
    for linux: ./bin/neo4j start 
    
    for more details see [here](http://neo4j.com/docs/stable/server-installation.html)  
  2. Switch between databases
```
#Change Grinn internal database by providing the database url, default location is "http://grinn.genomecenter.ucdavis.edu:7474/db/data/"
setGrinnDb("http://localhost:7474/db/data/")

#Check current Grinn internal database location
getGrinnDb()
```

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinn/)

Updates
=========
#### version 2.5 (12/01/16)
* Change parameter names in all correlation-related functions
* Update manual
* Bug fixed

see all [features](NEWS.md)

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinn/blob/master/LICENSE)
