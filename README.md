# Grinn
a Graph database and R package for omic data integration

Version: 2.8 (13 October 2016)

Description
=========
Grinn is a bioinformatics platform contains an internal graph database (Neo4j), and the R package for -omic studies.
Grinn databases incorporate data from several databases including KEGG, SMPDB, HMDB, REACTOME, CheBI, UniProt and ENSEMBL.
The R package allows reconstruction of different network types e.g. metabolite-protein-gene, metabolite-protein, metabolite-pathway, protein-gene, protein-pathway and gene-pathway.
Grinn applies different correlation-based network analyses to estimate relationships among different omics levels independently from domain knowledge, and with the internal graph database it provides rapid integration of domain knowledge i.e. to aid annotation of unknown metabolites.

Installation
=========
  1. Require [R software](https://www.r-project.org/)
  2. Require [shiny](http://shiny.rstudio.com/)
  3. Install grinn R package using the following commands
```
#Install devtools R package, if not exist
install.packages("devtools")

#Install grinn package
library(devtools)
devtools::install_github("kwanjeeraw/grinn")
library(grinn)
```
Graph databases
=========
The internal graph database is a part of the Grinn software to compute the networks. Graph databases are available for Human, Arabidopsis, Mouse, Rat, Saccharomyces cerevisiae and Escherichia coli k-12. The human database is provided by default and can be accessed directly after package installation. 

Alternatively the graph databases can be installed <b>locally</b>. The graph databases are available on request. 

<b>Local database installation</b>
  1. Require Neo4j-community 2.2.0 to 2.3.0 for the internal graph database

    - Download and then unzip [Neo4j server](http://neo4j.com/download/)

    - Extract and move the graph database files to the Neo4j server directory

    - Start the Neo4j server, 
    
    for windows: Double-click on %NEO4J_HOME%\bin\Neo4j.bat 
    
    for linux: ./bin/neo4j start 
    
    for more details see [here](http://neo4j.com/docs/stable/server-installation.html)  
  2. Switch between databases
```
#Change the internal database by providing the database url, e.g. "http://database.location:7474/db/data/"
setGrinnDb("http://localhost:7474/db/data/")

#Check current internal database location
getGrinnDb()
```

Documentation
=========
see [homepage](http://kwanjeeraw.github.io/grinn/)

Updates
=========
#### version 2.8 (13/10/16)
* Include partial correlation analysis

see all [features](NEWS.md)

References
=========
- Kanehisa M, Goto S, Sato Y, Furumichi M, Tanabe M: KEGG for integration and interpretation of large-scale molecular data sets. Nucleic acids research 2012, 40(Database issue):D109-114.
- Jewison T, Su Y, Disfany FM, Liang Y, Knox C, Maciejewski A, Poelzer J, Huynh J, Zhou Y, Arndt D et al: SMPDB 2.0: big improvements to the Small Molecule Pathway Database. Nucleic acids research 2014, 42(Database issue):D478-484.
- Croft D, Mundo AF, Haw R, Milacic M, Weiser J, Wu G, Caudy M, Garapati P, Gillespie M, Kamdar MR et al: The Reactome pathway knowledgebase. Nucleic acids research 2014, 42(Database issue):D472-477.
- Hastings J, de Matos P, Dekker A, Ennis M, Harsha B, Kale N, Muthukrishnan V, Owen G, Turner S, Williams M et al: The ChEBI reference database and ontology for biologically relevant chemistry: enhancements for 2013. Nucleic acids research 2013, 41(Database issue):D456-463.
- Wishart DS, Jewison T, Guo AC, Wilson M, Knox C, Liu Y, Djoumbou Y, Mandal R, Aziat F, Dong E et al: HMDB 3.0--The Human Metabolome Database in 2013. Nucleic acids research 2013, 41(Database issue):D801-807.
- UniProt C: UniProt: a hub for protein information. Nucleic acids research 2015, 43(Database issue):D204-212.
- Cunningham F, Amode MR, Barrell D, Beal K, Billis K, Brent S, Carvalho-Silva D, Clapham P, Coates G, Fitzgerald S et al: Ensembl 2015. Nucleic acids research 2015, 43(Database issue):D662-669.

License
=========
[GNU General Public License (v3)](https://github.com/kwanjeeraw/grinn/blob/master/LICENSE)

Citation
=========
- Wanichthanarak K, Fahrmann JF, Grapov D: Genomic, Proteomic, and Metabolomic Data Integration Strategies. Biomark Insights 2015, 10:1-6, doi:10.4137/BMI.S29511.
