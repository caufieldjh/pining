# pining
Protein Interaction Networks and Integration with Novel Graphs

A system for retrieving protein-protein interaction data for a set of proteins, integrating the resulting network with properties observed for the proteins (e.g., experimental observations), and identifying subgraphs within the interaction network most relevant to differences in the observed properties.

(Work in progress.)

This system includes several different modules.

## pining_for_new_data.py

This part of the system is intended for retrieving protein-protein interactions from databases of published interactions, given a set of protein IDs. Produces combined, annotated versions of the network of the retrieved interactions.

IntAct (http://www.ebi.ac.uk/intact/) and BioGRID (https://thebiogrid.org/) are searched as local databases.

Downloads OG mapping files from EggNOG 4.5 (http://eggnogdb.embl.de). Uses base NOGs and eukaryote NOGs (euNOGs or KOGs) by default for most general coverage, but can retrieve others.

### Requirements 
Requires [pandas](http://pandas.pydata.org/).
Also requires [Neo4j](https://neo4j.com/) - please install it and [py2neo](http://py2neo.org/v3/). Note that Neo4j requires the Java 8 runtime  - [see its documentation for more details](https://neo4j.com/docs/operations-manual/current/installation/linux/debian/). Setting up and accessing the Neo4j database will require administrative or sudo privileges on your system. 

The combined databases and conversion files require ~10 Gb of disk space and are used in other modules.

### Usage
Run as:
`pining_for_new_data.py -inputfile INPUT_FILE_NAME`

#### Input
A text file containing a list of UniprotAC IDs, one per line,
to be used as the interactors to search for.
Will accept other data if they are present but the file must be in long format (that is, one observation per line).

Uniprot IDs will be assumed to be in the first column.

For additional data fields, pining_for_new_data.py will determine whether the data category is qualitative or quantitative. Qualitative variables will be treated as categories and subcategories, such that an input that looks like the following (tab delimited):

`uniprot	group1	group2	value`

`M3W421	bags	pos	0.6133`

`M3W421	bags	neg	0.0182`

`M3WX82	yarn	neg	0.1232`

`...`

will be interpreted as "for M3W421, bags and pos is 0.6133" and so on.

#### Output
* A list of interactions, in PSI-MI TAB format, involving only interactors in the target set. See the following link for full details of PSI-MI TAB format:

  http://wiki.reactome.org/index.php/PSI-MITAB_interactions

* Also provides the same network, but as a simple edgelist.

* An annotation file containing all interactors from the input set with interactions found from the database search. These annotations include eggNOG OG assignments where available.

## pining_for_turnover.py

Will provide methods for comparing proteins and their interactions in terms of protein turnover dynamics, given proteomics data sets expressing protein dynamics in terms of median values across multiple experimental conditions and groups (e.g., in vivo proteomes from multiple animal strains).

## pining_for_hearts.py

Will provide methods for comparing proteins and their interactions in terms relevant to cardiac function and cardiovascular disease.
