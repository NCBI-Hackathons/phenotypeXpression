[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1210203.svg)](https://doi.org/10.5281/zenodo.1210203)

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/PhenoX.png" width="300" align= "middle"/>

# phenotypeXpression
## Synopsis

### Subclassification of disease states based on the intersection of literature and expression

PhenotypeXpression (PhenoX) provides a lightweight tool for clinicians to rapidly aggregate the publicly available gene expression data and literature available from the NCBI. It requires a minimum of resources and can quickly survey publicly available research. By combining the subtype HPO and DOID information with GEO dataset expression profiles, each sub-classification of a disease is given a characteristic fingerprint. This enables clinicians to (1) identify a set of genes that have altered expression in a phenotype and its related targets to aid in research, (2) use those target genes to inform treatment selection by affected pathway, and (3) rapidly identify public datasets that are available for comparison with their own gene expression data.

Flow Chart here

## Quick Start

```bash
Usage here
```

### Minimum Inputs

The primary inputs are a MeSH term and your email. The MeSH term should specify the parent condition for which you would like to derive subgroup information. The email is required to batch query NCBI databases.

### Optional Arguments

Additional arguements include a specification of and output file prefix (-o, --out), and two thresholds for clustering: (1) minimum count of gene X across all GDS sets and (2) minimum number of genes for a given GDS set. There two thresholds determine the sensitivity of the clustering. FURTHER EXPLANATION

## Motivation

Before personalized medicine can be fully implemented, disease classification ontologies will need to gain many orders of magnitude. However, this then complicates the role of clinicians seeking to use a more tailored diagnosis. Differential expression within a disease set, its subtypes, can have a significant impact on treatment effect.

USE CASES

## Installation

Implemented in python3 requires >=3.4. To install in git clone and unpack or... to run from command line, add to path:

```
git clone https://github.com/NCBI-Hackathons/phenotypeXpression.git
cd phenotypeXpression
./setup.sh
``` 

Examples of each input file type are provided in the test subfolder. For instance:

```
source activate phenoX
sample run here
```

Should produce the following files:

```
ls file result
```


### Dependencies

```
biopython
spacy
spacy-lookup
json
pickle
collections
```

Note on deps. CLUSTERING DEPS

### Memory/System requirements

Note on system requirements

## Details on Results...

For the given MeSH term, PhenoX will attempt to generate all subclassifications of the condition based on the available expression data in the NCBI Gene Expression Omnibus. Each subclass will have two major features, a word cloud of phenotypic terms associated with the seach tem via a literature search, and a gene expression profile of differentially expressed genes. The output will include several file types described in detail below to 

### Word Clouds

<a title="By Monikasj [CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0)], from Wikimedia Commons" href="https://commons.wikimedia.org/wiki/File:Wikipedia-word-cloud.jpg"><img width="512" alt="Wikipedia-word-cloud" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1e/Wikipedia-word-cloud.jpg/512px-Wikipedia-word-cloud.jpg"></a>

Word clouds are generated for each subcluster based on a combination of Human Phenotype Ontology and Human Disease Ontology search terms. Additionally, MeSH child terms of the parent query are counted to see if they correlate strongly with any given cluster. 

### Expression Profiles

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/heat.png" width="1000" align="middle"/>

Something about heatmap/dendrogram. PERHAPS A HEATMAP ALIGNED TO DENDROGRAM OF COUNTS (OR BOOTSTRAPS?). BE SURE TO INCLUDE DATA ON MULTISCALE BOOTSTRAP CLUSTERING. 

### Cluster map

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/Figure4ab.jpg" width="1000" align="middle"/>

An overall cluster map showing pairwise distances between GDS datasets based on differential gene expression. CLUSTERED BY ABOVE, COLORED BY CLUSTER, 

## Authors

Agnes Bao  
Robert R Butler III  
Huaiying Lin  
Subhajit Sengupta  
Lucy Lu Wang  

## References

<sup>1</sup>    
<sup>2</sup>    
<sup>3</sup>    
<sup>4</sup>    
<sup>5</sup>    
