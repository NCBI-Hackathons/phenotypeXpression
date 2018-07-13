[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1311476.svg)](https://doi.org/10.5281/zenodo.1311476)

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/PhenoX.png" width="300" align= "middle"/>

# phenotypeXpression

### Subclassification of disease states based on the intersection of literature and expression

PhenotypeXpression (PhenoX) is a prototype precision medicine tool for clinicians to rapidly aggregate the publicly available gene expression data and literature available from the NCBI. It requires a minimum of resources and can quickly survey publicly available research. By combining the subtype HPO and DOID information with GEO dataset expression profiles, each sub-classification of a disease is given a characteristic fingerprint. This enables clinicians to (1) identify a set of genes that have altered expression in a phenotype and its related targets to aid in research, (2) use those target genes to inform treatment selection by affected pathway, and (3) rapidly identify public datasets that are available for comparison with their own gene expression data.

Flow Chart here

### Minimum Inputs

The primary inputs are a MeSH term and your email. The MeSH term should specify the parent condition for which you would like to derive subgroup information. The email is required to batch query NCBI databases.

### Optional Arguments

Additional arguements include a specification of and output file prefix (-o, --out), and two thresholds for clustering: (1) minimum count of gene X across all GDS sets and (2) minimum number of genes for a given GDS set. There two thresholds determine the sensitivity of the clustering. FURTHER EXPLANATION

## Motivation

Before personalized medicine can be fully implemented, disease classification ontologies will need to gain many orders of magnitude. However, this then complicates the role of clinicians seeking to use a more tailored diagnosis. Differential expression within a disease set, its subtypes, can have a significant impact on treatment effect. In many cases, the subtypes for complex diseases have yet to be defined or named. These subtypes can be defined by gene expression profiles or phenotypic traits. This project aims to produce a quantitative and qualitative description of subtypes of common diseases, using gene expression signatures from the NCBI Gene Expression Omnibus (GEO), as well as phenotypic traits mined from the literature.

## Installation

Implemented in python3 requires >=3.4 and R >=3.5. To install:

```
git clone https://github.com/NCBI-Hackathons/phenotypeXpression.git
cd phenotypeXpression
sh setup.sh
``` 

Examples of each input file type are provided in the test subfolder. For instance:

```
source activate phenoX
python run_phenox.py A.N.Other@example.com "Psoriasis"
```

Should produce the following files:

* Dendgrogram of GDS with bottstraps and cluster boxes
* A word cloud for each cluster box
* Heatmap of differentially expressed genes per GDS (columns)
* A pairwise distance graph of GDS by differing gene set

### Dependencies

```
biopython
spacy
spacy-lookup
json
pickle
collections
tqdm
```

Note on deps. CLUSTERING DEPS

## Details on Results...

For the given MeSH term, PhenoX will attempt to generate all subclassifications of the condition based on the available expression data in the NCBI Gene Expression Omnibus. Each subclass will have two major features, a word cloud of phenotypic terms associated with the seach tem via a literature search, and a gene expression profile of differentially expressed genes. The output will include several file types described in detail below to 

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/Output.png" width="600" align="middle"/>

### Word Clouds

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/wordcloud.png" width="1000" align="middle"/>

Phenotype and disease ontological terms are enriched using occurence in PubMed abstract, title, and keywords. The spacy module scrapes terms which have been serialized by pickle for speed. Term frequency data is gathered for each GEO DataSet, and the counts are merged for the gene expression clusters. Word_cloud is then used to visualize the result, with a separate png image files for each cluster.

### Expression Profiles

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/hcluster-1.png" width="600" align="middle"/>
<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/heatmap-1.png" width="600" align="middle"/>

GEO profiles which are filtered for differential expression are matched to the query term by the Entrez eSearch method. Resulting GEO DataSets, each with differential NCBI GeneID sets and an associated PubMed ID are then used for downstream analysis. The GDS are then clustered using pvclust in R, which uses hclust for hierarchical clustering based on the binary gene presence or absence across all the GDS. GDS with only a single differential gene are excluded to reduce noise. Clustering generates P-values via approximately unbiased (AU) method and bootstrap probabilities (BP) at each node. Branches are defined as sub-classification clusters with shared AU node values greater than 95%. The GDS set in each sub-classification cluster is then used for PubMedID associated literature data.

### Cluster map

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/image.png" width="600" align="middle"/>

An overall cluster map showing pairwise distances between GDS datasets based on absoulte difference in number of genes per GDS. Based on the gene sets defined above.

## Authors

Agnes Bao  
Robert R Butler III  
Huaiying Lin  
Subhajit Sengupta  
Lucy Lu Wang  
