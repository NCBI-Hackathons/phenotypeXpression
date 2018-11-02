[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1311475.svg)](https://doi.org/10.5281/zenodo.1311475)

<div style="text-align:center"><img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/PhenoX.png" width="300"/></div>

# phenotypeXpression

### Subclassification of disease states using public gene expression data and literature

PhenotypeXpression (PhenoX) is a proof of principle precision medicine tool for clinicians to quickly discover potential subtypes of known diseases. PhenoX rapidly aggregates gene expression data from multiple public studies, and then mines PubMed literature to develop novel disease subclassifications and expression profiles. By combining literature occurrences of subtype HPO and DOID information with GEO dataset expression profiles, each sub-classification of a disease is given a characteristic terminology fingerprint. Here we demonstrate the use of this tool for a common disease—Psoriasis—and identify two clusters with unique gene expression profiles.

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/2018-11-02_Phenox_Fig1.png" width="1000" align="middle"/>

### Minimum Inputs

The primary inputs are a MeSH term and your email. The MeSH term should specify the parent condition for which you would like to derive subgroup information. The email is required to batch query NCBI databases.

### Optional Arguments

Additional arguements include a specification of and output file prefix (-o).

## Motivation

Before personalized medicine can be fully implemented, disease classification ontologies will need to gain many orders of magnitude. However, this then complicates the role of clinicians seeking to use a more tailored diagnosis. Differential expression within a disease set, its subtypes, can have a significant impact on treatment effect. In many cases, the subtypes for complex diseases have yet to be defined or named. These subtypes can be defined by gene expression profiles or phenotypic traits. This project aims to produce a quantitative and qualitative description of subtypes of common diseases, using gene expression signatures from the NCBI Gene Expression Omnibus (GEO), as well as phenotypic traits mined from the literature.

## Installation

Implemented in Python 3.6 and R. For simplicity, the install script setup.sh utilizes a conda environment for a self contained environment. Tested environments include Anaconda3 4.3 and Anaconda3 4.4. To install:

```
git clone https://github.com/NCBI-Hackathons/phenotypeXpression.git
cd phenotypeXpression
sh setup.sh
``` 

After installation, the program can be run from the same folder:

```
source activate phenoX
python run_phenox.py -e A.N.Other@example.com "Psoriasis"
```

Should produce the following files in the phenotypeXpression/output directory:

* Dendgrogram of GDS with node P-values and cluster boxes.
* A word cloud for each cluster box.
* Heatmap of differentially expressed genes per GDS (columns).
* A pairwise distance graph of each GDS by absolute difference in gene sets.
* A newick tree file
* A cluster batch statistics text file.
* A cluster HPO/DOID term frequency text file.

### Dependencies

python dependencies
```
biopython==1.72
spacy==2.0.12
spacy-lookup==0.0.2
tqdm==4.23.4
numpy==1.14.5
pandas==0.23.3
wordcloud==1.4.1
scikit-learn==0.19.2
pydendroheatmap==1.5
scipy==1.1.0
```

conda dependencies
```
rpy2
r-ape
r-pvclust
r-circlize
matplotlib
readline
r-gplots
```

## Details on Results

For the given MeSH term, PhenoX will attempt to generate all subclassifications of the condition based on the available expression data in the NCBI Gene Expression Omnibus. Each subclass will have two major features, a word cloud of phenotypic terms associated with the seach term via a literature search, and a gene expression profile of differentially expressed genes. The output will include several file types described in detail below: 

### Word Clouds

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/Psoriasis_GDS_wordcloud.png" width="1000" align="middle"/>

The PubMed identifier cited by each GDS is retained. Differential genes conserved across entire clusters are expanded using Hugo Gene Nomenclature (HGNC) gene names and aliases (Povey et al, 2001). This expanded gene list is used with the original MeSH term to query the PubMed database for additional PubMed articles potentially related to each cluster. For each cluster, we consolidate the titles, abstracts, and keywords corresponding to associated PubMed articles, and identify HPO and DOID terms within the article texts. The Spacy Lookup extension library is used for named entity recognition (NER) based on HPO and DOID dictionary terms. Term frequency data is gathered for each article, and the counts are merged for each cluster. Word clouds are generated and used to visualize the term frequency differences between clusters, with a separate term frequency list exported as a text file.

### Expression Profiles

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/phenox_Psoriasis_hierarchical_clusters.png" width="600" align="middle"/>
<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/phenox_Psoriasis_heatmap.png" width="600" align="middle"/>

GEO profiles matching the MeSH term are filtered using the “up down genes” filter. Resulting matches are parsed to return a list of differentially expressed gene names and corresponding GEO DataSets (GDS) (Barrett et al, 2013). All GDS are retrieved and hierarchically clustered using R-pvclust hclust (Suzuki et al, 2006); the presence or absence of each gene in a GDS is represented as a binary variable and used for clustering. GDS with only a single differential gene are excluded to reduce noise. Clustering is based on euclidean distance in order to generate p-values via an approximately unbiased method with bootstrap probabilities at each node. Clusters are defined when the derived p-value is above the threshold of alpha greater than 0.95. After cluster assignment, the GDS in each sub-classification cluster is used for PubMed literature retrieval. To check for batch effects during clustering, the number of samples, date of submission, and platform type for each GDS is collected. Each cluster is compared to the overall distribution for a significant difference (ɑ < 0.01) using the Kolmogorov-Smirnov test and chi-squared test for numerical and categorical data, respectively. Batch statistics are output to text file.

### Pairwise Distance Graph

<img src="https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/docs/phenox_Psoriasis_dist_graph.png" width="600" align="middle"/>

A pairwise distance graph showing the GDS datasets based on absoulte difference in number of genes per GDS is also produced. Based on the gene sets defined above.

## Authors

Lucy Lu Wang  
Huaiying Lin  
Agnes Bao  
Subhajit Sengupta  
Robert R Butler III  
