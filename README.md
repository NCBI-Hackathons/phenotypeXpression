[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1210203.svg)](https://doi.org/10.5281/zenodo.1210203)

![alt text](https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/PhenoX.png "PhenoX Logo")

# phenotypeXpression
## Synopsis

### Subclassification of disease states based on the intersection of literature and expression

Words, words, words

## Quick Start

```bash

```

### Minimum Inputs



### Optional Arguments

## Motivation

Abstract

## Installation

Implemented in python3 requires >=3.4. To install in git clone and unpack or... to run from command line, add to path:

```
export PATH=$PATH:path/to/folder/phenotypeXpression/phenox
``` 
or add to your ~/.bash_profile

Examples of each input file type are provided in the test subfolder. For instance:

```
sample run here
```

Should produce the following files:

```
result
```


### Dependencies

The following from my pipenv:

```
biopython==1.70
  - numpy [required: Any, installed: 1.14.0]
pandas==0.22.0
  - numpy [required: >=1.9.0, installed: 1.14.0]
  - python-dateutil [required: >=2, installed: 2.6.1]
    - six [required: >=1.5, installed: 1.11.0]
  - pytz [required: >=2011k, installed: 2018.3]
```

Note on deps

### Memory/System requirements

Note on system requirements

## Details on Results...

For the given MeSH term, PhenoX will attempt to generate all subclassification of the condition based on the available expression data in the NCBI Gene Expression Omnibus. Each subclass will have two major features, a word cloud of phenotypic terms associated with the seach tem via a literature search, and a gene expression profile of differentially expressed genes.

### Subsection A, Word Cloud

![alt text](https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/PhenoX.png "PhenoX Logo")

Something about wordcloud

### Subsection B, Expression Profiles

![alt text](https://github.com/NCBI-Hackathons/phenotypeXpression/blob/master/PhenoX.png "PhenoX Logo")

Something about heatmap

## Authors

Agnes Bao  
Robert R Butler III  
Huaiying Lin  
Subhajit Sengupta  
Lucy Lu Wang  

## Citation

Citation

```
@article{key ,
	author = {},
	title = {},
	journal = {F1000Research},
	volume = {},
	year = {2018},
	pages = {},
	doi = {},
	url = {}
}
```


## License

Copyright (C) 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License v3 as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

## References

<sup>1</sup>    
<sup>2</sup>    
<sup>3</sup>    
<sup>4</sup>    
<sup>5</sup>    
