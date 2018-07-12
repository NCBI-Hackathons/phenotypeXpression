import Entrez from Bio

#### Take a test case of Psoriasis in human samples with differentially expressed genes
print('Getting geoprofiles from NCBI: term = "up down genes"[filter] AND "Psoriasis"[Mesh] AND "Homo sapiens"[Organism]')
handle = Entrez.esearch(db="geoprofiles", term='"up down genes"[filter] AND "Psoriasis"[Mesh] AND "Homo sapiens"[Organism]', idtype="acc",retmax=9999)

#### Read the output
idrecord = Entrez.read(handle)

#### take a sample of geoprofiles
profileid = idrecord['IdList'][0:5000]

#### get docsums of the geoprofiles to exact gds:{gene:freq}
print("Getting docsums of geoprofiles... ")
query_results = Entrez.read(Entrez.efetch(db='geoprofiles', id=profileid, rettype='docsum'))

### output gdsdict and overall gene frequency
gdsdict,genefreq = gdsdict_from_profile(query_results)
print("Getting GEOdataset IDs from GEOprofiles... ")

### output gene ID => name table
genedict = genedict_from_profile(query_results)
print("GeneID to GeneName Table generated... ")

### output pid
#### get gdsIDs
gdslist = list(gdsdict.keys())
#### elink gds to pid
pidlist = Entrez.read(Entrez.elink(id=gdslist, db='pubmed', dbfrom='gds', linkname='gds_pubmed'))
#### extract pubmed IDs
pids = [el['LinkSetDb'][0]['Link'][0]['Id'] for el in pidlist]
print("Getting Pubmed IDs from GeoDataset... ")

