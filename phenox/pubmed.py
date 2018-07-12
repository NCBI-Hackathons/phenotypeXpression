# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 15:34:02 2018

@author: Xiaojun
"""

import Bio.Entrez as Entrez
import spacy
from spacy_lookup import Entity
import pickle

Entrez.emal = '390943746@qq.com'

'''
Takea a list of PubMed IDs, fetch abstract text and find DO and HPO entities
from the text 
'''
class pubmed:
    def __init__(self, pmid_list):
        '''
        input: list of PubMed IDs
        '''
        self.pmid_list = pmid_list
        self.pmid_abstracts = {}
        self.nlp = spacy.load('en')
        self.id2kw = json.load(open('../data/doid.json'))
        self.id2kw.update(json.load(open('../data/hp.json')))
        entity = Entity(keywords_list=['python', 'java platform'])
        
    # fetch abstract from a single pmid    
    def fetch_abstract(self,pmid):
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        record = Entrez.read(handle)
        abstract = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        self.pmid_abstracts[pmid] = abstract
    
    # find disease and phenotype NER
    def find_DNER(self,abstract_text):
        pass
    
abstract = u"The association studies of several immune-diseases with the 3' Regulatory Region 1 (3'RR1) increased interest on the role that the region plays in the immune-regulation. The 3'RR1 is a polymorphic region on human chromosome 14q32, acting as a cis-regulative element on the Immunoglobulin constant-gene locus recently considered as super-enhancer. The human 3'RR1 share large sequences with its paralogous 3'RR2, at high level of similarity. Thus, a focused investigation was necessary to discriminate each one of the duplicated components of the two regions and its specific contribution to the immunologic phenotype. One of the duplicated elements is the hs1.2 enhancer. The 3'RR1 alleles of this enhancer were demonstrated to play a role in autoimmune diseases, including Psoriasis. We sequenced a specific region internal to the 3'RR1 in hs1.2 homozygous subjects, to detect SNPs associated to the main alleles of the enhancer. We identified two alternative nine-SNPs haplotypes strictly linked to the allele *1 and *2 of hs1.2, that could be used as markers to further investigate the region and associations to pathology. Finally, we identified two haplotypes, namely E2A1 and E2A2, that strongly support the hypothesis of a relevant effect of the rs35216181 in the onset of Psoriasis when the *2 allele is present." 
nlp = spacy.load('en')
entity = Entity(keywords_file='../data/kwlist.txt',label='DO')
nlp.add_pipe(entity, last=True)

doc = nlp(abstract)
        
    
    
        
    
        
        