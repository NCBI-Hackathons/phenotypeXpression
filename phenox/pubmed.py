# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 15:34:02 2018

@author: Xiaojun
"""


import Bio.Entrez as Entrez
import spacy
from spacy_lookup import Entity
import pickle
from collections import Counter

Entrez.email = '390943746@qq.com'

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
        # disease and human phenotype NER
        self.pmid_dner = {}
        # raw entity text
        self.pmid_ent_text = {}
        self.dner_cluster = {}
        self.total_dner = []

        self.nlp = spacy.load('en')
        self.id2kw = pickle.load(open('../data/id2kw_dict.pkl','rb'))
        self.kw2id = pickle.load(open('../data/kw2id_dict.pkl','rb'))
        entity = Entity(keywords_list=list(self.kw2id.keys()),label='DO/HPO')
        self.nlp.add_pipe(entity, last=True)
        
    # fetch abstract from a single pmid    
    def fetch_abstract(self,pmid):
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        record = Entrez.read(handle)
        try:
            abstract = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        except:
            abstract = ''
        self.pmid_abstracts[pmid] = str(abstract)
        
    def extract_DNER(self,pmid):
        self.pmid_ent_text[pmid], self.pmid_dner[pmid] = self.find_DNER(self.pmid_abstracts[pmid])
    
    # find disease and phenotype NER
    def find_DNER(self,abstract_text):
        kw = []
        dner = []
        doc = self.nlp(abstract_text)
        for ent in doc.ents:
            if ent.label_=='DO/HPO' and ent.__getitem__(0).is_stop==False:
                kw.append(ent.text)
                # convert recognized keyword to the main term (first item) in DO/HPO
                dner.append(self.id2kw[self.kw2id[ent.text.lower()]][0])
        self.total_dner.extend(dner)
        return Counter(kw), dner
    
    # count word frequencies based on clustering results
    def cluster_count(self,cluster):
        '''
        cluster dtype: list of list
        '''
        n_cluster = 0
        for c in cluster:
            dner_cluster_list = []
            for pmid in c:
                dner_cluster_list.extend(self.pmid_dner[pmid])
            self.dner_cluster[n_cluster] = Counter(dner_cluster_list)
            n_cluster+=1
            
                
        

# pmid list for psoriasis
psoriasis = ['24646743', 
             '22479649', 
             '25129481', 
             '25129481', 
             '23633458', 
             '24391825', 
             '23771123', 
             '22277938', 
             '22348003', 
             '21483750', 
             '22479649', 
             '27667537', 
             '23407402', 
             '21388663', 
             '22908096', 
             '22677045', 
             '20829794', 
             '22752307', 
             '20688981', 
             '19718476', 
             '18648529', 
             '19052557', 
             '18716044', 
             '17277128', 
             '16858420', 
             '17020965', 
             '17947518', 
             '16885358', 
             '16618722', 
             '16505361', 
             '17075716', 
             '11121445']

test = pubmed(psoriasis)
for pmid in test.pmid_list:
    test.fetch_abstract(pmid)
    test.extract_DNER(pmid)
    
test.freq_list = Counter(test.total_dner)
        
    
        
        