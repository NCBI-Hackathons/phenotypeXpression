# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 15:34:02 2018

@author: Xiaojun
"""

import sys
import tqdm
import Bio.Entrez as Entrez
import spacy
import os
from spacy_lookup import Entity
import pickle
from collections import Counter
from phenox.paths import PhenoXPaths
import pandas as pd


# class for retrieving pubmed abstracts and finding disease/phenotype entities
class Pubmed:
    def __init__(self, email):
        Entrez.email = email
        self.paths = PhenoXPaths()
        self.pmid_abstracts = dict()
        # disease and human phenotype NER
        self.pmid_dner = {}
        # raw entity text
        self.pmid_ent_text = {}
        self.pmid_clustered = [] # a list of list: cluster of pmids
        self.dner_cluster = {}
        self.total_dner = []

        self.nlp = spacy.load('en')
        self.id2kw = pickle.load(open(os.path.join(self.paths.data_dir, 'id2kw_dict.pkl'), 'rb'))
        self.kw2id = pickle.load(open(os.path.join(self.paths.data_dir, 'kw2id_dict.pkl'), 'rb'))
        entity = Entity(keywords_list=list(self.kw2id.keys()), label='DO/HPO')
        self.nlp.add_pipe(entity, last=True)
        
    # fetch abstract from a single pmid    
    def fetch_abstract(self, pmid):
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        record = Entrez.read(handle)
        try:
            abstract = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        except:
            abstract = ''

        self.pmid_abstracts[pmid] = str(abstract)
        
    def extract_DNER(self, pmid):
        self.pmid_ent_text[pmid], self.pmid_dner[pmid] = self.find_DNER(self.pmid_abstracts[pmid])
    
    # find disease and phenotype NER
    def find_DNER(self, abstract_text):
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
    def cluster_count(self, cluster):
        '''
        cluster dtype: list of list
        '''
        n_cluster = 0
        for c in cluster:
            dner_cluster_list = []
            for pmid in c:
                dner_cluster_list.extend(self.pmid_dner[pmid])
            self.dner_cluster[n_cluster] = Counter(dner_cluster_list)
            n_cluster += 1
        return n_cluster

    # organize pmid into clusters from cluster analysis output (csv file)
    def get_cluster_from_csv(self, gdsid_pmid_map, csvfilepath):
        """
        Change read_csv path if csvfilepath changes in the future
        """
        cluster_res = pd.read_csv(csvfilepath)
        cluster_names = cluster_res['cluster'].unique()
        for n in cluster_names:
            c_pmid = []
            for g in cluster_res['gdsid'][cluster_res['cluster'] == n]:
                try:
                    c_pmid.append(gdsid_pmid_map[str(g)])
                except KeyError:
                    sys.stdout.write('Missing Pubmed key: {}\n'.format(g))
                    continue
            self.pmid_clustered.append(c_pmid)
        return self.pmid_clustered

    def get_term_frequencies(self, pmid_list):
        """
        Get all term frequencies
        :return:
        """
        for pmid in tqdm.tqdm(pmid_list):
            self.fetch_abstract(pmid)
            self.extract_DNER(pmid)

        freq_list = Counter(self.total_dner)
        return freq_list