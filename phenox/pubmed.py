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
from geo_data import GEOQuery


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
        abs_str = []
        # abstract text
        try:
            abstract = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
            for a in abstract:
                abs_str.append(a)
        except:
            pass
        # title
        try:
            title = record['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
            abs_str.append(title)
        except:
            pass
        # keyword list
        try:
            kwls = record['PubmedArticle'][0]['MedlineCitation']['KeywordList'][0]
            for w in kwls:
                abs_str.append(w)
        except:
            pass

        self.pmid_abstracts[pmid] = ' '.join(abs_str).strip()
        
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
                if ent.text.lower() in self.kw2id:
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
    def get_clusters(self, gdsid_pmid_map, cluster_dict):
        """
        Change read_csv path if csvfilepath changes in the future
        """
        pmid_clustered = dict()

        for cluster_name, cluster_genes in cluster_dict.items():
            c_pmid = []
            for gene_id in cluster_genes:
                try:
                    c_pmid.append(gdsid_pmid_map[gene_id])
                except KeyError:
                    sys.stdout.write('Missing Pubmed key: {}\n'.format(g))
                    continue
            pmid_clustered[cluster_name] = c_pmid

        print(pmid_clustered)
        return pmid_clustered

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
    
    
class PubmedQuery(GEOQuery):
    def get_ncbi_docsum(self, mesh_term, gene_name_list) -> List:
        """
        Get document summaries from NCBI database
        :param mesh_term:
        :param database:
        :return:
        """
        start = self.timing_tool()
        staging_list = []
        query_results = []
        
        gene_list = ''
        for i in range(len(gene_name_list)-1):
            gene_list+='"{}" OR '.format(gene_name_list[i])
        gene_list+='"{}"'.format(gene_name_list[-1])
        search_term = '"{}"[MeSH Terms] AND ({})'.format(mesh_term, gene_list)
        if len(search_term)>4000:
            print('Query term is longer than 4000 charaters!')
        
        # First, esearch using mesh_term, and keep results on ePost history server
        webenv1, query_key1, count1 = self.post_ncbi(
            'esearch', db='pubmed', usehistory='y',
            term=search_term)

        # Then, efetch in batches of efetch_batch
        self.batch_ncbi('efetch', staging_list, count1, db=database,
                        rettype='docsum', retmax=self.efetch_batch, webenv=webenv1,
                        query_key=query_key1)

        # append to passed list
        [query_results.append(batch) for batch in staging_list]

        # stop timing and log
        stop = self.timing_tool()
        logging.info('Download time: {} min, Batches run -> {}'
                     .format(((stop - start) / 60), len(query_results)))

        return query_results