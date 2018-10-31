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
import json
import pickle
from spacy_lookup import Entity
from collections import Counter
from typing import List, Dict, Tuple

from phenox.paths import PhenoXPaths


# class for retrieving pubmed abstracts and finding disease/phenotype entities
class Pubmed:
    def __init__(self, email: str, outprefix: str):
        Entrez.email = email
        self.paths = PhenoXPaths(outprefix)
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

        # read synonyms from HGNC
        with open(os.path.join(self.paths.data_dir, 'hgnc_synonyms.json'), 'r') as f:
            hgnc_syn = f.read()
            self.hgnc = json.loads(hgnc_syn)
        
    # fetch abstract from a single pmid    
    def fetch_abstract(self, pmid: str) -> None:
        """
        Fetch abstract from a single pmid   
        """

        self.pmid_abstracts[pmid] = ''

        try:
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
        except:
            pass
        
    def extract_DNER(self, pmid: str) -> None:
        """
        Store dner count and dner associated with each pmid
        :param pmid:
        """
        self.pmid_ent_text[pmid], self.pmid_dner[pmid] = self.find_DNER(self.pmid_abstracts[pmid])
    
    # find disease and phenotype NER
    def find_DNER(self, abstract_text: str) -> Tuple:
        """
        Identify disease and phenotype NER from abstract text
        :param abstract_text:
        :return:
        """
        kw = []
        dner = []
        try:
            doc = self.nlp(abstract_text)
            for ent in doc.ents:
                if ent.label_=='DO/HPO' and ent.__getitem__(0).is_stop==False:
                    kw.append(ent.text)
                    if ent.text.lower() in self.kw2id:
                        # convert recognized keyword to the main term (first item) in DO/HPO
                        dner.append(self.id2kw[self.kw2id[ent.text.lower()]][0])
            self.total_dner.extend(dner)
            return Counter(kw), dner
        except KeyboardInterrupt:
            raise
        except Exception:
            return Counter([]), []

    # count word frequencies based on clustering results
    def cluster_count(self, cluster: List) -> int:
        """
        :param cluster: list (number of clusters) of list (member of each cluster)
        :return:
        """
        n_cluster = 0
        for c in cluster:
            dner_cluster_list = []
            for pmid in c:
                dner_cluster_list.extend(self.pmid_dner[pmid])
            self.dner_cluster[n_cluster] = Counter(dner_cluster_list)
            n_cluster += 1
        return n_cluster

    # organize pmid into clusters from cluster analysis output (csv file)
    def get_clusters(self, gdsid_pmid_map: Dict, cluster_dict: Dict) -> Dict:
        """
        Organize pmid into clusters from cluster analysis output
        :param gdsid_pmid_map:
        :param cluster_dict:
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

        return pmid_clustered

    def get_term_frequencies(self, pmid_list: List) -> Dict:
        """
        Get all term frequencies
        :return:
        """
        for pmid in tqdm.tqdm(pmid_list):
            self.fetch_abstract(pmid)
            self.extract_DNER(pmid)

        freq_list = Counter(self.total_dner)

        return freq_list

    def construct_query_terms(self, mesh_term: str, gene_name_list: List) -> List:
        """
        Construct a list of query terms to make sure each term does not exceed
        4000 characters. 
        :param mesh_term:
        :param gene_name_list: list of gene names
        :return query_term_list: a list of query strings.
        """
        query_term_list = []

        # initialize search term
        search_prefix = mesh_term + " AND "
        search_term = "("

        # generate gene name synonym list from HGNC
        gene_syn_list = []
        for gene_name in gene_name_list:
            gene_syn_list.append(gene_name)
            if gene_name in self.hgnc:
                gene_syn_list += self.hgnc[gene_name]['names'] + self.hgnc[gene_name]['aliases']

        # iterate through gene name list
        for gene_name in tqdm.tqdm(gene_syn_list):
            if len(search_prefix + search_term + gene_name) > 4000:
                # remove trailing " OR " and append ")"
                search_term = search_term[:-4] + ")"
                # append to query term list
                query_term_list.append(search_term)
                # reinitialize search term
                search_term = "("
            else:
                # append gene name to search term
                search_term += gene_name + " OR "

        # add final search term if present
        if search_term != "(":
            search_term = search_term[:-4] + ")"
            query_term_list.append(search_term)

        return query_term_list
