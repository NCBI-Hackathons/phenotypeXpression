import os
import sys
from typing import List, Dict, Tuple
from phenox.paths import PhenoXPaths
from phenox.mesh_lookup import MeshSearcher
from phenox.geo_data import GEOQuery
from phenox.pubmed import Pubmed


# class for linking differential gene expression to disease
class PhenoX:
    def __init__(self, email: str, query_str: str) -> None:
        """
        Initialize class
        :param gene_list:
        """
        self.paths = PhenoXPaths()
        self.query_str = query_str
        self.email = email

    def _get_best_mesh_term(self) -> Tuple:
        """
        Retrieve the best MeSH term from search query
        :return:
        """
        mesh = MeshSearcher()
        mesh_entry = mesh.lookup(self.query_str)
        return mesh_entry['name'], [mesh.mesh[c]['name'] for c in mesh_entry['children']]

    def _get_geo_datasets(self, email: str, mesh_term: str) -> Tuple:
        """
        Given a MeSH term, fetch GEO datasets and corresponding PubMed IDs
        :param email:
        :param mesh_term:
        :return:
        """
        sys.stdout.write("Retrieving matching GEO datasets...\n")
        geo = GEOQuery(email=email)
        pubmed_ids, gene_freq, gene_dict = geo.get_all_geo_data(mesh_term)
        return pubmed_ids, gene_freq, gene_dict

    def _fetch_pubmed_abstracts(self, pubmed_ids: List) -> Dict:
        """
        Retrieve Pubmed abstracts from list of Pubmed IDs in the GEO data dict
        :param geo_data_dict:
        :return:
        """
        sys.stdout.write("Retrieving matching PubMed abstracts...\n")
        pubmed = Pubmed(self.email, pubmed_ids)
        term_freq = pubmed.get_term_frequencies()
        return term_freq

    def _analyze_aggregate_geo_data(self, geo_data_dict: Dict) -> Dict:
        """
        Analyze aggregate gene expression data
        :param geo_data_dict:
        :return:
        """
        sys.stdout.write("Analyzing aggregating gene expression data...\n")
        return dict()

    def _merge_and_combine(self, term_freq: Dict, clusters: Dict) -> Dict:
        """
        Merge term frequency data into cluster data
        :param term_freq:
        :param clusters:
        :return:
        """
        sys.stdout.write("Merging literature and gene expression data...\n")
        return clusters

    def _visualize(self, clusters: Dict) -> None:
        """
        Visualize clusters and labels
        :param clusters:
        :return:
        """
        return

    def subtype(self):
        """
        Run pipeline
        :return:
        """
        # get best mesh term from user query
        mesh_term, mesh_children = self._get_best_mesh_term()

        # retrieve GEO datasets and pubmed ids using MeSH disease term
        pubmed_ids, geo_gene_freq, geo_gene_dict = self._get_geo_datasets(self.email, mesh_term)

        # get pubmed abstracts
        term_frequency = self._fetch_pubmed_abstracts(pubmed_ids)

        # perform gene expression aggregation analysis
        geo_clusters = self._analyze_aggregate_geo_data(geo_gene_dict)

        # combine gene expression and phenotype data
        geo_clusters = self._merge_and_combine(term_frequency, geo_clusters)

        # visualize data
        self._visualize(geo_clusters)