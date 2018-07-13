import os
import sys
import math
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
from wordcloud import WordCloud

from phenox.paths import PhenoXPaths
from phenox.mesh_lookup import MeshSearcher
from phenox.geo_data import GEOQuery
from phenox.pubmed import Pubmed
from phenox.wordcloud import WordcloudPlotter

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

    def _get_geo_datasets(self, email: str, mesh_term: str) -> List:
        """
        Given a MeSH term, fetch GEO datasets and corresponding PubMed IDs
        :param email:
        :param mesh_term:
        :return:
        """
        sys.stdout.write("Retrieving matching GEO datasets...\n")
        geo = GEOQuery(email=email)
        pubmed_clusters = geo.get_all_geo_data(mesh_term)
        return pubmed_clusters

    def _fetch_pubmed_abstracts(self, pubmed_ids: List) -> Dict:
        """
        Retrieve Pubmed abstracts from list of Pubmed IDs in the GEO data dict
        :param geo_data_dict:
        :return:
        """
        sys.stdout.write("Retrieving matching PubMed abstracts...\n")

        wordcloud_data = dict()

        for i, id_list in enumerate(pubmed_ids):
            pubmed = Pubmed(self.email, id_list)
            term_freq = pubmed.get_term_frequencies()
            wordcloud_data[i + 1] = term_freq

        return wordcloud_data

    def _visualize(self, clusters: Dict) -> None:
        """
        Visualize clusters and labels
        :param clusters:
        :return:
        """
        plotter = WordcloudPlotter()
        plotter.generate_wordclouds(clusters)
        return

    def subtype(self):
        """
        Run pipeline
        :return:
        """
        # get best mesh term from user query
        mesh_term, mesh_children = self._get_best_mesh_term()

        # retrieve GEO datasets, generate clusters, visualize in R,
        # and output clustered pubmed abstracts
        pubmed_clusters = self._get_geo_datasets(self.email, mesh_term)

        # NER and count term frequency in pubmed clusters
        term_frequency = self._fetch_pubmed_abstracts(pubmed_clusters)

        # visualize everything
        self._visualize(term_frequency)