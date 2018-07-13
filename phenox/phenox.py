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
        return mesh_entry, [mesh.mesh[c]['name'] for c in mesh_entry['children']]

    def _get_geo_datasets(self, mesh_term: str) -> Dict:
        """
        Given a MeSH term, fetch GEO datasets and corresponding PubMed IDs
        :param mesh_term:
        :return:
        """
        sys.stdout.write("Retrieving matching GEO datasets...\n")
        geo = GEOQuery(email=self.email)
        pubmed_dict = geo.get_all_geo_data(mesh_term)
        return pubmed_dict

    def _fetch_pubmed_abstracts(self, pubmed_dict: Dict) -> Dict:
        """
        Retrieve Pubmed abstracts from list of Pubmed IDs in the GEO data dict
        :param geo_data_dict:
        :return:
        """
        sys.stdout.write("Retrieving matching PubMed abstracts...\n")

        cluster_csv_file = os.path.join(self.paths.data_dir, 'memb.csv')

        pubmed = Pubmed(self.email)
        pubmed_ids = pubmed.get_cluster_from_csv(pubmed_dict, cluster_csv_file)

        wordcloud_data = dict()

        for i, id_list in enumerate(pubmed_ids):
            term_freq = pubmed.get_term_frequencies(id_list)
            wordcloud_data[i + 1] = term_freq

        return wordcloud_data

    def _visualize(self, clusters: Dict) -> None:
        """
        Visualize clusters and labels
        :param clusters:
        :return:
        """
        output_file = os.path.join(self.paths.output_dir, 'wordcloud.png')
        plotter = WordcloudPlotter()
        plotter.generate_wordclouds(clusters, output_file)
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
        pubmed_dict = self._get_geo_datasets(mesh_term['name'])

        # NER and count term frequency in pubmed clusters
        term_frequency = self._fetch_pubmed_abstracts(pubmed_dict)

        # visualize everything
        self._visualize(term_frequency)