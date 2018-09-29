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
    def __init__(self, email: str, query_str: str, outprefix: str) -> None:
        """
        Initialize class
        :param gene_list:
        """
        self.paths = PhenoXPaths(outprefix)
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
        geo = GEOQuery(outprefix=self.paths.outprefix, term=mesh_term, email=self.email)
        pubmed_dict = geo.get_all_geo_data(mesh_term)
        return pubmed_dict

    def _fetch_pubmed_abstracts(self, pubmed_dict: Dict, cluster_dict: Dict) -> Dict:
        """
        Retrieve Pubmed abstracts from list of Pubmed IDs in the GEO data dict
        :param pubmed_dict: key=cluster_ids, value=list of pubmed ids
        :param cluster_dict: key=cluster_ids, value=list of gene ids
        :return:
        """
        sys.stdout.write("Retrieving matching PubMed abstracts...\n")

        pubmed = Pubmed(self.email)
        pubmed_ids = pubmed.get_clusters(pubmed_dict, cluster_dict)

        wordcloud_data = dict()

        for cluster_name, cluster_ids in pubmed_ids.items():
            term_freq = pubmed.get_term_frequencies(cluster_ids)
            wordcloud_data[cluster_name] = term_freq

        return wordcloud_data

    def _visualize(self, clusters: Dict) -> None:
        """
        Visualize clusters and labels
        :param clusters:
        :return:
        """
        output_file = os.path.join(
            self.paths.output_dir,
            '{}_GDS_wordcloud.png'.format(self.query_str.replace(' ', '-'))
        )
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
        pubmed_dict, cluster_dict = self._get_geo_datasets(mesh_term['name'])

        # NER and count term frequency in pubmed clusters
        term_frequency = self._fetch_pubmed_abstracts(pubmed_dict, cluster_dict)

        # visualize everything
        self._visualize(term_frequency)