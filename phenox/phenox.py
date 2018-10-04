import os
import sys
import math
import tqdm
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
from wordcloud import WordCloud

from phenox.paths import PhenoXPaths
from phenox.mesh_lookup import MeshSearcher
from phenox.geo_data import GEOQuery
from phenox.pubmed import Pubmed
from phenox.wordcloud import WordcloudPlotter
import phenox.utils.base_utils as base_utils

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

    def _get_geo_datasets(self, mesh_term: str) -> Tuple:
        """
        Given a MeSH term, fetch GEO datasets and corresponding PubMed IDs
        :param mesh_term:
        :return:
        """
        sys.stdout.write("Retrieving matching GEO datasets...\n")
        geo = GEOQuery(term=mesh_term, email=self.email)
        gene_dict, pubmed_dict, gds_dict, cluster_dict = geo.get_all_geo_data(mesh_term)
        return geo, gene_dict, pubmed_dict, gds_dict, cluster_dict

    def _fetch_pubmed_abstracts(
            self,
            geo: GEOQuery,
            mesh_term: str,
            gene_dict: Dict,
            pubmed_dict: Dict,
            gds_dict: Dict,
            cluster_dict: Dict
        ) -> Dict:
        """
        Retrieve Pubmed abstracts from list of Pubmed IDs in the GEO data dict
        :param geo: GEOQuery class used for querying
        :param gene_dict: dictionary key=gene_id, value=gene_name
        :param pubmed_dict: key=cluster_ids, value=list of pubmed ids
        :param gds_dict: key=GDS id, value=list of gene ids
        :param cluster_dict: key=cluster_ids, value=list of GDS ids
        :return:
        """

        # function to map GDS ids to a set of gene ids
        map_gds_to_gid = lambda gds_ids: list(set(base_utils.flatten(
            [gds_dict[gds_id] for gds_id in gds_ids]
        )))

        # function to map gene ids to gene names
        map_gid_to_gname = lambda gene_ids: [
            gene_dict[gid] for gid in gene_ids if gid in gene_dict
        ]

        sys.stdout.write("Retrieving matching PubMed abstracts...\n")

        # initialize pubmed clusters
        pubmed = Pubmed(self.email)
        pubmed_ids = pubmed.get_clusters(pubmed_dict, cluster_dict)

        sys.stdout.write("Retrieving matching PubMed abstracts for genes...\n")

        # fetch more pubmed ids based on gene name searches
        cluster_to_gene_name = {
            k: map_gid_to_gname(map_gds_to_gid(gds_list)) for k, gds_list in cluster_dict.items()
        }

        for clust, gene_names in cluster_to_gene_name.items():
            query_terms = pubmed.construct_query_terms(mesh_term['name'], gene_names)
            pubmed_clust = []
            sys.stdout.write('Querying PubMed for genes from {}\n'.format(clust))
            for query_term in tqdm.tqdm(query_terms, desc='PubMed batches'):
                pubmed_clust += geo.get_ncbi_docsum(mesh_term['name'], "pubmed", query_term)
            pubmed_ids[clust] += base_utils.flatten(pubmed_clust)

        # perform wordcloud NER
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
        geo, gene_dict, pubmed_dict, gds_dict, cluster_dict = self._get_geo_datasets(
            mesh_term['name']
        )

        # NER and count term frequency in pubmed clusters
        term_frequency = self._fetch_pubmed_abstracts(
            geo, mesh_term, gene_dict, pubmed_dict, gds_dict, cluster_dict
        )

        # visualize everything
        self._visualize(term_frequency)