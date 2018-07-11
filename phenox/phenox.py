import os
import sys
from typing import List, Dict, Tuple
from phenox.paths import PhenoXPaths


# class for linking differential gene expression to disease
class PhenoX:
    def __init__(self, query_str: str) -> None:
        """
        Initialize class
        :param gene_list:
        """
        self.paths = PhenoXPaths()
        self.query_str = query_str

    def _get_best_mesh_term(self) -> Tuple:
        """
        Retrieve the best MeSH term from search query
        :return:
        """
        sys.stdout.write("Mapping to closest MeSH term...\n")
        return "", list()

    def _get_geo_datasets(self, mesh_term: str) -> Dict:
        """
        Given a MeSH term, fetch GEO datasets and corresponding PubMed IDs
        :param mesh_term:
        :return:
        """
        sys.stdout.write("Retrieving matching GEO datasets...\n")
        return dict()

    def _fetch_pubmed_abstracts(self, geo_data_dict: Dict) -> Dict:
        """
        Retrieve Pubmed abstracts from list of Pubmed IDs in the GEO data dict
        :param geo_data_dict:
        :return:
        """
        sys.stdout.write("Retrieving matching PubMed abstracts...\n")
        return dict()

    def _extract_phenotype_disease(self, pubmed_dict: Dict) -> Dict:
        """
        Extract phenotyp and disease terms and sort by frequency
        :param pubmed_dict:
        :return:
        """
        sys.stdout.write("Extracting phenotypes from abstracts...\n")
        return dict()

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
        mesh_term, mesh_tree = self._get_best_mesh_term()

        # retrieve GEO datasets and pubmed ids using MeSH disease term
        geo_data_dict = self._get_geo_datasets(mesh_term)

        # get pubmed abstracts
        pubmed_dict = self._fetch_pubmed_abstracts(geo_data_dict)

        # get phenotype and disease terms from pubmed abstracts
        term_frequency = self._extract_phenotype_disease(pubmed_dict)

        # perform gene expression aggregation analysis
        geo_clusters = self._analyze_aggregate_geo_data(geo_data_dict)

        # combine gene expression and phenotype data
        geo_clusters = self._merge_and_combine(term_frequency, geo_clusters)

        # visualize data
        self._visualize(geo_clusters)