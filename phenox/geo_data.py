import os
import sys
import logging
import time
import tqdm
from typing import List, Dict, Tuple
from collections import defaultdict
import Bio.Entrez as Entrez
from urllib.error import HTTPError

import numpy as np
import pandas as pd
import pydendroheatmap as pdh
from datetime import datetime

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.conversion import Converter, localconverter
from rpy2.robjects.packages import importr
from rpy2 import rinterface

from phenox.paths import PhenoXPaths


# class for querying GEO databases
class GEOQuery:
    def __init__(self, outprefix, term, email, tool="phenotypeXpression", efetch_batch=5000, elink_batch=100):
        Entrez.email = email
        Entrez.tool = tool
        self.efetch_batch = efetch_batch
        self.elink_batch = elink_batch
        self.db = 'geoprofiles'

        self.paths = PhenoXPaths(outprefix)
        term_name = term.replace(' ', '-')
        self.hcluster_file = os.path.join(
            self.paths.output_dir, "{}_{}_hierarchical_clusters.pdf".format(self.paths.outprefix, term_name)
        )
        self.tree_file = os.path.join(
            self.paths.output_dir, "{}_{}_newick_tree.txt".format(self.paths.outprefix, term_name)
        )
        self.heatmap_file = os.path.join(
            self.paths.output_dir, "{}_{}_heatmap.pdf".format(self.paths.outprefix, term_name)
        )
        self.dist_graph_file = os.path.join(
            self.paths.output_dir, "{}_{}_dist_graph.pdf".format(self.paths.outprefix, term_name)
        )

    @staticmethod
    def timing_tool():
        return time.perf_counter()

    @staticmethod
    def http_attempts(query_type, **kwargs):
        """
        HTTP attempts
        :param query_type:
        :param kwargs:
        :return:
        """
        fetch_handle = None

        for attempt in range(4):  # HTTPError and retry three times
            try:
                fetch_handle = getattr(Entrez, query_type)(**kwargs)
            except ValueError as oops:
                logging.warning('Empty set: Likely total number = batch size')
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    logging.info("Attempt {} of 4".format(attempt + 1))
                    logging.warning(err)
                    time.sleep(15)
                elif err.code == 400:
                    logging.info("Attempt {} of 4".format(attempt + 1))
                    logging.warning(err)
                    sys.tracebacklimit = None
                    time.sleep(15)
                else:
                    raise
            except KeyboardInterrupt:
                raise KeyboardInterrupt()
            else:
                break
        else:
            logging.critical('Failed to connect to NCBI 4 times exiting...')
            sys.exit()
        return fetch_handle

    def post_ncbi(self, query_type, **kwargs):
        """
        Post request to NCBI
        :param query_type:
        :param kwargs:
        :return:
        """
        logging.debug('{} to ncbi using kwargs: {}'.format(query_type, kwargs))
        fetch_handle = self.http_attempts(query_type, **kwargs)
        query = Entrez.read(fetch_handle)
        fetch_handle.close()

        webenv = query['WebEnv']
        query_key = query['QueryKey']
        count = int(query['Count'])

        logging.debug('returned webenv: {} and query key: {}'
                      .format(webenv, query_key))
        return webenv, query_key, count

    # Utilized for batch efetch/esummary from ncbi using history server.
    def batch_ncbi(self, query_type, query_results, count1, **kwargs):
        # number in history list
        count = count1

        # sliding window of history list in batch size
        for start in range(0, count, self.efetch_batch):
            end = min(count, start + self.efetch_batch)
            logging.info("Going to download record {} to {}"
                         .format(start + 1, end))

            logging.debug('{} run with {}'
                          .format(query_type, dict(retstart=start, **kwargs)))
            fetch_handle = self.http_attempts(
                query_type, **dict(retstart=start, **kwargs)
            )
            query_results.append(Entrez.read(fetch_handle))
            fetch_handle.close()
        return

    def get_ncbi_docsum(self, mesh_term, database, query_term='"up down genes"[filter]') -> List:
        """
        Get document summaries from NCBI database
        :param mesh_term:
        :param database:
        :param query_term:
        :return:
        """
        start = self.timing_tool()
        staging_list = []
        query_results = []

        # First, esearch using mesh_term, and keep results on ePost history server
        webenv1, query_key1, count1 = self.post_ncbi(
            'esearch', db=database, usehistory='y',
            term='{} AND {}'.format(mesh_term, query_term)
        )

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

    def batch_local(self, query_type, id_list, **kwargs) -> Dict:
        """
        Batch local for elink,einfo queries w/ no retstart parameter
        :param query_type:
        :param id_list:
        :param kwargs:
        :return:
        """
        count = len(id_list)
        result_dict = dict()

        for start in range(0, count, self.elink_batch):
            end = min(count, start + self.elink_batch)
            logging.info("Looking up PMIDs for GDSs {} to {}"
                         .format(start + 1, end))

            # elink each batch
            fetch_handle = self.http_attempts(
                query_type, **dict(id=','.join(id_list[start:end]), **kwargs)
            )
            record = Entrez.read(fetch_handle)
            fetch_handle.close()
            # for each Link'Id', dig down through dict add to list
            for el in record:
                result_dict[el['IdList'][0]] = el['LinkSetDb'][0]['Link'][0]['Id']
        return result_dict

    def gdsdict_from_profile(self, query_results) -> Dict:
        """
        Generate GDS data dict from profile data
        using <gene name> as dict key
        :param query_results: Dict from ncbi query
        :return gds_dict: k=gds_id, v={k=geneName, v=gene freq}
        """
        gds_dict = dict()

        for batch in query_results:
            for docsum in tqdm.tqdm(batch['DocumentSummarySet']['DocumentSummary'],desc="GEO Datasets"):
                gdsid = docsum['GDS']
                gnamelist = docsum['geneName'].split("<:>")
                for gname in gnamelist:
                    if len(gname) > 0:
                        if gdsid not in gds_dict:
                            gds_dict[gdsid] = defaultdict(int)
                        gds_dict[gdsid][gname.upper()] += 1
        return gds_dict

    def meta_from_gds(self, gds_dict: Dict) -> Dict:
        """
        Get meta information, such as gds submission time, n_samples, platform from gds.
        :param gds_dict: k=gds_id, v={k=geneName, v=gene freq}
        :return: meta_dict: key=gds_ids, value=list of gds metrics (n_samples, dates, GPLs)
        """
        gds_list = list(gds_dict.keys())
        qout = Entrez.read(Entrez.esummary(db="gds", id=",".join(gds_list)))
        meta_dict = {sm['Id']:[sm['n_samples'],
                               datetime.strptime(sm['PDAT'],'%Y/%m/%d').timestamp(),
                               sm['GPL']]
                     for sm in qout}
        return meta_dict

    def get_pubmed_ids(self, data_dict) -> Dict:
        """
        Fetch PubMed terms from GDS data
        :param gds_dict: k=gds_id, v={k=geneName, v=gene freq}
        :return pids:
        """
        gds_list = list(data_dict.keys())
    
        # elink gds to pid
        pid_list = Entrez.read(
            Entrez.elink(
                id=gds_list, db='pubmed', dbfrom='gds', linkname='gds_pubmed'
            )
        )
    
        # extract pubmed IDs
        pids = {
            el['IdList'][0]: el['LinkSetDb'][0]['Link'][0]['Id'] for el in pid_list
        }
    
        return pids
    
    def gds_to_pd_dataframe(self, gds_dict: Dict):
        """
        Convert GDS to pandas dataframe
        :param gds_dict: k=gds_id, v={k=geneName, v=gene freq}
        :return pdd:
        """
        pdd = pd.DataFrame(list(gds_dict.values()))
        pdd.index = gds_dict.keys()
        pdd = pdd.fillna(0)
        pdd[pdd > 1] = 1
        col_select = pdd.columns[pdd.sum() > 1]
        pdd = pdd.loc[:, col_select]
        pdd = pdd[pdd.sum(axis=1) > 1]
        return pdd

    def _generate_dist_graph(self, gds_py):
        """
        Generate a distance plot from GDS matrix
        :param gds_py:
        :return:
        """
        exp_list = list(gds_py.index)
        gds_array = np.array(gds_py)
        gds_array -= gds_array.mean()

        similarities = euclidean_distances(gds_array)
        mds = manifold.MDS(n_components=2, max_iter=100, eps=1e-9,
                           dissimilarity="precomputed", n_jobs=1)
        pos = mds.fit(similarities).embedding_

        clf = PCA(n_components=2)
        pos = clf.fit_transform(pos)

        fig = plt.figure(1)
        ax = plt.axes([0., 0., 1., 1.])

        plt.scatter(pos[:, 0], pos[:, 1], color='navy', s=100, lw=0, label='MDS')
        for i, txt in enumerate(exp_list):
            ax.annotate(txt, (pos[i, 0] + 0.1, pos[i, 1] + 0.1), color='black')

        segments = [[pos[i, :], pos[j, :]] for i in range(len(pos)) for j in range(len(pos))]
        values = np.abs(similarities)
        lc = LineCollection(segments,
                            zorder=0, cmap=plt.cm.Blues,
                            norm=plt.Normalize(0, values.max()))
        lc.set_array(similarities.flatten())
        lc.set_linewidths(0.5 * np.ones(len(segments)))
        ax.add_collection(lc)

        xmin, xmax = plt.xlim()
        plt.xlim((xmin, xmax + 0.5))

        ymin, ymax = plt.ylim()
        plt.ylim((ymin - 1., ymax + 1.))

        fig.suptitle("Distances")
        fig.savefig(self.dist_graph_file, bbox_inches='tight')
        print("Distance graph written to {}".format(self.dist_graph_file))
        return

    def call_r_clustering(self, gds_py):
        """
        Call R script
        :param gds: dataframe with GDS data
        :return:
        """
        # load all R libraries
        base = importr("base")
        pvclust = importr("pvclust")
        graphics = importr("graphics")
        grdevices = importr("grDevices")
        ape = importr("ape")
        gplots = importr("gplots")

        # convert pandas df to R df
        with localconverter(default_converter + pandas2ri.converter) as cv:
            gds = pandas2ri.py2ri(gds_py)

        matrix = base.as_matrix(gds)
        mat_trans = matrix.transpose()

        # cluster over studies
        print("Clustering on Studies...")
        fit = pvclust.pvclust(mat_trans, nboot=5000, method_hclust="ward.D2", method_dist="euclidean")

        # write clustering output to pdf
        grdevices.pdf(self.hcluster_file, paper="a4")
        graphics.plot(fit)
        pvclust.pvrect(fit, alpha=.95)
        grdevices.dev_off()
        print("Clustering diagram written to {}".format(self.hcluster_file))

        # write heatmap with pvclust to pdf 
        grdevices.pdf(self.heatmap_file, paper="a4")
        gplots.heatmap_2(mat_trans,Rowv = fit,dendrogram = "col",col=gplots.bluered(100), labRow = False, trace="none")
        grdevices.dev_off()
        print("Heatmap written to {}".format(self.heatmap_file))


        # write cluster tree output to tree file
        hc = fit.rx2("hclust")
        tree = ape.as_phylo(hc)
        ape.write_tree(phy=tree, file=self.tree_file)
        print("Cluster tree written to {}".format(self.tree_file))

        # extract cluster membership
        pvp = pvclust.pvpick(fit)
        clusters = pvp.rx2("clusters")

        cluster_members = defaultdict(list)

        if clusters[0] == rinterface.NULL:
            print("WARNING: Only one cluster!")
            cluster_members['cluster0'] = list(gds_py.index)
        else:
            for i, clust in enumerate(clusters):
                clust_name = "cluster{}".format(i)
                for pmid in clust:
                    cluster_members[clust_name].append(pmid)

        # generate distance graph
        self._generate_dist_graph(gds_py)
        return cluster_members

    def get_all_geo_data(self, mesh_term: str) -> Tuple:
        """
        Link all functions together to retrieve GEO data
        :param mesh_term:
        :return:
        """
        # query GEO using MeSH term and retrieve datasets
        query_results = self.get_ncbi_docsum(mesh_term, self.db)

        # Map GEO datasets to genes names and count gene frequency
        gds_dict = self.gdsdict_from_profile(query_results)

        # Map GEO datasets to Pubmed IDs
        geo_to_pid_dict = self.get_pubmed_ids(gds_dict)

        # Retrieve gene names from GEO gene profiles
        # gid_to_gname_dict = self.genedict_from_profile(query_results)

        # Export GDS
        gds = self.gds_to_pd_dataframe(gds_dict)

        # Run clustering algorithm
        clusters = self.call_r_clustering(gds)
        
        # Batch effect data
        meta_dict = self.meta_from_gds(gds_dict)

        return geo_to_pid_dict, gds_dict, clusters, meta_dict

