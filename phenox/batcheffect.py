import os
import numpy as np
import pandas as pd
import scipy.stats as stats
from typing import List, Dict, Tuple
from collections import Counter, defaultdict

from phenox.paths import PhenoXPaths


class BatchEffect:
    def __init__(self, cluster_dict: Dict, meta_dict: Dict, outprefix: str) -> None:
        """
        Initialize class
        :param cluster_dict: key=cluster_ids, value=list of GDS ids
        :param meta_dict: key=gds_ids, value=list of gds metrics (n_samples, dates, GPLs)
        :param meta_list: index list of gds categories
        """
        self.paths = PhenoXPaths(outprefix)
        self.clusters = cluster_dict
        self.num_clusters = len(cluster_dict)
        self.meta_dict = meta_dict
        self.meta_list = ('Sample N', 'Submission Age', 'GPL')
    
    def _total_stats(self) -> List:
        """
        Generate total distributions for each gds metadata value
        :return total_stats: list of lists (distributions of all gds for each
            set: n_samples, dates, GPLs)
        """
        total_stats = [[] for i in range(3)]
        [[total_stats[i].append(self.meta_dict[gds][i]) for gds in self.meta_dict] for i in range(3)]
        return total_stats

    def _generate_ks_test(self, meta_value: int, total_dist: List, clust_stats=None) -> Dict:
        """
        Kolmogrov-Smirnov test for cluster vals and total vals from same dist
        :param meta_value: integer index for gds metric
        :param total_dist: list distribution of gds metric for ONE gds set
        :param clust_stats: pass mutable dict for recursive additions
            key=cluster_ids, value=ordered list of stats for all gds metrics
            (KS/chi-stat, p-value, mean/mode, s.d./category_n)
        :return:
        """
        # Create if empty
        if clust_stats is None:
            clust_stats = defaultdict(list)
        
        #total output values
        total_array = np.array(total_dist)
        clust_stats['Overall'].extend([None, None, np.mean(total_array), np.std(total_array, ddof=1)])    
        
        #per cluster outputs
        for cluster_id, cluster_list in self.clusters.items():
            cluster_dist = [self.meta_dict[gds][meta_value] for gds in cluster_list]
            clust_array = np.array(cluster_dist)
            KS_stat, p_val = stats.ks_2samp(clust_array, total_array)
            clust_stats[cluster_id].extend([KS_stat, p_val, np.mean(clust_array), np.std(clust_array, ddof=1)])
            if p_val <= 0.05:
                print('{}\'s {} sample is drawn from a significantly different distribution than the overall sample'.format(cluster_id, self.meta_list[meta_value]))
        return
    
    # For chi-squared test of platform types different from overall
    def _generate_chisq_test(self, total_dist: List, clust_stats: Dict) -> Dict:
        """
        Chi Squared test for cluster distribution independence from total dist
        :param total_dist: list distribution of gds metric for total gds set
        :param clust_stats: pass mutable dict for recursive additions
            key=cluster_ids, value=ordered list of stats for all gds metrics
            (KS/chi-stat, p-value, mean/mode, s.d./category_n)
        :return:
        """
        # total distribution to counts of categories
        total_frame = pd.DataFrame.from_dict(Counter(total_dist), orient='index')
        clust_stats['Overall'].extend([None, None, total_frame[0].idxmax(), len(total_frame.index)]) 
        
        # per cluster counts
        for cluster_id, cluster_list in self.clusters.items():
            cluster_dist = [self.meta_dict[gds][2] for gds in cluster_list] # GPL index in meta_list is 2
            clust_frame = pd.DataFrame.from_dict(Counter(cluster_dist), orient='index')
            chi_frame = clust_frame.merge(total_frame, how='right', right_index=True, left_index=True).fillna(0)
            chi_stat, p_val = stats.chisquare(f_obs=chi_frame['0_x'], f_exp=chi_frame['0_y'])
            clust_stats[cluster_id].extend([chi_stat, p_val, clust_frame[0].idxmax(), len(clust_frame.index)])
            if p_val <= 0.05:
                print('{}\'s platform distribution is significantly different from the overall distribution'.format(cluster_id))
        return

    def _stat_outfile(self, clust_stats: Dict) -> None:
        """
        Generate output tab separated text file
        :param clust_stats: pass mutable dict for recursive additions
            key=cluster_ids, value=ordered list of stats for all gds metrics
            (KS/chi-stat, p-value, mean/mode, s.d./category_n)
        :param outprefix: 
        :return:
        """
        # define order for column labels
        stat_list = ('KS-stat', 'P-value', 'Mean', 'Std Dev', 'Chi-stat', 'P-value', 'Mode', 'Category N')
        columns = [' '.join([i, j]) for i in (self.meta_list[0], self.meta_list[1]) for j in stat_list[:4]]
        columns.extend([' '.join([self.meta_list[2], j]) for j in stat_list[4:]])
        
        # transfer dict to dataframe
        out_tbl = pd.DataFrame.from_dict(clust_stats, orient='index', columns=columns)
        outfile = os.path.join(self.paths.output_dir, '{}_cluster_stats.txt'.format(self.paths.outprefix))
        out_tbl.to_csv(outfile, sep='\t', na_rep='.', index_label='Cluster ID')
        return

    def cluster_stats(self):
        """
        Run pipeline
        :return:
        """
        total_stats = self._total_stats()
        clust_stats = defaultdict(list)
        
        # For n_samples and date, populate clust_stats recursively
        [self._generate_ks_test(i, total_stats[i], clust_stats) for i in range(2)]
        
        # For GPL, chi squared stats in clust_stats
        self._generate_chisq_test(total_stats[2], clust_stats)
        
        #generate outfile
        self._stat_outfile(clust_stats)