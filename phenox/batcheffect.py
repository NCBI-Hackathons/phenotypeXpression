import numpy as np
import pandas as pd
import scipy.stats as stats

from phenox.paths import PhenoXPaths


class BatchEffect:
    def __init__(self, clusters):
        self.paths = PhenoXPaths()
        self.clusters = clusters
        self.num_clusters = len(clusters)
    
    #from geneName.py -> to be sent to geo_data.py
    def meta_from_gds(self, gds_list: list) -> Dict:
        """
        Get meta information, such as gds submission time, n_samples, platform from gds.
        gds_list is from gdsdict.keys()
        - This can replace the get_pubmed_ids()
        concerns:
        >> do not know how to get rid of the IntergerElement
        >> now the publication time is in seconds, we can convert to days or years?
        """
        qout = Entrez.read(Entrez.esummary(db="gds", id=",".join(gds_list)))
        return {sm['Id']:[sm['n_samples'],sm['PubMedIds'],sm['GPL'],datetime.strptime(sm['PDAT'],'%Y/%m/%d').timestamp()] for sm in qout}

    # For chi-squared test of platform types different from overall
    def generate_t_test(self, meta_dict):
        """
        Generate t-test for sample_n/date batch effects
        :param meta_dict dict of meta_from_gds:
        :adds chi-square stat & pvalue to meta_dict:
        :return clust_stats dict:
        """
        clust_stats = {}
        
        for cluster in self.clusters:
            clust_frame = pd.DataFrame()
            ####do something to count categories for platform
            test_stat, p_val = stats.chisquare(f_obs=clust_count, f_exp=overall_count)
            clust_stats[cluster] = {'GPL': [chi_sq, p_val]}
            if p_val <= 0.05:
                print('{} is significantly different from overall platform distribution'.format(cluster))
        return clust_stats
    
    
    
        # num_columns = 3
        # num_rows = int(math.ceil(num_clusters/num_columns))
        # fig, axes = plt.subplots(num_rows, num_columns, figsize=(18, 5*num_rows), dpi= 80)
        # fig.subplots_adjust(hspace=0, wspace=0)

        # single_row = False
        # if num_rows == 1:
        #     single_row = True

        # turn axes off for everything and remove ticks
        # if single_row:
        #     for i in range(num_columns):
        #         axes[i].get_xaxis().set_ticks([])
        #         axes[i].get_yaxis().set_ticks([])
        #         axes[i].axis('off')
        # else:
        #     for i in range(num_rows):
        #         for j in range(num_columns):
        #             axes[i, j].get_xaxis().set_ticks([])
        #             axes[i, j].get_yaxis().set_ticks([])
        #             axes[i, j].axis('off')
        # 
        # clusters_sort = [(k, v) for k, v in clusters.items()]
        # clusters_sort.sort(key=lambda x: x[0])
        # 
        # # plot clusters
        # for i, cluster in enumerate(clusters_sort):
        #     cluster_name, freqs = cluster
        #     row = int(math.floor(i / num_columns))
        #     col = i % num_columns
        #     wordcloud = WordCloud(background_color="white").generate_from_frequencies(freqs)
        #     if single_row:
        #         axes[col].imshow(wordcloud, interpolation="bilinear")
        #         axes[col].set_title('{}'.format(cluster_name))
        #         axes[col].axis('on')
        #     else:
        #         axes[row, col].imshow(wordcloud, interpolation="bilinear")
        #         axes[row, col].set_title('{}'.format(cluster_name))
        #         axes[row, col].axis('on')
        # 
        # plt.savefig(output_file, bbox_inches='tight')
        # print('Wordcloud saved to %s' % output_file)
