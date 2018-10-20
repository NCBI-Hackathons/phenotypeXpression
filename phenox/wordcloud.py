from wordcloud import WordCloud
import matplotlib.pyplot as plt
import math

from phenox.paths import PhenoXPaths


class WordcloudPlotter:
    def __init__(self, outprefix: str):
        self.paths = PhenoXPaths(outprefix)

    def generate_wordclouds(self, clusters, output_file):
        """
        Generate plot with subplots of wordclouds
        :param clusters:
        :return:
        """
        num_clusters = len(clusters)
        num_columns = 3
        num_rows = int(math.ceil(num_clusters/num_columns))
        fig, axes = plt.subplots(num_rows, num_columns, figsize=(18, 5*num_rows), dpi= 80)
        fig.subplots_adjust(hspace=0, wspace=0)

        single_row = False
        if num_rows == 1:
            single_row = True

        # turn axes off for everything and remove ticks
        if single_row:
            for i in range(num_columns):
                axes[i].get_xaxis().set_ticks([])
                axes[i].get_yaxis().set_ticks([])
                axes[i].axis('off')
        else:
            for i in range(num_rows):
                for j in range(num_columns):
                    axes[i, j].get_xaxis().set_ticks([])
                    axes[i, j].get_yaxis().set_ticks([])
                    axes[i, j].axis('off')

        clusters_sort = [(k, v) for k, v in clusters.items()]
        clusters_sort.sort(key=lambda x: x[0])

        # plot clusters
        for i, cluster in enumerate(clusters_sort):
            cluster_name, freqs = cluster
            row = int(math.floor(i / num_columns))
            col = i % num_columns
            wordcloud = WordCloud(background_color="white").generate_from_frequencies(freqs)
            if single_row:
                axes[col].imshow(wordcloud, interpolation="bilinear")
                axes[col].set_title('{}'.format(cluster_name))
                axes[col].axis('on')
            else:
                axes[row, col].imshow(wordcloud, interpolation="bilinear")
                axes[row, col].set_title('{}'.format(cluster_name))
                axes[row, col].axis('on')

        plt.savefig(output_file, bbox_inches='tight')
        print('Wordcloud saved to %s' % output_file)