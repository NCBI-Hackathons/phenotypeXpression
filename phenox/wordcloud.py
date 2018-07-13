from wordcloud import WordCloud
import matplotlib.pyplot as plt
import math

from phenox.paths import PhenoXPaths


class WordcloudPlotter:
    def __init__(self):
        self.paths = PhenoXPaths()

    def generate_wordclouds(self, clusters):
        """
        Generate plot with subplots of wordclouds
        :param clusters:
        :return:
        """
        num_clusters = len(clusters)
        num_columns = 3
        num_rows = int(math.ceil(num_clusters/num_columns))
        fig, axes = plt.subplots(num_rows, num_columns)
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

        # plot clusters
        for k, freqs in clusters.items():
            row = int(math.floor((k - 1) / num_columns))
            col = (k - 1) % num_columns
            wordcloud = WordCloud(background_color="white").generate_from_frequencies(freqs)
            if single_row:
                axes[col].imshow(wordcloud, interpolation="bilinear")
                axes[col].set_title('Group %i' % k)
                axes[col].axis('on')
            else:
                axes[row, col].imshow(wordcloud, interpolation="bilinear")
                axes[row, col].set_title('Group %i' % k)
                axes[row, col].axis('on')

        plt.show()