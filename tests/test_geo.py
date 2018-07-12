import os
import json
import unittest
import Bio.Entrez as Entrez
from phenox.paths import PhenoXPaths
from phenox.get_docsum import GEOQuery


paths = PhenoXPaths()
test_data_path = os.path.join(paths.test_dir, 'data', 'test_docsums.json')

class TestGEO(unittest.TestCase):
    def test_geo_data_fetcher(self):
        geo = GEOQuery("test@test.xyz")
        query_results = geo.get_ncbi_docsum("Muscular Dystrophies", "geoprofiles")
        assert(len(query_results) == 1)
        assert(len(query_results[0]['DocumentSummarySet']['DocumentSummary']) >= 844)

    def test_gdsdict_from_profile(self):
        geo = GEOQuery("test@test.xyz")

        data = []
        with open(test_data_path, 'r') as f:
            for l in f:
                data.append(json.loads(l))

        # output gdsdict and overall gene frequency
        gds_dict, gene_freq = geo.gdsdict_from_profile(data)

        assert(len(gds_dict) == 15)
        assert(len(gene_freq) == 2413)

        gds_list = list(gds_dict.keys())
        assert('5420' in gds_list)

        # elink gds to pid
        pid_list = Entrez.read(
            Entrez.elink(
                id=gds_list, db='pubmed', dbfrom='gds', linkname='gds_pubmed'
            )
        )

        # extract pubmed IDs
        pids = [el['LinkSetDb'][0]['Link'][0]['Id'] for el in pid_list]

        assert(len(pids) == 15)
        assert(len(pids) == len(gds_dict))
        assert('24646743' in pids)

    def test_genedict_from_profile(self):
        geo = GEOQuery("test@test.xyz")

        data = []
        with open(test_data_path, 'r') as f:
            for l in f:
                data.append(json.loads(l))

        # output gene ID => name table
        gene_dict = geo.genedict_from_profile(data)

        assert(len(gene_dict) == 2404)
        assert('7057' in gene_dict)
        assert(gene_dict['7057'] == 'THBS1')






