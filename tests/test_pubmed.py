import unittest
from collections import Counter
from phenox.pubmed import Pubmed

# test: pmid list for psoriasis
psoriasis = ['24646743',
             '22479649',
             '25129481',
             '25129481',
             '23633458']


class TestPubMed(unittest.TestCase):
    def test_pubmed(self):
        pubmed = Pubmed("test@test.xyz", psoriasis)

        for pmid in pubmed.pmid_list:
            pubmed.fetch_abstract(pmid)
            pubmed.extract_DNER(pmid)

        pubmed.freq_list = Counter(pubmed.total_dner)

        assert(len(pubmed.pmid_abstracts) <= len(psoriasis))
        assert(set(pubmed.pmid_abstracts.keys()).difference(set(psoriasis)) == set([]))
        assert(len(pubmed.freq_list) == 12)