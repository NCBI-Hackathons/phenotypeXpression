import unittest

from phenox.mesh_lookup import MeshSearcher


class TestMeshSearcher(unittest.TestCase):

    def test_search_correct_term(self):
        mesh = MeshSearcher()
        match = mesh.lookup('psoriasis')

        assert(match['name'] == 'Psoriasis')
        assert('C17.800.859.675' in match['ids'])
