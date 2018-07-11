import os
import sys
import json
import glob
import difflib

from phenox.paths import PhenoXPaths


class MeshSearcher:
    def __init__(self):
        """
        Initialize MeSH search tool
        """
        paths = PhenoXPaths()
        mesh_json_path = os.path.join(paths.data_dir, 'mesh.json')
        self.mesh = dict()

        if not os.path.exists(mesh_json_path):
            mesh_bin_file = glob.glob(os.path.join(paths.data_dir, '*.bin'))
            if mesh_bin_file:
                self._parse_mesh_bin(mesh_bin_file[0], mesh_json_path)

        self.mesh = json.load(open(mesh_json_path, 'r'))

    def _parse_mesh_bin(self, bin_file, json_file):
        """
        Parse MeSH tree from bin file
        :param bin_file:
        :param json_file:
        :return:
        """

        def _chunks(filename, start):
            """
            Split file into chunks
            :param filename:
            :param start:
            :return:
            """
            with open(filename, 'r') as f:
                buffer = []
                for line in f:
                    if line.startswith(start):
                        if buffer:
                            yield buffer
                            buffer = []
                    else:
                        buffer.append(line.strip())

        mesh = dict()

        # iterate through each chunk and parse MeSH record
        for c in _chunks(bin_file, '*NEWRECORD'):
            name = None
            ids = []
            aliases = []

            for l in c:
                # parse mesh id
                if l.startswith('MN = '):
                    id = l.split('=')[1].strip()
                    if id.startswith('C'):
                        ids.append(id)
                # parse main heading
                elif l.startswith('MH = '):
                    name = l.split('=')[1].strip()
                # parse aliases
                elif l.startswith('ENTRY = '):
                    aliases.append(l.split('=')[1].strip().split('|')[0].lower())

            # if disease ID, create a record
            if ids:
                record = {'ids': ids,
                          'name': name,
                          'aliases': aliases,
                          'parents': [],
                          'children': []}
                if name and name.lower() not in mesh:
                    mesh[name.lower()] = record
                else:
                    sys.stdout.write('Duplicate name! %s\n' % name)

        # create MeSH id to MeSH MH dictionary
        id_to_name_dict = dict()

        for name, data in mesh.items():
            for id in data['ids']:
                id_to_name_dict[id] = name

        # get parent child relationships and add back into MeSH tree
        par_map = list()

        for name, data in mesh.items():
            for id in data['ids']:
                par_id = '.'.join(id.split('.')[:-1])
                if par_id and name:
                    par_map.append((name, id_to_name_dict[par_id]))

        for ent, par in par_map:
            mesh[ent]['parents'].append(par)
            mesh[par]['children'].append(ent)

        json.dump(mesh, open(json_file, 'w'))

        return

    def lookup(self, query_text):
        """
        Find closest MeSH term
        :param query_text:
        :return:
        """
        query = query_text.lower()
        if query in self.mesh.keys():
            return self.mesh[query]
        else:
            closest = difflib.get_close_matches(query, self.mesh.keys())
            print('Did you mean?')
            for ind, match in enumerate(closest):
                print('%i) %s' % (ind + 1, self.mesh[match]['name']))
            selection = input()
            try:
                if int(selection) <= len(closest):
                    return self.mesh[closest[0]]
                else:
                    sys.stdout.write('Not a known selection, exiting...\n')
                    sys.exit(0)
            except Exception:
                sys.stdout.write('Unknown exception, exiting...\n')
                sys.exit(0)












