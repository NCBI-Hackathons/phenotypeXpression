import os
import json

from phenox.paths import PhenoXPaths


paths = PhenoXPaths()

hgnc_json_file = os.path.join(paths.data_dir, 'hgnc.json')
hgnc_syn_file = os.path.join(paths.data_dir, 'hgnc_synonyms.json')

with open(hgnc_json_file, 'r') as f:
    text = f.read()
    hgnc = json.loads(text)

num_entries = hgnc['response']['numFound']
data = hgnc['response']['docs']

hgnc_dict = dict()

for entry in data:
    try:
        symbol = entry['symbol']
        alias_symbol = entry['alias_symbol'] if 'alias_symbol' in entry else []
        prev_symbol = entry['prev_symbol'] if 'prev_symbol' in entry else []

        name = entry['name']
        prev_name = entry['prev_name'] if 'prev_name' in entry else []

        hgnc_dict[symbol] = {'names': [name] + prev_name,
                             'aliases': alias_symbol + prev_symbol}
    except KeyError:
        print('Entry has no symbol or name! ({})'.format(entry['hgnc_id']))
        continue

with open(hgnc_syn_file, 'w') as outf:
    json.dump(hgnc_dict, outf)






