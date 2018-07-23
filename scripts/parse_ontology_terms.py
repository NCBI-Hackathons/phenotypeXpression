import tqdm
import json
import rdflib
from rdflib import Namespace
from rdflib.namespace import RDF, RDFS, OWL

doid_path = '/Users/lwang/git/ontologies/doid.owl'
hpo_path = '/Users/lwang/git/ontologies/hp.owl'

doid_out_path = '/Users/lwang/git/phenotypeXpression/data/doid.json'
hpo_out_path = '/Users/lwang/git/phenotypeXpression/data/hp.json'

oboInOwl = Namespace('http://www.geneontology.org/formats/oboInOwl#')

def owl_to_json(owl_file, json_file):
    ent_dict = dict()
    g = rdflib.Graph()
    g.parse(owl_file)
    for s in tqdm.tqdm(g.subjects(RDF.type, OWL.Class)):
        labels = list(g.objects(s, RDFS.label))
        all_labels = [str(l) for l in labels]
        synonyms = list(g.objects(s, oboInOwl['hasExactSynonym']))
        all_labels += [str(s) for s in synonyms]
        if all_labels:
            ent_dict[s] = all_labels
    with open(json_file, 'wb') as f:
        f.write(json.dumps(ent_dict).encode('utf-8'))

owl_to_json(doid_path, doid_out_path)
owl_to_json(hpo_path, hpo_out_path)