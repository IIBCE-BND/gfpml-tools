import re

ONTOLOGIES = ['molecular_function', 'biological_process', 'cellular_component']
RELATIONSHIP_TYPES = ['part_of', 'regulates', 'positively_regulates', 'negatively_regulates']
ANNOTATION_RE = '(?m)^id: GO:.*(?:\r?\n(?!\[(?:Typedef|Term)\]).*)*'


def parse_obo(path):
    gos = {}  # Map between Go terms and their ontologies, it contains alternative GO ids.
    alt_ids = {}  # Dictionary to keep track of alternative ontologies ids.
    
    with open(path, 'r') as f:
        content = f.read()

        for row in re.findall(ANNOTATION_RE, content):
            if 'is_obsolete: true' in row:
                continue

            cols = row.split()
            ontology = cols[cols.index('namespace:') + 1]
            go_id = cols[1]
            gos[go_id] = ontology

            go_alt_ids = [cols[idx + 1] for idx, x in enumerate(cols) if x == 'alt_id:']
            for alt_id in go_alt_ids:
                alt_ids[alt_id] = go_id

    return gos, alt_ids