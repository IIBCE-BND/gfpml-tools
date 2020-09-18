import networkx as nx
import re

ONTOLOGIES = ['molecular_function', 'biological_process', 'cellular_component']
RELATIONSHIP_TYPES = ['part_of', 'regulates', 'positively_regulates', 'negatively_regulates']
ANNOTATION_RE = '(?m)^id: GO:.*(?:\r?\n(?!\[(?:Typedef|Term)\]).*)*'

ontology_graphs = {ontology: nx.DiGraph() for ontology in ONTOLOGIES}
ontology_gos = {ontology:[] for ontology in ONTOLOGIES}

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
            ontology_gos[ontology].append(go_id)
            ontology_graphs[ontology].add_node(go_id)

            go_alt_ids = [cols[idx + 1] for idx, x in enumerate(cols) if x == 'alt_id:']
            for alt_id in go_alt_ids:
                alt_ids[alt_id] = go_id

            is_a = [cols[idx + 1] for idx, x in enumerate(cols) if x == 'is_a:']
            for related_GO in is_a:
                ontology_graphs[ontology].add_edge(go_id, related_GO)

    return gos, ontology_gos, alt_ids, ontology_graphs
