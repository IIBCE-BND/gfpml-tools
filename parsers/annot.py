EVIDENCE_BLACKLIST = ['IEA']


def add_to_map(dmap, key, value):
    if key in dmap and value not in dmap[key]:
        dmap[key].append(value)
    elif key not in dmap:
        dmap[key] = [value]

def parse_annot(path, alt_ids):
    gos_genes = {}
    genes_gos = {}

    with open(path, 'r') as f:
        for line in f.readlines():
            split = line.split()
            if len(split) > 6:
                gene_id = split[1 + (split[0] == 'UniProtKB')]
                go_id = split[3 + ('GO:' not in split[3])]
                evidence = split[5 + ('GO:' not in split[3])]
                
                # Replace alternative GO id with the original one.
                if go_id in alt_ids:
                    go_id = alt_ids[go_id]

                if evidence in EVIDENCE_BLACKLIST:
                    add_to_map(gos_genes, go_id, gene_id)
                    add_to_map(genes_gos, gene_id, go_id)

    return gos_genes, genes_gos