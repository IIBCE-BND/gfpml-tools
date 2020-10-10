import pandas as pd

EVIDENCE_BLACKLIST = ['IEA']
COLUMNS_NAMES = ['DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID', 'DB:Reference (|DB:Reference)', 'Evidence Code', 'With (or) From', 'Aspect', 'DB Object Name', 'DB Object Synonym (|Synonym)', 'DB Object Type', 'Taxon(|taxon)', 'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID']

get_name = lambda attr: attr.split('|')[0]

def parse_annot(path, alt_ids):
    '''
    Given the path to the annotations file and alt_ids a dictionary whose keys are obsoletes GO terms and their values
    are the replacement terms returns a dataframe with all annotations
    '''
    df = pd.read_csv(
        path,
        delimiter='\t',
        names=COLUMNS_NAMES,
        usecols=['DB', 'DB Object Symbol', 'GO ID', 'Evidence Code', 'DB Object Synonym (|Synonym)'],
        engine='c',
        low_memory=True,
        compression='infer',
        comment='!'
    )

    df.dropna(subset=['DB Object Symbol'], inplace=True) # when load .tar.gz first row is nan
    df = df[~df['Evidence Code'].isin(EVIDENCE_BLACKLIST)]
    df.drop_duplicates(subset=['DB Object Symbol', 'GO ID'], inplace=True)
    df['GO ID'].replace(alt_ids, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # recover gene_id used in .gtf file
    if df.at[0, 'DB'] == 'SGD':
        df['DB Object Symbol'] = df['DB Object Synonym (|Synonym)'].apply(get_name)

    df = df[['DB Object Symbol', 'GO ID']]
    df.rename(columns={'DB Object Symbol': 'gene_id', 'GO ID': 'go_id'}, inplace=True)

    return df

def get_ancestors(go_id, ontology_graph):
    # return the set of ancestors of GO terms of go_id in ontology_graph
    ans = set([go_id])
    aux = set([go_id])
    while len(aux) > 0:
        for node in aux:
            ans = ans.union(set(ontology_graph.successors(node)))
            aux = aux.union(set(ontology_graph.successors(node)))
            aux = aux - {node}
    return ans

def expand_annots(annots, ontology_graph):
    # Given a set of annotations and the ontology graph return the hierarchical annotations
    expanded_annots = []
    if 'pos' in annots.columns:
        for go_id, annots_go in annots.groupby('go_id'):
            ancestors = get_ancestors(go_id, ontology_graph)
            for _, row in annots_go.iterrows():
                gene_id = row['gene_id']
                pos = row['pos']
                seqname = row['seqname']
                df = pd.DataFrame({'go_id': list(ancestors), 'gene_id': gene_id, 'pos': pos, 'seqname': seqname})
                expanded_annots.append(df)
    else:
        for go_id, annots_go in annots.groupby('go_id'):
            ancestors = get_ancestors(go_id, ontology_graph)
            for gene_id in annots_go.gene_id.values:
                df = pd.DataFrame({'go_id': list(ancestors), 'gene_id': gene_id})
                expanded_annots.append(df)
    expanded_annots = pd.concat(expanded_annots)
    expanded_annots = expanded_annots.drop_duplicates()
    return expanded_annots
