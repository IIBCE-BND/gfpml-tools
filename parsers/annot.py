import pandas as pd

EVIDENCE_BLACKLIST = ['IEA']
COLUMNS_NAMES = ['DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID', 'DB:Reference (|DB:Reference)', 'Evidence Code', 'With (or) From', 'Aspect', 'DB Object Name', 'DB Object Synonym (|Synonym)', 'DB Object Type', 'Taxon(|taxon)', 'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID']

get_name = lambda attr: attr.split('|')[0]

def parse_annot(path, alt_ids):
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
