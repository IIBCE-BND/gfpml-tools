import pandas as pd
import numpy as np

import re

CHROMOSOMES = {'Caenorhabditis_elegans': ['I','II','III','IV','V','X'],
               'Drosophila_melanogaster': ['2L','2R','3L','3R','4','X','Y'],
               'Homo_sapiens': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'],
               'Mus_musculus': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y'],
               'Saccharomyces_cerevisiae': ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
               }

GTF_HEADER = [
    'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'
]
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA = re.compile(r'\s*,\s*')
R_KEYVALUE = re.compile(r'(\s+|\s*=\s*)')

FEATURE_GENE_VALUES = ['gene']

# List of alternatives columns names for each field.
COLUMNS = {
    'id': ['gene_id', 'id'],
    'name': ['gene_name', 'name'],
    'biotype': ['gene_biotype', 'biotype'],
}

ATTRIBUTES_FIELDS = [
    'gene_biotype', 'biotype', 'gene_name', 'name', 'gene_id', 'id',
    'gene_type'
]

def parse_gtf(path, centromere_path):
    df = pd.read_csv(
        path,
        delimiter='\t',
        names=GTF_HEADER,
        usecols=GTF_HEADER,
        dtype={
            'seqname':str, 'source':str, 'feature':str, 'start':int, 'end':int, 'score':str,
            'strand':str, 'frame':str, 'attributes':str
        },
        engine='c',
        low_memory=True,
        compression='infer',
        comment='#'
    )

    organism_name = re.findall(r'\/(\w+?)\.', path)[0]
    df = df[df.seqname.isin(CHROMOSOMES[organism_name])]
    df = centromere_split(df, organism_name, centromere_path)

    # Exclude rows that are not genes.
    df = df[df.feature.isin(FEATURE_GENE_VALUES)]

    # this is to compute gene position in lea implementation
    df['strand'] = df['strand'].replace({'+':1, '-':-1})
    df['size'] = np.abs(df.end - df.start)

    if df.empty:
        raise Exception

    df = pd.concat([
        df, df['attributes'].apply(attr_to_series)
    ], axis=1)

    df.rename(
        index=str,
        columns={
            'gene_id': 'id', 'gene_name': 'name', 'gene_biotype': 'biotype',
            'gene_type': 'biotype'
        },
        inplace=True
    )

    return df


def centromere_split(genome, organism_name, centromere_path):
    if organism_name == 'Homo_sapiens':
        hg_centromeres_regions = pd.read_csv(centromere_path, delimiter='\t')
        for chromosome, chromosome_genes in genome.groupby('seqname'):
            # get chromsome region
            centromere = hg_centromeres_regions[hg_centromeres_regions.chromosome == chromosome]
            if centromere.empty:
                continue

            # filter chromosome regions
            region_1 = chromosome_genes[chromosome_genes.start < centromere.iloc[0].start]
            region_2 = chromosome_genes[chromosome_genes.start > centromere.iloc[0].end]

            # add sufixb
            region_q = region_1 if len(region_1) > len(region_2) else region_2
            region_p = region_2 if len(region_1) > len(region_2) else region_1
            genome.loc[region_q.index, 'seqname'] = '{}q'.format(chromosome)
            genome.loc[region_p.index, 'seqname'] = '{}p'.format(chromosome)


    if organism_name == 'Saccharomyces_cerevisiae':
        scer_centromeres_regions = pd.read_csv(centromere_path, delimiter='\t')
        # groupby genes by chromomsome
        for chromosome, chromosome_genes in genome.groupby('seqname'):

            # get chromsome region
            centromere = scer_centromeres_regions[scer_centromeres_regions.chromosome == chromosome]
            if centromere.empty:
                continue

            # filter chromosome region
            region_1 = chromosome_genes[chromosome_genes.start < centromere.iloc[0].start]
            region_2 = chromosome_genes[chromosome_genes.start > centromere.iloc[0].end]

            # add sufix
            genome.loc[region_1.index, 'seqname'] = '{}L'.format(chromosome)
            genome.loc[region_2.index, 'seqname'] = '{}R'.format(chromosome)

    return genome


def attr_to_series(attr):
    keys, values = zip(
        *[(dct['label'], dct['value']) for dct in parse_attributes(attr)]
    )
    return pd.Series(values, index=keys)


def parse_attributes(attr):
    infos = [x for x in attr.split(';') if x.strip()]

    result = []
    for i, info in enumerate(infos, 1):
        info = info.strip()

        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
            key = key.lower()

            if key not in ATTRIBUTES_FIELDS:
                continue

        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result.append({'label': key.lower(), 'value': _get_value(value)})

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)

    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value
