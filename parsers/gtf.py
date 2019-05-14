import pandas as pd
import numpy as np

import re


GTF_HEADER = [
    'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame'
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


def parse_gtf(path):
    df = pd.read_csv(
        path,
        delimiter='\t',
        names=[
            'seqname', 'source', 'feature', 'start', 'end', 'score',
            'strand', 'frame', 'attributes'
        ],
        usecols=[
            'seqname', 'source', 'feature', 'start', 'end', 'score',
            'strand', 'frame', 'attributes'
        ],
        # sometime a chromosome is loaded as 1(int) and sometime as '1'(str), ensure than start and end are integer and not strings
        dtype={
            'seqname':str, 'source':str, 'feature':str, 'start':int, 'end':int, 'score':str,
            'strand':str, 'frame':str, 'attributes':str
        },
        engine='c',
        low_memory=True,
        compression='infer',
        comment='#'
    )

    ###############################################################
    ## TODO: ignore mito, MT, GL456233.1 chromosomes
    ###############################################################

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
