import click
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd
import os

from itertools import chain

from joblib import Parallel, delayed
from scipy.stats import binom_test
from tqdm import tqdm

from parsers import parse_gtf, parse_obo, parse_annot
from parsers.annot import expand_annots

from sklearn.model_selection import train_test_split


def calculate_p_value(N, B, n, b):
    """
    Apply Binomial-test to calculate the p-value.
    """
    x = b
    n = n # ???
    p = float(B) / N
    return binom_test(x, n, p, alternative='greater')


def reshape_lea(gene_pos, window_size):
    n = len(gene_pos)
    w = window_size
    ans = np.zeros((n, 2 * w + 1))
    for i in range(2 * w + 1):
        ans[max(w - i, 0):min(n + w - i, n),i] = gene_pos[max(i - w, 0):min(n, n + i - w)]

    return ans

def calculate_enrichment(gene_pos, window_size):
    w = window_size
    windows_size = np.repeat(2 * window_size + 1, len(gene_pos))
    windows_size[:window_size] = np.arange(window_size + 1, 2 * window_size + 1)
    windows_size[-window_size:] = np.arange(2 * window_size, window_size, -1)

    lea_array = (reshape_lea(gene_pos, window_size).sum(axis=1) / windows_size) / (gene_pos.mean())

    return lea_array

def calculate_seq_lea2(seq_genes, seq_annots, window_size, save_path, th=100):
    if 'name' in seq_genes.columns: gene_identifier = 'name'
    elif 'id' in seq_genes.columns: gene_identifier = 'id'
    else: raise Exception

    chromosome = seq_genes.seqname.values[0]

    seq_genes = seq_genes.sort_values(['start', 'strand', 'size'], ascending=True) # correct order
    seq_genes['pos'] = range(len(seq_genes))
    seq_len = len(seq_genes)

    seq_lea = {}
    for go_id, seq_annots_go in seq_annots.groupby('go_id'):
        if len(seq_annots_go) < th:
            continue

        seq_annots_go['pos'] = seq_annots_go['gene_id'].replace(seq_genes.set_index(gene_identifier)['pos'].to_dict())

        train_size = 0.8
        X_train, X_test, pos_train, pos_test = train_test_split(seq_annots_go.gene_id.values, seq_annots_go.pos.values, train_size=train_size)

        train_data = {'pos':pos_train, 'gene_id':X_train}
        train_df = pd.DataFrame(data=train_data)
        test_data = {'pos':pos_test, 'gene_id':X_test}
        test_df = pd.DataFrame(data=test_data)

        if not os.path.exists('{}'.format(save_path)):
            os.mkdir('{}'.format(save_path))
        if not os.path.exists('{}/{}'.format(save_path, chromosome)):
            os.mkdir('{}/{}'.format(save_path, chromosome))

        train_df.to_csv('{}/{}/{}_train.csv'.format(save_path, chromosome, go_id), sep='\t', index=False)
        test_df.to_csv('{}/{}/{}_test.csv'.format(save_path, chromosome, go_id), sep='\t', index=False)

        mask_train = np.isin(range(seq_len), pos_train)
        seq_lea[go_id] = calculate_enrichment(mask_train, window_size)

    if len(seq_lea) > 0:
        seq_lea = pd.DataFrame(data=seq_lea)
        seq_lea.to_csv('{}/{}/seq_lea.csv'.format(save_path, chromosome), sep='\t', index=False)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--genome', 'genome_path')
@click.option('--centromeres', 'centromeres_path')
@click.option('--ontology', 'ontology_path')
@click.option('--annotations', 'annotations_path')
@click.option('--window', 'window_size', default=10)
@click.option('--save-path')
def generate(genome_path, centromeres_path, ontology_path, annotations_path, window_size, save_path):
    click.echo('Loading genome: {}'.format(genome_path))
    genome = parse_gtf(genome_path, centromeres_path)

    click.echo('Loading ontology: {}'.format(ontology_path))
    gos, go_alt_ids, ontology_graphs = parse_obo(ontology_path)

    click.echo('Loading annotations: {}'.format(annotations_path))
    annots = parse_annot(annotations_path, go_alt_ids)
    print('ahora agreagr ontology')
    annots['ontology'] = annots['go_id'].replace(gos)
    print('termina agreagr ontology')

    # 'name' is not in genome.columns for scer, instead we must use 'id'
    if 'name' in genome.columns: gene_identifier = 'name'
    elif 'id' in genome.columns: gene_identifier = 'id'
    else: raise Exception

    def select_annots(annots, seq_genes):
        #  return annotations of genes in genes
        return annots[annots['gene_id'].isin(seq_genes[gene_identifier].values)]

    click.echo('Calculating genes local enrichment for {} genes'.format(len(genome)))
    for ontology, annots_ontology in annots.groupby('ontology'):
        print(ontology)
        ontology_graph = ontology_graphs[ontology]
        print('expanding')
        expanded_annots = expand_annots(annots_ontology, ontology, ontology_graph)
        print('finished expanding')

        Parallel(n_jobs=-1, verbose=10)(
            delayed(calculate_seq_lea2)(seq_genes, select_annots(expanded_annots, seq_genes), window_size, save_path)
            for _, seq_genes in genome.groupby('seqname')
        )


@cli.command()
@click.argument('path')
@click.option('--threshold', default=.005)
@click.option('--save-path')
def transform(path, threshold, save_path):
    lea = pd.read_csv(path)
    go_idx_map = {go: idx for idx, go in enumerate(lea.go.sort_values().unique())}

    def generate_input(gene_lea):
        x = np.zeros(len(go_idx_map))
        for _, row in gene_lea[gene_lea.p_value > threshold].iterrows():
            x[go_idx_map[row.go]] = row.enrichment
        return x

    X = lea.groupby('gene').apply(generate_input)

    if save_path:
        X.to_pickle(save_path)


if __name__ == '__main__':
    cli()
