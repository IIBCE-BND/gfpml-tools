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
        ans[min(n, max(w - i, 0)):max(0, min(n + w - i, n)),i] = gene_pos[min(n, max(i - w, 0)):max(0, min(n, n + i - w))]
    return ans


def calculate_enrichment(gene_pos, window_size):
    n = len(gene_pos)
    w = window_size
    if w < n - 1:
        m = min(2 * w + 1, n)
        h = m - (w + 1)
        window_sizes = np.repeat(m, n)
        window_sizes[:h] = np.arange(w + 1, m)
        window_sizes[-h:] = np.arange(m - 1, w, -1)
    else:
        window_sizes = np.repeat(n, n)

    lea_array = (reshape_lea(gene_pos, window_size).sum(axis=1) / window_sizes) / (gene_pos.mean())
    return lea_array


def calculate_seq_lea(seq_genes, window_size, gos_genes, genes_gos, target):
    seq_len = len(seq_genes)
    if 'name' in seq_genes.columns: gene_identifier = 'name'
    elif 'id' in seq_genes.columns: gene_identifier = 'id'
    else: raise Exception

    seq_gene_names = seq_genes.sort_values(['start', 'strand', 'size'], ascending=True)[gene_identifier].values # correct order

    seq_gos_genes = {
        go_id: set(go_genes) & set(seq_gene_names)
        for go_id, go_genes in gos_genes.items()
    }

    seq_lea = []
    for pos, name in enumerate(seq_gene_names):
        start = pos - min(pos, window_size)
        end = pos + min(window_size, seq_len - 1 - pos) # porque pos < seq_len

        window_genes = [g for g in seq_gene_names[start:end]]
        window_gos = set(chain(*[genes_gos[g] for g in window_genes])) # maybe some of this GO_id are repeated

        N = seq_len
        n = len(window_genes)
        e = []
        for go_id in window_gos:
            go_genes = set(seq_gos_genes[go_id])
            if go_id == target and name in go_genes:
                go_genes.remove(name)
            B = len(go_genes)
            b = len(go_genes & set(window_genes))
            if B:
                p_value = calculate_p_value(N, B, n, b)
                go_e = calculate_term_enrichment(N, B, n, b)
            else:
                p_value = 0
                go_e = 0
            seq_lea.append((name, go_id, go_e, p_value))

    return seq_lea


def calculate_seq_lea2(seq_genes, seq_annots, ontology, window_size, save_path, th=100):
    if 'name' in seq_genes.columns: gene_identifier = 'name'
    elif 'id' in seq_genes.columns: gene_identifier = 'id'
    else: raise Exception

    chromosome = seq_genes.seqname.values[0]

    seq_genes = seq_genes.sort_values(['start', 'strand', 'size'], ascending=True)
    seq_genes['pos'] = range(len(seq_genes))
    seq_len = len(seq_genes)

    seq_lea_ontology = {}
    for go_id, seq_annots_go in seq_annots.groupby('go_id'):
        if len(seq_annots_go) < th:
            continue

        seq_annots_go['pos'] = seq_annots_go['gene_id'].replace(seq_genes.set_index(gene_identifier)['pos'].to_dict())

        train_size = 0.8
        X_train, X_test, pos_train, pos_test = train_test_split(seq_annots_go.gene_id.values, seq_annots_go.pos.values, train_size=train_size)

        train_df = pd.DataFrame(data={'pos':pos_train, 'gene_id':X_train})
        test_df = pd.DataFrame(data={'pos':pos_test, 'gene_id':X_test})

        if not os.path.exists('{}'.format(save_path)):
            os.mkdir('{}'.format(save_path))
        if not os.path.exists('{}/{}'.format(save_path, chromosome)):
            os.mkdir('{}/{}'.format(save_path, chromosome))

        train_df.to_csv('{}/{}/{}_train.csv'.format(save_path, chromosome, go_id), sep='\t', index=False)
        test_df.to_csv('{}/{}/{}_test.csv'.format(save_path, chromosome, go_id), sep='\t', index=False)

        mask_train = np.isin(range(seq_len), pos_train)
        seq_lea_ontology[go_id] = calculate_enrichment(mask_train, window_size)

    if len(seq_lea_ontology) > 0:
        seq_lea_ontology = pd.DataFrame(data=seq_lea_ontology)
        seq_lea_ontology.to_csv('{}/{}/seq_lea_{}.csv'.format(save_path, chromosome, ontology), sep='\t', index=False)


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
    gos, ontology_gos, go_alt_ids, ontology_graphs = parse_obo(ontology_path)

    click.echo('Loading annotations: {}'.format(annotations_path))
    annots = parse_annot(annotations_path, go_alt_ids)

    # 'name' is not in genome.columns for scer, instead we must use 'id'
    if 'name' in genome.columns: gene_identifier = 'name'
    elif 'id' in genome.columns: gene_identifier = 'id'
    else: raise Exception

    def select_annots(annots, seq_genes):
        #  return annotations of genes in genes
        return annots[annots['gene_id'].isin(seq_genes[gene_identifier].values)]

    click.echo('Calculating genes local enrichment for {} genes'.format(len(genome)))
    for ontology in ontology_gos:
        annots_ontology = annots[annots['go_id'].isin(ontology_gos[ontology])]
        ontology_graph = ontology_graphs[ontology]
        expanded_annots = expand_annots(annots_ontology, ontology_graph)

        # Parallel(n_jobs=-1, verbose=10)(
        #     (seq_genes, window_size, gos_genes, genes_gos, target)
        #     delayed(calculate_seq_lea)(seq_genes, select_annots(expanded_annots, seq_genes), ontology, window_size, save_path)
        #     for _, seq_genes in genome.groupby('seqname')
        # )

        Parallel(n_jobs=-1, verbose=10)(
            delayed(calculate_seq_lea2)(seq_genes, select_annots(expanded_annots, seq_genes), ontology, window_size, save_path)
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
