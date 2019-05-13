import click
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd

from itertools import chain

from joblib import Parallel, delayed
from scipy.stats import binom_test
from tqdm import tqdm

from parsers import parse_gtf, parse_obo, parse_annot


def calculate_p_value(N, B, n, b):
    """
    Apply Binomial-test to calculate the p-value.
    """
    x = b
    n = n # ???
    p = float(B) / N
    return binom_test(x, n, p, alternative='greater')


def calculate_term_enrichment(N, B, n, b):
    """
    Calculate the GO term enrichment in a target set.
    """
    return (float(b) / float(n)) / (float(B) / float(N))


def calculate_seq_lea(seq_genes, window_size, gos_genes, genes_gos, target):
    seq_len = len(seq_genes)
    if 'name' in seq_genes.columns:
        seq_gene_names = seq_genes.sort_values(['start', 'size'], ascending=True).name.values # correct order
    elif 'id' in seq_genes.columns:
        seq_gene_names = seq_genes.sort_values(['start', 'size'], ascending=True).id.values # correct order
    else:
        raise Exception

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


@click.group()
def cli():
    pass


@cli.command()
@click.option('--genome', 'genome_path')
@click.option('--ontology', 'ontology_path')
@click.option('--annotations', 'annotations_path')
@click.option('--window', 'window_size', default=10)
@click.option('--target', default=None)
@click.option('--save-path')
def generate(genome_path, ontology_path, annotations_path, window_size, target, save_path):
    click.echo('Loading genome: {}'.format(genome_path))
    genome = parse_gtf(genome_path)

    click.echo('Loading ontology: {}'.format(ontology_path))
    gos, go_alt_ids = parse_obo(ontology_path)

    click.echo('Loading annotations: {}'.format(annotations_path))
    gos_genes, genes_gos = parse_annot(annotations_path, go_alt_ids)

    # Ingore genes without annotations.
    # 'name' is not in genome.columns for scer
    if 'name' in genome.columns:
        genome = genome[genome.name.isin(genes_gos.keys())]
    elif 'id' in genome.columns:
        genome = genome[genome.id.isin(genes_gos.keys())]
    else:
        raise Exception

    click.echo('Calculating genes local enrichment for {} genes'.format(len(genome)))
    results = Parallel(n_jobs=-1, verbose=10)(
        delayed(calculate_seq_lea)(seq_genes, window_size, gos_genes, genes_gos, target)
        for _, seq_genes in genome.groupby('seqname')
    )

    if save_path:
        lea = pd.DataFrame(
            chain(*results),
            columns=['gene', 'go', 'enrichment', 'p_value']
        )
        lea.to_csv(save_path, index=False)


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
