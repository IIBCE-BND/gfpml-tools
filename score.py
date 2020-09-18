import click
import numpy as np
import pandas as pd
import os

from joblib import Parallel, delayed

from parsers import parse_gtf, parse_obo, parse_annot

from sklearn.model_selection import train_test_split
from sklearn.metrics import pairwise_distances

def score_function(num_genes_in_chromosome, positions):
#     return score for every gen position
    genes_positions_GO = np.array(positions).reshape(-1, 1)
    distances = pairwise_distances(genes_positions_GO, genes_positions_GO, metric='l1')
    distances = np.sort(distances, axis=1)
    mean = distances.mean(axis=0)
    # median = np.median(distances, axis=0)

    genes_positions = np.arange(num_genes_in_chromosome).reshape(-1, 1)
    distances = pairwise_distances(genes_positions, genes_positions_GO, metric='l1')
    norma = np.linalg.norm(np.sort(distances, axis=1) - mean, axis=1)
    # norma = np.linalg.norm(np.sort(distances, axis=1) - median, axis=1)

    score = num_genes_in_chromosome / (num_genes_in_chromosome + norma)
    
    return score


def calculate_score(seq_genes, seq_annots, save_path, th=100):
    if 'name' in seq_genes.columns: gene_identifier = 'name'
    elif 'id' in seq_genes.columns: gene_identifier = 'id'
    else: raise Exception

    chromosome = seq_genes.seqname.values[0]

    seg_genes = seq_genes.sort_values(['start', 'strand', 'size'], ascending=True)
    seq_genes['pos'] = range(len(seq_genes))
    seq_len = len(seq_genes)

    seq_score = {}
    for go_id, seq_annots_go in seq_annots.groupby('go_id'):
        seq_annots_go['pos'] = seq_annots_go['gene_id'].replace(seq_genes.set_index(gene_identifier)['pos'].to_dict())
        seq_score[go_id] = score_function(seq_len, seq_annots_go.pos.values)

    seq_score = pd.DataFrame(data=seq_score)
    seq_score.to_csv('{}/{}/seq_score.csv'.format(save_path, chromosome), sep='\t', index=False)


def calculate_score2(seq_genes, seq_annots, save_path, th=100):
    if 'name' in seq_genes.columns: gene_identifier = 'name'
    elif 'id' in seq_genes.columns: gene_identifier = 'id'
    else: raise Exception

    chromosome = seq_genes.seqname.values[0]

    seg_genes = seq_genes.sort_values(['start', 'strand', 'size'], ascending=True)
    seq_genes['pos'] = range(len(seq_genes))
    seq_len = len(seq_genes)

    seq_score = {}
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
        seq_score[go_id] = score_function(seq_len, pos_train)

    if len(seq_score) > 0:
        seq_score = pd.DataFrame(data=seq_score)
        seq_score.to_csv('{}/{}/seq_score.csv'.format(save_path, chromosome), sep='\t', index=False)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--genome', 'genome_path')
@click.option('--centromeres', 'centromeres_path')
@click.option('--ontology', 'ontology_path')
@click.option('--annotations', 'annotations_path')
@click.option('--save-path')
def generate(genome_path, centromeres_path, ontology_path, annotations_path, save_path):
    click.echo('Loading genome: {}'.format(genome_path))
    genome = parse_gtf(genome_path, centromeres_path)

    click.echo('Loading ontology: {}'.format(ontology_path))
    gos, go_alt_ids = parse_obo(ontology_path)

    click.echo('Loading annotations: {}'.format(annotations_path))
    annots = parse_annot(annotations_path, go_alt_ids)

    # 'name' is not in genome.columns for scer, instead we must use 'id'
    if 'name' in genome.columns: gene_identifier = 'name'
    elif 'id' in genome.columns: gene_identifier = 'id'
    else: raise Exception

    def select_annots(seq_genes):
        #  return annotations of genes in genes
        return annots[annots['gene_id'].isin(seq_genes[gene_identifier].values)]

    click.echo('Calculating genes score input for {} genes'.format(len(genome)))
    Parallel(n_jobs=-1, verbose=10)(
        delayed(calculate_score)(seq_genes, select_annots(seq_genes), save_path)
        for _, seq_genes in genome.groupby('seqname')
    )

    # Parallel(n_jobs=-1, verbose=10)(
    #     delayed(calculate_score2)(seq_genes, select_annots(seq_genes), save_path)
    #     for _, seq_genes in genome.groupby('seqname')
    # )

if __name__ == '__main__':
    cli()
