import click
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np

from scipy.stats import binom_test


def calculate_p_value(N, B, n, b):
    """
    Apply Binomial-test to calculate the p-value.
    """
    x = b
    n = n
    p = float(B) / N
    return binom_test(x, n, p, alternative='greater')


def calculate_term_enrichment(N, B, n, b):
    """
    Calculate the GO term enrichment in a target set.
    """
    return (float(b) / float(n)) / (float(B) / float(N))


def calculate_windows_enrichment(terms, seq_len, target_set, p_value_treshold,
                                 start, end):
    """
    Given a background set and target set it calculates the probability that
    the target set have b or more genes associated with a GO term.
    Those values that are not statistically significant are discarded.
    """
    N = seq_len
    n = len(target_set)
    window_enritchment = []
    for term in terms:
        B = len(term['genes'])
        b = len([
            g['gene_id'] for g in term['genes'] if g['gene_id'] in target_set])
        p_value = calculate_p_value(N, B, n, b)
        enrichment = 0
        if p_value < p_value_treshold:
            enrichment = calculate_term_enrichment(N, B, n, b)
        window_enritchment.append((term['_id'], enrichment, p_value))
    print('start: {}, end: {}'.format(start, end))
    return window_enritchment


def calculate_enrichment(terms, seq_len, chromosome_genes, window_size,
                         p_value_treshold):
    """
    For each gene on the sequence calculate the enrichment for every GO term
    of the window cenetered on the gene.
    """
    for idx, gene in enumerate(chromosome_genes):
        start = idx - min(idx, window_size)
        end = idx + min(window_size, len(chromosome_genes) - idx)
        target = [g['_id'] for g in chromosome_genes[start:end]]
        result = calculate_windows_enrichment(
            terms, seq_len, target, p_value_treshold, start, end)
        yield result


def calculate_enrichment_parallel(terms, seq_len, chromosome_genes,
                                  window_size, p_value_treshold):
    """
    For each gene on the sequence calculate the enrichment for every GO term
    of the window cenetered on the gene.
    """
    pool = multiprocessing.Pool()
    results = []
    for idx, gene in enumerate(chromosome_genes):
        start = idx - min(idx, window_size)
        end = idx + min(window_size, len(chromosome_genes) - idx)
        target = [g['_id'] for g in chromosome_genes[start:end]]
        result = pool.apply_async(calculate_windows_enrichment,
                                  args=(terms, seq_len, target,
                                        p_value_treshold, start, end))
        results.append(result)
    pool.close()
    pool.join()
    return [r.get() for r in results]


def plot_window_enrichment(term, background_len, windows_size, values):
    """
    Plots the genoma enrichment for a specific GO term.
    """
    fig, ax = plt.subplots()
    ind = np.arange(0, background_len)
    ax.bar(ind, values, 1, color='blue')
    # ax.set_xticks(np.arange(0, background_len, background_len / 10))
    ax  .set_title('{} - window: {}'.format(term, windows_size))
    ax.set_xlabel('Sequence')
    ax.set_ylabel('Enrichment')
    plt.savefig('plots/{}-{}.png'.format(term, windows_size))
