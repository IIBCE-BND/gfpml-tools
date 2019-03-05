import aiohttp
import asyncio
import click
import json
import pandas as pd
import time

from parsers import parse_gtf
from tqdm import tqdm

tqdm().pandas()

BASE_URL = 'https://rest.ensembl.org'

MAX_CONCURRENT_TASKS = 15

HOMOLOGY_KEYS = ['type', 'method_link_type']
HOMOLOGY_TARGET_KEYS = ['id', 'species']


@click.group()
def cli():
    pass


async def fetch(session, url):
    async with session.get(url) as response:
        return await response.text()


def parse_homology(h):
    data = {}
    for k in HOMOLOGY_TARGET_KEYS:
        data[k] = h['target'][k]
    for k in HOMOLOGY_KEYS:
        data[k] = h[k]
    return data


async def get_homology(gene_id):
    url = '{}/homology/id/{}'.format(BASE_URL, gene_id)
    url += '?content-type=application/json'

    async with aiohttp.ClientSession() as session:
        res = await fetch(session, url)

    try:
        data = [
            {'source_id': gene_id, **parse_homology(h)}
            for h in json.loads(res)['data'][0]['homologies']
        ]
    except Exception as e:
        click.echo('Error extracting homologies: {} {} {}'.format(e, url, res))
        return []

    return data


async def get_homologies(loop, gene_ids, fast_download):
    homologies = []
    current_tasks = set()
    ts_queue = []
    for gene_id in tqdm(gene_ids):
        if fast_download:
            current_time = time.time()
            if len(ts_queue) == 15:
                ts1, *ts_queue = ts_queue
                elapsed_time = current_time - ts1
                if elapsed_time < 1:
                    done, current_tasks = await asyncio.wait(
                        current_tasks, return_when=asyncio.FIRST_COMPLETED)
                    for t in done:
                        res = t.result()
                        if res:
                            homologies.extend(res)
            ts_queue.append(current_time)
        else:
            if len(current_tasks) >= MAX_CONCURRENT_TASKS:
                done, current_tasks = await asyncio.wait(
                    current_tasks, return_when=asyncio.FIRST_COMPLETED)
                for t in done:
                    res = t.result()
                    if res:
                        homologies.extend(res)

        task = loop.create_task(get_homology(gene_id))
        current_tasks.add(task)

    if current_tasks:
        done, _ = await asyncio.wait(current_tasks)
        for t in done:
            res = t.result()
            if res:
                homologies.extend(res)

    return homologies


@cli.command()
@click.argument('path')
@click.argument('output_path')
@click.option('--fast', 'fast_download', is_flag=True)
def homologies(path, output_path, fast_download):
    click.echo('Parsing GTF file')
    df = parse_gtf(path)
    click.echo('{} genes loaded'.format(len(df)))

    click.echo('Extracting genome homologies')
    loop = asyncio.get_event_loop()
    homologies = loop.run_until_complete(
        get_homologies(loop, df.id.values, fast_download))

    click.echo('Storing homologies to {}'.format(output_path))
    pd.DataFrame(homologies).to_csv(output_path, index=False)


if __name__ == '__main__':
    cli()
