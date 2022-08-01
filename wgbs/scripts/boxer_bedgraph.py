#!/usr/bin/env python3

import gzip
import click
import logging
import coloredlogs
import numpy as np

np.seterr(divide='ignore')
logger = logging.getLogger("CX BEDGRAPH")
coloredlogs.install(level="INFO")

@click.command()
@click.argument('input_file', type=click.File('rb'))
@click.option('--coverage', default=5, help="Coverage Threshold")
@click.option('--output_prefix', default='output', help='Provide output prefix')
def cx_bedgraph(input_file, coverage, output_prefix):
    wH = gzip.open("{}.bedGraph.gz".format(output_prefix), "wt")
    for i,j in enumerate(input_file, start=1):
        sca, start, end, m, c = j.strip().split(b"\t")
        meth = int(m)
        cov = int(c)
        if cov >= int(coverage):
            print(sca.decode(), start.decode(), end.decode(), round((meth/cov), 4), sep="\t", file=wH)
        if i % 10000000 == 0:
            logger.info("{:,} lines processed".format(i))
    logger.info("Finished Reading {}".format(input_file))
    wH.close()

if __name__ == '__main__':
    cx_bedgraph()