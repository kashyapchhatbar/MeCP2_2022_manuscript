#!/usr/bin/env python3

# Convert bsmap custom output to bedGraph
# Remember bsmap coordinate is 1-based while bed start is zero-based
# https://bedtools.readthedocs.io/en/latest/content/overview.html#bed-starts-are-zero-based-and-bed-ends-are-one-based

# Output format: tab delimited txt file with the following columns:
#     1) chromorome
#     2) coordinate (1-based)
#     3) strand
#     4) sequence context (2nt upstream to 2nt downstream in Watson strand direction)
#     5) methylation ratio
#     6) number of converted and unconverted Cs covering this locus 
#     7) number of unconverted Cs covering this locus
#     8) lower bound of 95% confidence interval of methylation ratio
#     9) upper bound of 95% confidence interval of methylation ratio

# For some reason, the output files are in this format
#     1) chromosome
#     2) coordinate (1-based) # hopefully
#     3) strand
#     4) sequence context # CHH, CHG or CG (useless as far as we are concerned)
#     5) number of converted and unconverted Cs covering this locus 
#     6) number of unconverted Cs covering this locus

import numpy as np
import pandas as pd
import logging
import coloredlogs
import click
import gzip

np.seterr(divide='ignore')
logger = logging.getLogger("CX BEDGRAPH")
coloredlogs.install(level="INFO")
       
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--coverage', default=5, help="Coverage Threshold")
@click.option('--output_prefix', default='output', help='Provide output prefix')
def cx_bedgraph(input_file, coverage, output_prefix):
    wH = gzip.open("{}.bedGraph.gz".format(output_prefix), "wt")
    logger.info("Reading {}".format(input_file))    
    with gzip.open(input_file, "rb") as gO:
        for i,j in enumerate(gO, start=1):
            _sca, _start, _end, perc, m, u = j.strip().split(b"\t")            
            sca = _sca.decode()
            cov = int(u) + int(m)         
            meth = int(m)
            end = int(_end)
            start = end - 1
            if cov >= coverage:                    
                print(sca, start, end, round((meth/cov), 2), sep="\t", file=wH)
            if i % 10000000 == 0:
                logger.info("{:,} lines processed".format(i))
    logger.info("Finished Reading {}".format(input_file))
    wH.close()    

if __name__ == '__main__':
    cx_bedgraph()
