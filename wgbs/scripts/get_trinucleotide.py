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
from pyfaidx import Faidx

np.seterr(divide='ignore')
logger = logging.getLogger("CX BEDGRAPH")
coloredlogs.install(level="INFO")

chr_conversion = {'mm9': 'https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCm37_ensembl2UCSC.txt',
                 'mm10': 'https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCm38_ensembl2UCSC.txt'}
        
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--coverage', default=5, help="Coverage Threshold")
@click.option('--output_prefix', default='output', help='Provide output prefix')
@click.option('--genome', default='mm9', help='Provide genome assembly')
@click.option('--context', default='CAC', help='Trinucleotide context')
def cx_bedgraph(input_file, coverage, output_prefix, genome, context):    
    context = context.encode()
    wH = gzip.open("{}.bedGraph.gz".format(output_prefix), "wt")
    logger.info("Reading {}".format(input_file))
    genome_chromosomes = pd.read_csv(chr_conversion[genome], sep="\t", header=None, index_col=0).dropna()
    if genome == "mm10":
        genome_chromosomes.index = genome_chromosomes[1].astype(bytes)
    else:
        genome_chromosomes.index = genome_chromosomes.index.astype(bytes)
    genome_chromosomes = genome_chromosomes[1].to_dict()
    
    mm9 = Faidx(f"/datastore/homes3/genomes/mouse/{genome}/{genome}.fa")
    
    with gzip.open(input_file, "rb") as gO:
        for i,j in enumerate(gO, start=1):
            _sca, _end, strand, fcontext, m, c = j.strip().split(b"\t")
            sca = genome_chromosomes.get(_sca, None) # convert 1 to chr1 for example            
            if sca:
                cov = int(c)
                meth = int(m)
                if strand == b"+":                    
                    start = int(_end) + 1
                    end = start + 2
                    _context = str(mm9.fetch(sca, start, end)).encode().upper()
                elif strand == b"-":
                    start = int(_end) - 1
                    end = start + 2
                    _context = str(mm9.fetch(sca, start, end).reverse.complement).encode().upper()                
                if cov >= int(coverage) and _context.startswith(context.upper()):
                    print(sca, start, end, round((meth/cov), 4), sep="\t", file=wH)
            if i % 10000000 == 0:
                logger.info("{:,} lines processed".format(i))
    logger.info("Finished Reading {}".format(input_file))
    wH.close()    

if __name__ == '__main__':
    cx_bedgraph()
