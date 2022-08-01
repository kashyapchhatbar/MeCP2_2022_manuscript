#!/usr/bin/env python3

import numpy as np
import pandas as pd
import logging
import coloredlogs
import tarfile
import click
import gzip
from pyfaidx import Faidx

np.seterr(divide='ignore')
logger = logging.getLogger("CX BEDGRAPH")
coloredlogs.install(level="INFO")

chr_conversion = {'mm9': 'https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCm37_ensembl2UCSC.txt',
                 'mm10': 'https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCm38_ensembl2UCSC.txt'}
        
@click.command()
@click.argument('input_files', nargs=-1, type=click.Path(exists=True))
@click.option('--coverage', default=5, help="Coverage Threshold")
@click.option('--output_prefix', default='output', help='Provide output prefix')
@click.option('--genome', default='mm9', help='Provide genome assembly')
@click.option('--context', default='CA', help='Trinucleotide context')
def cx_bedgraph(input_files, coverage, output_prefix, genome, context):
    
    context = context.encode()
    wH = gzip.open("{}.preprocessed.gz".format(output_prefix), "wt")    
    genome_chromosomes = pd.read_csv(chr_conversion[genome], sep="\t", header=None, index_col=0).dropna()
    if genome == "mm10":
        genome_chromosomes.index = genome_chromosomes[1].astype(bytes)
    else:
        genome_chromosomes.index = genome_chromosomes.index.astype(bytes)
    genome_chromosomes = genome_chromosomes[1].to_dict()
    
    mm9 = Faidx(f"/datastore/homes3/genomes/mouse/{genome}/{genome}.fa")
    
    for input_file in input_files:
        logger.info("Reading {}".format(input_file))
        with tarfile.open(input_file, "r:gz") as tar:
            file_names = tar.getnames()
            for file_name in file_names:
                with tar.extractfile(file_name) as gO:
                    header = next(gO)
                    for i,j in enumerate(gO, start=1):
                        _sca, _end, strand, fcontext, m, c, _ = j.strip().split(b"\t")
                        sca = genome_chromosomes.get(_sca, None) # convert 1 to chr1 for example            
                        if sca:
                            cov = int(c)
                            meth = int(m)
                            if strand == b"+":                    
                                start = int(_end) + 1
                                end = start + 2                    
                            elif strand == b"-":
                                start = int(_end) - 1
                                end = start + 2                    
                            if cov >= int(coverage) and fcontext.startswith(context.upper()):
                                # print(sca, start, end, round((meth/cov), 4), sep="\t", file=wH)
                                print(sca, start, end, meth, cov, fcontext.decode(), sep="\t", file=wH)
                        if i % 10000000 == 0:
                            logger.info("{:,} lines processed".format(i))
            logger.info("Finished Reading {}".format(input_file))
    wH.close()    

if __name__ == '__main__':
    cx_bedgraph()
