#!/usr/bin/env python
'''
genetic correlation estimation with individual data and summary stats

GENJI

Created on 2020-10-25

@author: Yiliang Zhang
'''

import argparse, os.path, sys
import pandas as pd
import numpy as np
from prep import prep
from ggrscore import ggrscore
from calculate import calculate


try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('GENJI requires pandas version > 0.15.2')


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


# returns whether the parent directory of path exists
def parent_dir_exists(path):
    return os.path.exists(os.path.abspath(os.path.join(path, os.pardir)))

def pipeline(args):
    pd.options.mode.chained_assignment = None

    # Sanity check args
    if not parent_dir_exists(args.out):
        raise ValueError('--out flag points to an invalid path.')

    print('Preparing files for analysis...')
    gwas_snps, ggr_df, N2 = prep(args.bfile, args.genotype, args.sumstats, args.N2, args.phenotype, args.covariates, args.chr, args.start, args.end)
    m = len(gwas_snps)
    print('{} SNPs included in our analysis...'.format(m))
    if args.ovpunknown and args.ovp is not None:
        raise ValueError('--ovp and --ovpunknown should not be simutaneously set.')
    if not args.ovpunknown and args.intercept is not None:
        raise ValueError('When the information of overlapped samples is known, --intercept should not be provided.')
    intercept, intercept_se = ggrscore(args.bfile, args.genotype, gwas_snps, args.ovp, ggr_df, args.ovpunknown, args.intercept, N2, args.h1, args.h2)
    print('Calculating genetic covariance...')
    out = calculate(ggr_df, args.h1, args.h2, intercept, intercept_se, N2, m)
    out.to_csv(args.out, sep=' ', na_rep='NA', index=False)


parser = argparse.ArgumentParser()

parser.add_argument('phenotype',
    help='The individual-level phenotype data for the first trait.')
parser.add_argument('genotype',
    help='Prefix for Plink .bed/.bim/.fam file of individual-level genotype data for the first trait.')
parser.add_argument('sumstats',
    help='The sumstats file for the second trait.')

parser.add_argument('--bfile', required=True, type=str,
    help='Prefix for Plink .bed/.bim/.fam file of reference panel for the second trait.')
parser.add_argument('--h1', required=True, type=float,
    help='The estimated heritability of the first trait.')
parser.add_argument('--h2', required=True, type=float,
    help='The estimated heritability of the second trait.')
parser.add_argument('--N2', type=int,
    help='N of the sumstats file for the second trait. If not provided, this value will be inferred from the sumstats arg.')
parser.add_argument('--ovp', type=str,
    help='Text file indicating the overlapping samples between the two GWASs. If not provided, the method will assume no sample overlap.')
parser.add_argument('--covariates', type=str,
    help='Text file for the covariates of the individuals included in the first trait if provided.')
parser.add_argument('--chr', type=int,
    help='')
parser.add_argument('--start', type=int,
    help='')
parser.add_argument('--end', type=int,
    help='')
parser.add_argument('--ovpunknown', action='store_true',
    help='If the information of the overlapped samples is unknown.')
parser.add_argument('--intercept', type=float,
    help='If the information of the overlapped samples is unknown, you can provide the intercept estimation of LDSC. If not provided, GENJI will estimate it.')
parser.add_argument('--out', required=True, type=str,
    help='Location to output results.')

if __name__ == '__main__':
    pipeline(parser.parse_args())

