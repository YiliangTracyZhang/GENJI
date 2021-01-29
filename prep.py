import os

import numpy as np
import pandas as pd

def allign_alleles(df):
    """Look for reversed alleles and inverts the z-score for one of them.

    Here, we take advantage of numpy's vectorized functions for performance.
    """
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    a = []  # array of alleles
    for colname in ['A1_ref', 'A2_ref', 'A1_gen', 'A2_gen', 'A1_y', 'A2_y']:
        tmp = np.empty(len(df[colname]), dtype=int)
        for k, v in d.items():
            tmp[np.array(df[colname]) == k] = v
        a.append(tmp)
    matched_alleles_gen = (((a[0] == a[2]) & (a[1] == a[3])) |
        ((a[0] == 3 - a[2]) & (a[1] == 3 - a[3])))
    reversed_alleles_gen = (((a[0] == a[3]) & (a[1] == a[2])) |
        ((a[0] == 3 - a[3]) & (a[1] == 3 - a[2])))
    matched_alleles_y = (((a[0] == a[4]) & (a[1] == a[5])) |
        ((a[0] == 3 - a[4]) & (a[1] == 3 - a[5])))
    reversed_alleles_y = (((a[0] == a[5]) & (a[1] == a[4])) |
        ((a[0] == 3 - a[5]) & (a[1] == 3 - a[4])))
    df['Z_y'] *= -2 * reversed_alleles_y + 1
    df['reversed'] = reversed_alleles_gen
    df = df[((matched_alleles_y|reversed_alleles_y)&(matched_alleles_gen|reversed_alleles_gen))]


def get_files(file_name):
    if '@' in file_name:
        valid_files = []
        for i in range(1, 23):
            cur_file = file_name.replace('@', str(i))
            if os.path.isfile(cur_file):
                valid_files.append(cur_file)
            else:
                raise ValueError('No file matching {} for chr {}'.format(
                    file_name, i))
        return valid_files
    else:
        if os.path.isfile(file_name):
            return [file_name]
        else:
            ValueError('No files matching {}'.format(file_name))


def prep(bfile, genotype, sumstats2, N2, phenotype):
    bim_files = get_files(bfile + '.bim')
    genotype_files = get_files(genotype + '.bim')
    # read in bim files
    bims = [pd.read_csv(f,
                        header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                        delim_whitespace=True) for f in bim_files]
    bim = pd.concat(bims, ignore_index=True)

    genotype_bims = [pd.read_csv(f,
                        header=None,
                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                        delim_whitespace=True) for f in genotype_files]
    genotype_bim = pd.concat(genotype_bims, ignore_index=True)

    summary_stats = pd.read_csv(sumstats2, delim_whitespace=True)

    # rename cols
    bim.rename(columns={'CHR': 'CHR_ref', 'CM': 'CM_ref', 'BP':'BP_ref', 'A1': 'A1_ref', 'A2': 'A2_ref'}, inplace=True)
    genotype_bim.rename(columns={'CHR': 'CHR_gen', 'CM': 'CM_gen', 'BP':'BP_gen', 'A1': 'A1_gen', 'A2': 'A2_gen'}, inplace=True)
    summary_stats.rename(columns={'A1': 'A1_y', 'A2': 'A2_y', 'N': 'N_y', 'Z': 'Z_y'}, inplace=True)

    # take overlap between output and ref genotype files
    df = pd.merge(bim, genotype_bim, on=['SNP']).merge(summary_stats, on=['SNP'])
    # flip sign of z-score for allele reversals
    allign_alleles(df)
    df = df.drop_duplicates(subset='SNP', keep=False)
    
    if N2 is not None:
        N2 = N2
    else:
        N2 = summary_stats['N_y'].max()
    df.rename(columns={'CHR_ref':'CHR'}, inplace=True)

    ggr_df = pd.read_csv(phenotype, header=None, names=['IID'], delim_whitespace=True, usecols=[1])
    fam_files = get_files(genotype + '.fam')
    for i in range(len(fam_files)):
        fam_data = pd.read_csv(fam_files[i], header=None, names=['IID'], delim_whitespace=True, usecols=[1])
        ggr_df = pd.merge(fam_data, ggr_df, on=['IID'])
    ggr_df['gg'] = 0
    ggr_df['grg'] = 0
    ggr_df['ggg'] = 0
    ggr_df['gz'] = 0
    ggr_df['grrg'] = 0
    return df[['CHR', 'SNP', 'Z_y', 'reversed']], ggr_df, N2
