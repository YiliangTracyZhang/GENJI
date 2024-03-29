import os

import numpy as np
import pandas as pd
from sklearn import linear_model

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
    df.where(pd.Series(((matched_alleles_y|reversed_alleles_y)&(matched_alleles_gen|reversed_alleles_gen))), inplace=True)
    df.dropna(inplace=True)


def get_files(file_name, chr):
    if '@' in file_name:
        valid_files = []
        if chr is None:
            for i in range(1, 23):
                cur_file = file_name.replace('@', str(i))
                if os.path.isfile(cur_file):
                    valid_files.append(cur_file)
                else:
                    raise ValueError('No file matching {} for chr {}'.format(
                        file_name, i))
        else:
            cur_file = file_name.replace('@', chr)
            if os.path.isfile(cur_file):
                valid_files.append(cur_file)
            else:
                raise ValueError('No file matching {} for chr {}'.format(
                    file_name, chr))
        return valid_files
    else:
        if os.path.isfile(file_name):
            return [file_name]
        else:
            ValueError('No files matching {}'.format(file_name))


def prep(bfile, genotype, sumstats2, N2, phenotype, covariates, chr, start, end):
    bim_files = get_files(bfile + '.bim', chr)
    genotype_files = get_files(genotype + '.bim', chr)
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
    if chr is not None:
        if start is None:
            start = 0
        if end is None:
            end = float('inf')
        genotype_bim = genotype_bim[np.logical_and(np.logical_and(genotype_bim['CHR']==chr, genotype_bim['BP']<=end), genotype_bim['BP']>=start)].reset_index(drop=True)
        bim = bim[np.logical_and(np.logical_and(bim['CHR']==chr, bim['BP']<=end), bim['BP']>=start)].reset_index(drop=True)
    summary_stats = pd.read_csv(sumstats2, delim_whitespace=True)

    # rename cols
    bim.rename(columns={'CHR': 'CHR_ref', 'CM': 'CM_ref', 'BP':'BP_ref', 'A1': 'A1_ref', 'A2': 'A2_ref'}, inplace=True)
    genotype_bim.rename(columns={'CHR': 'CHR_gen', 'CM': 'CM_gen', 'BP':'BP_gen', 'A1': 'A1_gen', 'A2': 'A2_gen'}, inplace=True)
    summary_stats.rename(columns={'A1': 'A1_y', 'A2': 'A2_y', 'N': 'N_y', 'Z': 'Z_y'}, inplace=True)

    # take overlap between output and ref genotype files
    df = pd.merge(bim, genotype_bim, on=['SNP']).merge(summary_stats, on=['SNP'])
    df = df[df['CHR_ref']==df['CHR_gen']]
    # flip sign of z-score for allele reversals
    allign_alleles(df)
    df = df.drop_duplicates(subset='SNP', keep=False).reset_index(drop=True)
    if N2 is not None:
        N2 = N2
    else:
        N2 = summary_stats['N_y'].max()
    df.rename(columns={'CHR_ref':'CHR'}, inplace=True)

    ggr_df = pd.read_csv(phenotype, header=None, names=['IID', 'Phenotype'], delim_whitespace=True, usecols=[1, 2])
    fam_files = get_files(genotype + '.fam', chr)
    for i in range(len(fam_files)):
        fam_data = pd.read_csv(fam_files[i], header=None, names=['IID'], delim_whitespace=True, usecols=[1])
        ggr_df = pd.merge(fam_data, ggr_df, on=['IID'])
    ii = ggr_df['Phenotype'] != 9
    pheno_avg = np.mean(ggr_df['Phenotype'][ii])
    ggr_df['Phenotype'][np.logical_not(ii)] = pheno_avg
    ggr_df['Phenotype'] = ggr_df['Phenotype'] - pheno_avg
    if covariates is not None:
        covariates_df = pd.read_csv(covariates, header=None, delim_whitespace=True)
        covariates_df = covariates_df.iloc[:, 1:]
        covariates_df = covariates_df.rename(columns={1:'IID'})
        regression_df = pd.merge(ggr_df[['IID']], covariates_df, on=['IID'])
        ggr_df = pd.merge(ggr_df, covariates_df[['IID']], on=['IID'])
        colnames = regression_df.columns
        for jj in range(1, len(colnames)):
            cur_col = colnames[jj]
            iii = regression_df[cur_col] != 9
            temp_avg = np.mean(regression_df[cur_col][iii])
            regression_df[cur_col][np.logical_not(iii)] = temp_avg
        lm = linear_model.LinearRegression().fit(regression_df.iloc[:, 1:], ggr_df[['Phenotype']])
        ggr_df[['Phenotype']] = ggr_df[['Phenotype']] - lm.predict(regression_df.iloc[:,1:])


    phenotype_denom = np.std(ggr_df['Phenotype'])
    ggr_df['Phenotype'] = ggr_df['Phenotype'] / phenotype_denom

    ggr_df['gg'] = 0.0
    ggr_df['grg'] = 0.0
    ggr_df['ggg'] = 0.0
    ggr_df['gz'] = 0.0
    ggr_df['grrg'] = 0.0
    return df[['CHR', 'SNP', 'Z_y', 'reversed']], ggr_df, N2
