#!/usr/bin/env python

from __future__ import division, print_function
import ld.ldscore as ld
import ld.parse as ps
import numpy as np
import pandas as pd


def _remove_dtype(x):
    '''Removes dtype: float64 and dtype: int64 from pandas printouts'''
    x = str(x)
    x = x.replace('\ndtype: int64', '')
    x = x.replace('\ndtype: float64', '')
    return x


def loj_bim(filter_df, array):
    r = filter_df.columns[1]
    l = array.IDList.columns[0]
    merge_df = filter_df.iloc[:,[1]]
    merge_df.loc[:,'keep'] = True
    z = pd.merge(array.IDList, merge_df, how='left', left_on=l, right_on=r, sort=False)
    ii = z['keep'] == True
    return ii.to_numpy().nonzero()[0]


def loj_fam(filter_df, array):
    r = filter_df.columns[0]
    l = array.IDList.columns[0]
    merge_df = filter_df.iloc[:,[0]]
    merge_df.loc[:,'keep'] = True
    z = pd.merge(array.IDList, merge_df, how='left', left_on=l, right_on=r, sort=False)
    ii = z['keep'] == True
    return ii.to_numpy().nonzero()[0]


def __filter_bim__(filter_df, array):
    merged_list = loj_bim(filter_df, array)
    len_merged_list = len(merged_list)
    if len_merged_list > 0:
        c = 'After merging, {0} SNPs remain'
    else:
        error_msg = 'No SNPs retained for analysis'
        raise ValueError(error_msg)
    return merged_list


def __filter_fam__(filter_df, array):
    merged_list = loj_fam(filter_df, array)
    len_merged_list = len(merged_list)
    if len_merged_list > 0:
        c = 'After merging, {0} individuals in GWAS 1 remain'
    else:
        error_msg = 'No individuals retained for analysis'
        raise ValueError(error_msg)
    return merged_list


def subset_annot_file(a_df, GWAS_df, kept_cols):
    GWAS_df.loc[:,'idx'] = pd.Series(range(len(GWAS_df.SNP.values)))
    a_df = pd.merge(a_df, GWAS_df, on=['SNP'])
    a_df = a_df.sort_values(['idx'])
    a_df.drop('idx', axis=1, inplace=True)
    a_df.rename(columns={'CHR_x':'CHR', 'BP_x':'BP', 'CM_x':'CM'}, inplace=True)
    a_df = a_df.iloc[:,0:kept_cols]
    return a_df


def remove_brackets(x):
    return x.replace('[', '').replace(']', '').strip()


def _ggrscore(bfile, genotype, phenotype_data, gwas_snps):

    snp_file, fsnp_file, snp_obj = bfile+'.bim', genotype+'.bim', ps.PlinkBIMFile
    ind_file, find_file, ind_obj = bfile+'.fam', genotype+'.fam', ps.PlinkFAMFile
    array_file, farray_file, array_obj = bfile+'.bed', genotype+'.bed', ld.PlinkBEDFile

    # read bim/snp
    array_snps = snp_obj(snp_file)
    farray_snps = snp_obj(fsnp_file)

    # snp list

    keep_snps_ref = __filter_bim__(gwas_snps, array_snps)
    keep_snps_genotype = __filter_bim__(gwas_snps, farray_snps)

    # read fam
    array_indivs = ind_obj(ind_file)
    genotype_indivs = ind_obj(find_file)

    # phenotype data
    phenotype_info = pd.merge(genotype_indivs.IDList, phenotype_data, on='IID')
    ii = phenotype_info['Phenotype'] != 9
    pheno_avg = np.mean(phenotype_info['Phenotype'][ii])
    phenotype_info['Phenotype'][np.logical_not(ii)] = pheno_avg
    phenotype_denom = np.std(phenotype_info['Phenotype'])
    if phenotype_denom == 0:
        phenotype_denom = 1
    phenotype_info['Phenotype'] = (phenotype_info['Phenotype'] - pheno_avg) / phenotype_denom

    n = len(array_indivs.IDList)
    m = len(genotype_indivs.IDList)
    # read keep_indivs
    keep_indivs_ref = None
    keep_indivs_genotype = __filter_fam__(phenotype_info, genotype_indivs)
    # read genotype array
    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps_ref,
        keep_indivs=keep_indivs_ref, mafMin=None)
    geno_farray = array_obj(farray_file, m, farray_snps, keep_snps=keep_snps_genotype,
        keep_indivs=keep_indivs_genotype, mafMin=None)

    #determine block widths

    max_dist = 1
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    block_left = ld.getBlockLefts(coords, max_dist)

    y, ggr, cor_sum = geno_array.ggrscoreVarBlocks(geno_farray, phenotype_info, gwas_snps, block_left, 50)

    # print .ldscore. Output columns: CHR, BP, RS, [LD Scores]
    df = pd.DataFrame({'y':y, 'ggr':ggr, 'l2':cor_sum})

    return df, len(phenotype_info)


def ggrscore(bfile, genotype, phenotype, gwas_snps):
    
    # read phenotype data
    phenotype_data = pd.read_csv(phenotype, header=None, names=['FID', 'IID', 'Phenotype'], delim_whitespace=True)

    df = None
    if '@' in bfile:
        all_dfs = []
        N1 = float('inf')
        for i in range(1, 23):
            cur_bfile = bfile.replace('@', str(i))
            cur_gwas_snps = gwas_snps[gwas_snps.iloc[:,0]==i].reset_index(drop=True)
            if len(cur_gwas_snps) == 0:
                continue
            if '@' in genotype:
                cur_genotype = genotype.replace('@', str(i))
                cur_df, cur_N = _ggrscore(cur_bfile, cur_genotype, phenotype_data, cur_gwas_snps)
                if cur_N < N1:
                    N1 = cur_N
                all_dfs.append(cur_df)
            else:
                cur_df, cur_N = _ggrscore(cur_bfile, genotype, phenotype_data, cur_gwas_snps)
                if cur_N < N1:
                    N1 = cur_N
                all_dfs.append(cur_df)
            print('Computed LD scores for chromosome {}'.format(i))
        df = pd.concat(all_dfs)
    else:
        if '@' in genotype:
            all_dfs = []
            N1 = float('inf')
            for i in range(1, 23):
                cur_gwas_snps = gwas_snps[gwas_snps.iloc[:,0]==i].reset_index(drop=True)
                if len(cur_gwas_snps) == 0:
                    continue
                cur_genotype = genotype.replace('@', str(i))
                cur_df, cur_N = _ggrscore(bfile, cur_genotype, phenotype_data, cur_gwas_snps)
                if cur_N < N1:
                    N1 = cur_N
                all_dfs.append(cur_df)
                print('Computed LD scores for chromosome {}'.format(i))
            df = pd.concat(all_dfs)
        else:
            df, N1 = _ggrscore(bfile, genotype, phenotype_data, gwas_snps)
    return df, N1
