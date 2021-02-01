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


def _ggrscore(bfile, genotype, gwas_snps, ggr_df, ovp_sample, N2):

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

    n = len(array_indivs.IDList)
    m = len(genotype_indivs.IDList)
    # read keep_indivs
    keep_indivs_ref = None
    keep_indivs_genotype = __filter_fam__(ggr_df, genotype_indivs)
    # read genotype array
    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps_ref,
        keep_indivs=keep_indivs_ref, mafMin=None)
    geno_farray = array_obj(farray_file, m, farray_snps, keep_snps=keep_snps_genotype,
        keep_indivs=keep_indivs_genotype, mafMin=None)
    ovp_index = pd.merge(genotype_indivs.IDList, ovp_sample, how='left', on='IID')['ovp'] == True
    #determine block widths

    max_dist = 1
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    block_left = ld.getBlockLefts(coords, max_dist)

    tmp_ggr = pd.DataFrame({'IID': genotype_indivs.IDList['IID'][geno_farray.keep_indivs], 'gg': 0, 
        'grg': 0, 'ggg': 0, 'gz': 0, 'grrg': 0})
    geno_array.ggrscoreVarBlocks(geno_farray, gwas_snps, tmp_ggr, ovp_index, N2, block_left, 50)
    tmp_ggr = pd.merge(ggr_df[['IID']], tmp_ggr, on='IID')
    ggr_df['gg'] += tmp_ggr['gg']
    ggr_df['grg'] += tmp_ggr['grg']
    ggr_df['ggg'] += tmp_ggr['ggg']
    ggr_df['gz'] += tmp_ggr['gz']
    ggr_df['grrg'] += tmp_ggr['grrg']

def ggrscore(bfile, genotype, gwas_snps, ovp, ggr_df, N2):
    if ovp is None:
        ovp_sample = pd.DataFrame({'IID':[]})
    else:
        ovp_sample = pd.read_csv(ovp, header=None, names=['IID'], delim_whitespace=True)
    
    ovp_sample = pd.merge(ovp_sample, ggr_df[['IID']], on='IID')
    ovp_sample['ovp'] = True
    if '@' in bfile:
        for i in range(1, 23):
            cur_bfile = bfile.replace('@', str(i))
            cur_gwas_snps = gwas_snps[gwas_snps.iloc[:,0]==i].reset_index(drop=True)
            if len(cur_gwas_snps) == 0:
                continue
            if '@' in genotype:
                cur_genotype = genotype.replace('@', str(i))
                _ggrscore(cur_bfile, cur_genotype, cur_gwas_snps, ggr_df, ovp_sample, N2)
            else:
                _ggrscore(cur_bfile, genotype, cur_gwas_snps, ggr_df, ovp_sample, N2)
            print('Done with SNPs in chromosome {}'.format(i))
    else:
        if '@' in genotype:
            for i in range(1, 23):
                cur_gwas_snps = gwas_snps[gwas_snps.iloc[:,0]==i].reset_index(drop=True)
                if len(cur_gwas_snps) == 0:
                    continue
                cur_genotype = genotype.replace('@', str(i))
                _ggrscore(bfile, cur_genotype, cur_gwas_snps, ggr_df, ovp_sample, N2)
                print('Done with SNPs in chromosome {}'.format(i))
        else:
            _ggrscore(bfile, genotype, gwas_snps, ggr_df, ovp_sample, N2)
    ggr_df = ggr_df.merge(ovp_sample, on='IID', how='left')
    ggr_df['ovp'] = ggr_df['ovp'] == True
