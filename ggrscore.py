#!/usr/bin/env python

from __future__ import division, print_function
from code import interact
import ld.ldscore as ld
import ld.parse as ps
import numpy as np
import pandas as pd

from sklearn import linear_model

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


def _ggrscore(bfile, genotype, gwas_snps, ggr_df, ovp_sample, N2, intercept):

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

    tmp_ggr = pd.DataFrame({'IID': genotype_indivs.IDList['IID'][geno_farray.keep_indivs], 'gg': 0.0, 
        'grg': 0.0, 'ggg': 0.0, 'gz': 0.0, 'grrg': 0.0})
    tmp_ggr = tmp_ggr.merge(ggr_df[['IID', 'Phenotype']], on='IID')
    geno_array.ggrscoreVarBlocks(geno_farray, gwas_snps, tmp_ggr, ovp_index, N2, block_left, 50, intercept)
    tmp_ggr = pd.merge(ggr_df[['IID']], tmp_ggr, on='IID')
    ggr_df['gg'] += tmp_ggr['gg']
    ggr_df['grg'] += tmp_ggr['grg']
    ggr_df['ggg'] += tmp_ggr['ggg']
    ggr_df['gz'] += tmp_ggr['gz']
    ggr_df['grrg'] += tmp_ggr['grrg']

def __intercept__(zz, rxx, xxxx, l, h1, h2, N1, N2):
    m = len(zz)
    rho = np.sum(zz) / np.sqrt(N1 * N2)
    w1 = 1 - h1 + N1 * h1 * xxxx / m
    w2 = 1 + N2 * h2 * l / m
    w3 = np.sqrt(N1 * N2) * rho * rxx / m
    w = 1 / ((w1 * w2 + w3 * w3) * rxx ** 1.5)
    w[(w < 0) | (w == np.inf) | (w == -np.inf)] = 0
    ldsc_model = linear_model.LinearRegression().fit(pd.DataFrame(rxx), pd.DataFrame(zz), sample_weight=w)
    intercept = ldsc_model.intercept_[0]
    
    rho = ldsc_model.coef_[0][0] * m / np.sqrt(N1 * N2)
    w3 = np.sqrt(N1 * N2) * rho * rxx / m + intercept
    w = 1 / ((w1 * w2 + w3 * w3) * rxx ** 1.5)
    w[(w < 0) | (w == np.inf) | (w == -np.inf)] = 0
    nblock = 200
    intercept_block = np.empty(nblock)
    rho_block = np.empty(nblock)
    bsize = m // nblock
    for j in range(nblock):
        rxx_b = np.hstack([rxx[:j*bsize], rxx[(j+1)*bsize:]])
        zz_b = np.hstack([zz[:j*bsize], zz[(j+1)*bsize:]])
        w_b = np.hstack([w[:j*bsize], w[(j+1)*bsize:]])
        ldsc_model = linear_model.LinearRegression().fit(pd.DataFrame(rxx_b), pd.DataFrame(zz_b), sample_weight=w_b)
        intercept_block[j] = ldsc_model.intercept_[0]
        rho_block[j] = ldsc_model.coef_[0] * m /np.sqrt(N1 * N2)
    return np.mean(intercept_block), np.sqrt(np.var(rho_block) * (nblock - 1))

def ggrscore(bfile, genotype, gwas_snps, ovp, ggr_df, ovpunknown, intercept, N2, h1, h2):
    if ovp is None:
        ovp_sample = pd.DataFrame({'IID':[]})
        if not ovpunknown:
            intercept = 0.0
    else:
        ovp_sample = pd.read_csv(ovp, header=None, names=['IID'], delim_whitespace=True)
        intercept = 0.0
    
    ovp_sample = pd.merge(ovp_sample, ggr_df[['IID']], on='IID')
    ovp_sample['ovp'] = True
    if intercept is None:
        gwas_snps['Z_x'] = 0.0
        gwas_snps['rxx'] = 0.0
        gwas_snps['xxxx'] = 0.0
        gwas_snps['l'] = 0.0
    if '@' in bfile:
        for i in range(1, 23):
            cur_bfile = bfile.replace('@', str(i))
            cur_gwas_snps = gwas_snps[gwas_snps.iloc[:,0]==i].reset_index(drop=True)
            if len(cur_gwas_snps) == 0:
                continue
            if '@' in genotype:
                cur_genotype = genotype.replace('@', str(i))
                _ggrscore(cur_bfile, cur_genotype, cur_gwas_snps, ggr_df, ovp_sample, N2, intercept)
            else:
                _ggrscore(cur_bfile, genotype, cur_gwas_snps, ggr_df, ovp_sample, N2, intercept)
            if intercept is None:
                gwas_snps['Z_x'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['Z_x'].to_list()
                gwas_snps['rxx'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['rxx'].to_list()
                gwas_snps['xxxx'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['xxxx'].to_list()
                gwas_snps['l'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['l'].to_list()
            print('Done with SNPs in chromosome {}'.format(i))
    else:
        if '@' in genotype:
            for i in range(1, 23):
                cur_gwas_snps = gwas_snps[gwas_snps.iloc[:,0]==i].reset_index(drop=True)
                if len(cur_gwas_snps) == 0:
                    continue
                cur_genotype = genotype.replace('@', str(i))
                _ggrscore(bfile, cur_genotype, cur_gwas_snps, ggr_df, ovp_sample, N2, intercept)
                if intercept is None:
                    gwas_snps['Z_x'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['Z_x'].to_list()
                    gwas_snps['rxx'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['rxx'].to_list()
                    gwas_snps['xxxx'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['xxxx'].to_list()
                    gwas_snps['l'][gwas_snps.iloc[:,0]==i] = cur_gwas_snps['l'].to_list()
                print('Done with SNPs in chromosome {}'.format(i))
        else:
            _ggrscore(bfile, genotype, gwas_snps, ggr_df, ovp_sample, N2, intercept)

    ldsc_se = 0.0
    if intercept is None:
        zz = gwas_snps['Z_x'] * gwas_snps['Z_y']
        intercept, ldsc_se = __intercept__(zz, gwas_snps['rxx'], gwas_snps['xxxx'], gwas_snps['l'], h1, h2, len(ggr_df), N2)
    ggr_df['ovp'] = ggr_df[['IID']].merge(ovp_sample, on='IID', how='left')['ovp']
    ggr_df['ovp'] = ggr_df['ovp'] == True
    return intercept, ldsc_se
