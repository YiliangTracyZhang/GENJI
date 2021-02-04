# GENJI

GENJI (Genetic-covariance EstimatioN Jointly using Individual-level and summary data) is a statistical framework to perform genetic covariance analysis. GENJI needs individual-level data for one trait and summary data for another to estimate genetic covariance.

## Requirements

The software is developed and tested in Linux and Mac OS environments. The following softwares and packages are required:

1. **Python 3**
2. **numpy**
3. **scipy**
4. **pandas**
6. **bitarray**

## Tutorial

You can download GENJI by:

```
$ git clone https://github.com/YiliangTracyZhang/GENJI
$ cd ./GENJI
```

```
python3 genji.py  \[phenotype data for trait 1\] genotype_data_for_trait_1 summary_data_for_trait_2\
--bfile ref_panel\
--h1 di_yi_ge_gwas_de_heritability_ke_yi_yong_reml_huo_zhe_ldsc_gu_ji\
--h2 di_er_ge_gwas_de_heritability_ke_yi_yong_ldsc_xian_gu_ji_hao\
--N2 zhi_xu_yao_di_er_ge_gwas_de_sample_size\
--out shu_chu_de_wei_zhi
```

plink bfile ru guo shi an zhao chromosome fen kai de, na me ke yi yong '@' lai dai ti shu zi. bi ru: eur_chr@_SNPmaf5.