# GENJI

GENJI (Genetic-covariance EstimatioN Jointly using Individual-level and summary data) is a statistical framework to perform genetic covariance analysis. GENJI needs individual-level data for one trait and summary data for another to estimate genetic covariance.

```
python3 genji.py  phenotype_data_for_trait_1 genotype_data_for_trait_1 summary_data_for_trait_2\
--bfile ref_panel\
--h1 di_yi_ge_gwas_de_heritability_ke_yi_yong_reml_huo_zhe_ldsc_gu_ji\
--h2 di_er_ge_gwas_de_heritability_ke_yi_yong_ldsc_xian_gu_ji_hao\
--N2 zhi_xu_yao_di_er_ge_gwas_de_sample_size\
--out shu_chu_de_wei_zhi
```

plink bfile ru guo shi an zhao chromosome fen kai de, na me ke yi yong '@' lai dai ti shu zi. bi ru: eur_chr@_SNPmaf5.