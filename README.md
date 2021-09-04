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

You can download SUPERGNOVA by:

```
$ git clone https://github.com/YiliangTracyZhang/GENJI
$ cd ./GENJI
```

You can download GENJI by:

```
$ git clone https://github.com/YiliangTracyZhang/GENJI
$ cd ./GENJI
```

We'll need a few types of files to implement GENJI:

- **Summary statistics file:** We assume that the file is in the standard format that ``ldsc`` understands. If not, please make sure to run them through the ``munge_sumstats.py`` file under Python2.7 included in ``ldsc`` (see [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics) for instructions). Before using GENJI, please remove all the SNPs with missing values in your GWAS summary data.

- **Plink bfiles:** These are files in .bed/.bim/.fam format containing genotype information of reference panel and individual-level for one of the GWASs.

- **Phenotype file:** Three-column space-separated txt file with the first, second, and third column the family ID, sample ID, and phenotype values, respectively. The file should not have header.

- **Overlap sample IDs:** This txt file only have one columns which contains the sample IDs of the overlapped samples.

You may run the following command:

```
python3 genji.py  [phenotype data for trait 1] [genotype data for trait 1] [summary stats for trait 2]\
--bfile [genotype data for reference panel]\
--h1 [heritability of trait 1]\
--h2 [heritability of trait 2]\
--N2 [sample size of trait 2] (optinal)\
--ovp [overlapped samples in the two studies] (optional)\
--out [file location for the results]
```

As a toy example, at current directory, you can run the following command

```
python3 genji.py ./input/BMI.phen ./input/nfbc_genotype_chr22 ./input/T2D.txt\ 
--bfile ./input/eur_chr22\ 
--covarites ./input/covariates_nfbc.txt\ 
--h1 0.171038\ 
--h2 0.2462\
--out ./output/BMI_T2D_chr22.txt
```

It should take about 30 seconds for this toy example to run. The correct output should be:



### Explanation of Command-Line Arguments

- The first two arguments denote the locations of the phenotype and genotype data. 

- Summary stats for trait 2 may be compressed using gzip, bz2, zip, xz, or not compressed at all. The program will infer the compression method if the files end with .gz, .bz2, .zip, xz, respectively. As previously mentioned, we assume that the files are in the standard format that `ldsc` understands.

- `h1` and `h2` arguments denote the heritability estimation for trait 1 and trait 2, respectively.

- The `N2` argument (optional) denote the sample size of the study only having summary stats. If they are not provided, they will be inferred from the summary statistics file.

- The `bfile` argument denotes the prefix of the `.bed/.bim/.fam` genotypic data files. Note the '@', which denotes a wildcard character that GENJI will be replaced with 1-22. Alternatively, if you have one set of genotypic data files with 22 chromosome combined, you can just specify one bfile.

- The `out` flag denotes the file location for the results to be outputted to.

### Explanation of Output
The output will be a whitespace-delimited text file, with the rows corresponding to different annotations and the columns as such:

- `rho`: The estimation of genetic covariance.
- `se`: The standard error of the estimate.
- `pvalue`: The p value of genetic covariance.
- `corr`: The estimation of genetic correlation.
- `h2_1`: The estimation of heritability of the first trait.
- `h2_2`: The estimation of heritability of the second trait.
- `m`: The number of SNPs involved in the estimation of genetic covariance.
- `N1`: The sample size of study 1.
- `N2`: The sample size of study 2.


## Credits

Those who use GENJI should cite [Zhang, Y.L. et al. Estimating genetic correlation jointly using individual-level and summary-level GWAS data. 2021.](https://www.biorxiv.org/content/10.1101/2021.08.18.456908v1)

The LD score calculation  and the estimation of phenotypic covariance are adapted from `ldsc.py` in  `ldsc` and `ldsc_thin.py` in `GNOVA`. See [Bulik-Sullivan, B. et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.](https://www.nature.com/articles/ng.3406) and [Lu, Q.S. et al. A powerful approach to estimating annotation-stratified genetic covariance using GWAS summary statistics. The American Journal of Human Genetics, 2017.](https://www.cell.com/ajhg/fulltext/S0002-9297(17)30453-6)


