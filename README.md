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

- **Plink bfiles:** These are files in .bed/.bim/.fam format containing genotype information of reference panel and individual-level for study 2.

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



