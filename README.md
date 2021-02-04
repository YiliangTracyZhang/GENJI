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
python3 genji.py  [phenotype data for trait 1] [genotype data for trait 1] [summary stats for trait 2]\
--bfile [genotype data for reference panel]\
--h1 [heritability of trait 1] (optinal)\
--h2 [heritability of trait 2] (optinal)\
--N2 [sample size of trait 2] (optinal)\
--ovp [overlapped samples in the two studies] (optional)\
--re [covariance of environmental errors between trait 1 and trait 2] (optional)\
--chr [chromosome of the genomic region] (optional)\
--start [start point of the genomic region]\
--end [end point of the genomic region]\
--out [file location for the results]