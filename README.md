# easyGWAS

R interface to perform GWAS using Factored Spectrally Transformed Linear Mixed Models (FaST-LMM). FaST-LMM, which stands for Factored Spectrally Transformed Linear Mixed Models is a program for performing single-SNP and SNP-set genome-wide association studies (GWAS) on extremely large data sets. It runs on both Windows and Linux systems, and has been tested on data sets with over 120,000 individuals.


# How to Install

To install this package, use devtools:

devtools::install_github("camult/easyGWAS")


# Overview

To see the full list of exported functions:

ls(package:easyGWAS)

# `fastlmmGWAS`: Factored Spectrally Transformed Linear Mixed Model GWAS

## Description


 R interface to perform GWAS using Factored Spectrally Transformed Linear Mixed Models (FaST-LMM).
 FaST-LMM, which stands for Factored Spectrally Transformed Linear Mixed Models is a
 program for performing single-SNP and SNP-set genome-wide association studies
 (GWAS) on extremely large data sets. It runs on both Windows and Linux systems,
 and has been tested on data sets with over 120,000 individuals.


## Usage

```r
fastlmmGWAS(formula = NULL, genoFileName, phenFileName, MarkerType, IDname,
  mapFileName = NULL, nPC = 0, useG = FALSE, maf = 0.01,
  covariate = NULL, rmNonP = TRUE, rmMAF = TRUE, rSrcDir = NULL,
  phenName = NULL, MarkerRow = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
```formula```     |     A formula specifying the model.
```genoFileName```     |     It is the name of markers file with its extension, i.e., "datafile.txt"
```phenFileName```     |     It is the name of phenotypes file with its extension, i.e., "phenfile.txt"
```MarkerType```     |     Use "SNP" or "AFLP".
```IDname```     |     The name of the trait.
```mapFileName```     |     It is the name of map file with its extension, i.e., "mapfile.txt"
```nPC```     |     Number of principal components if any.
```useG```     |     A logical value indicating whether the kinship must be used.
```maf```     |     A value indicating the minor allele frequency.
```covariate```     |     Name of covariate(s) if there is any.
```rmNonP```     |     A logical value indicating whether the non-polimorfic markers must be removed.
```rmMAF```     |     A logical value indicating whether the non-polimorfic markers must be removed.
```rSrcDir```     |     Optional path to the folder where the FaST-LMM programm are.
```phenName```     |     It is the name of the phenotype to be analyzed.
```MarkerRow```     |     A logical value indicating whether the markers are in rows.

## Value


 A text file with GWAS statistics. Manhattan plot and QQ-Plot.


## References


 C. Lippert, J. Listgarten, Y. Liu, C.M. Kadie, R.I. Davidson, and D. Heckerman.
 FaST Linear Mixed Models for Genome-Wide Association Studies.
 Nature Methods 8: 833-835, Oct 2011 (doi:10.1038/nmeth.1681).
 
 C. Widmer, C. Lippert, O. Weissbrod, N. Fusi, C.M. Kadie, R.I. Davidson, J.
 Listgarten, and D. Heckerman. Further Improvements to Linear Mixed Models for
 Genome-Wide Association Studies.
 Scientific Reports 4, 6874, Nov 2014 (doi:10.1038/srep06874).


