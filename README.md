# `fastlmmGWAS`: Factored Spectrally Transformed Linear Mixed Model GWAS

## Description


 R interface to perform GWAS using Factored Spectrally Transformed Linear Mixed Models (FaST-LMM).
 FaST-LMM, which stands for Factored Spectrally Transformed Linear Mixed Models is a
 program for performing single-SNP and SNP-set genome-wide association studies
 (GWAS) on extremely large data sets. It runs on both Windows and Linux systems,
 and has been tested on data sets with over 120,000 individuals.


## Usage

```r
fastlmmGWAS(formula = NULL, geno, phen, IDname, map, nPC = 0,
  useG = FALSE, MAF = 0.01, HWE = 1e-10, HZ = 0.01, SNPcall = 0.9,
  INDcall = 0.9, rmMAF = TRUE, rmHWE = TRUE, rmHZ = TRUE,
  rmSNPCall = TRUE, rmINDCall = FALSE, rSrcDir = NULL, phenName = NULL,
  MarkerRow = TRUE)
```


## Arguments

Argument      |Description
------------- |----------------
```formula```     |     A formula specifying the model.
```geno```     |     It is the name of markers file.
```phen```     |     It is the name of phenotypes file.
```IDname```     |     It is the individual's name.
```map```     |     It is the map file with colukns in the following order: Marker, Chromosome and Position.
```nPC```     |     Number of principal components if any.
```useG```     |     A logical value indicating whether the kinship must be used.
```MAF```     |     A value indicating the Minor Allele Frequency.
```HWE```     |     A value indicating the Hardy-Weinberg Equilibrium.
```HZ```     |     A value indicating the SNP heterozygosity.
```SNPcall```     |     A value indicating the SNP call rate.
```INDcall```     |     A value indicating the  individual call rate for autosomal SNPs.
```rmMAF```     |     A logical value indicating whether markers sould be removed based on MAF.
```rmHWE```     |     A logical value indicating whether markers sould be removed based on HWE.
```rmHZ```     |     A logical value indicating whether markers sould be removed based on HZ.
```rmSNPCall```     |     A logical value indicating markers sould be removed based on SNP call rate.
```rmINDCall```     |     A logical value indicating markers sould be removed based on individual call rate.
```rSrcDir```     |     Optional path to the folder where the FaST-LMM programm are.
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


