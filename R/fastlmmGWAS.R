#' @title Factored Spectrally Transformed Linear Mixed Model GWAS
#' 
#' @description R interface to perform GWAS using Factored Spectrally Transformed Linear Mixed Models (FaST-LMM). 
#'              FaST-LMM, which stands for Factored Spectrally Transformed Linear Mixed Models is a 
#'              program for performing single-SNP and SNP-set genome-wide association studies 
#'              (GWAS) on extremely large data sets. It runs on both Windows and Linux systems, 
#'              and has been tested on data sets with over 120,000 individuals.
#' 
#' @param formula A formula specifying the model.
#' @param geno It is the name of markers file.
#' @param phen It is the name of phenotypes file.
#' @param map It is the map file with colukns in the following order: Marker, Chromosome and Position.
#' @param IDname It is the individual's name.
#' @param nPC Number of principal components if any.
#' @param useG A logical value indicating whether the kinship must be used. 
#' @param MAF A value indicating the Minor Allele Frequency.
#' @param HWE A value indicating the Hardy-Weinberg Equilibrium.
#' @param HZ A value indicating the SNP heterozygosity.
#' @param SNPcall A value indicating the SNP call rate.
#' @param INDcall A value indicating the  individual call rate for autosomal SNPs.
#' @param rmMAF A logical value indicating whether markers sould be removed based on MAF. 
#' @param rmHWE A logical value indicating whether markers sould be removed based on HWE.
#' @param rmHZ A logical value indicating whether markers sould be removed based on HZ. 
#' @param rmSNPCall A logical value indicating markers sould be removed based on SNP call rate. 
#' @param rmINDCall A logical value indicating markers sould be removed based on individual call rate. 
#' @param rSrcDir Optional path to the folder where the FaST-LMM programm are. 
#' @param MarkerRow A logical value indicating whether the markers are in rows. 
#' 
#' @return A text file with GWAS statistics. Manhattan plot and QQ-Plot. 
#' 
#' @references C. Lippert, J. Listgarten, Y. Liu, C.M. Kadie, R.I. Davidson, and D. Heckerman. 
#'             FaST Linear Mixed Models for Genome-Wide Association Studies. 
#'             Nature Methods 8: 833-835, Oct 2011 (doi:10.1038/nmeth.1681).
#'             
#'             C. Widmer, C. Lippert, O. Weissbrod, N. Fusi, C.M. Kadie, R.I. Davidson, J. 
#'             Listgarten, and D. Heckerman. Further Improvements to Linear Mixed Models for 
#'             Genome-Wide Association Studies. 
#'             Scientific Reports 4, 6874, Nov 2014 (doi:10.1038/srep06874).
#'             
#' @export fastlmmGWAS
#' @import data.table lattice qqman gaston
fastlmmGWAS <- function(formula=NULL, geno, phen, IDname, map, nPC=0, useG=FALSE,
                        MAF=0.01, HWE=1e-10, HZ=0.01, SNPcall=0.90, INDcall=0.90, 
                        rmMAF=TRUE, rmHWE=TRUE, rmHZ=TRUE, rmSNPCall=TRUE, rmINDCall=FALSE,
                        rSrcDir=NULL, phenName=NULL, MarkerRow=TRUE){
  cat("\n")
  centerText <- function() {
    width <- getOption("width")
    A <- ("                                                       ._____.    \n")
    B <- ("    _/////_            Fernando Brito Lopes           _|_____|_   \n")
    C <- ("   (' o o ')     Animal Scientist (Zootechnician)     (' o o ')   \n")
    D <- ("__ooO_(_)_Ooo_____ Animal Breeding and Genetics _____ooO_(_)_Ooo__\n")
    E <- ("                    e-mail: <camult@gmail.com>                    \n")
    ws <- rep(" ", floor((width - nchar(A))/2))
    cat(ws,A,ws,B,ws,C,ws,D,ws,E,ws,sep = "")
  }
  centerText()
  cat("\n")
  if (substr(version$platform, 1, 11) == 'x86_64-w64-'){
    if(is.null(rSrcDir)){
      if(dir.exists(paste0(unique(c(.Library.site, .Library))[1],"/easyGWAS/bin"))){
        rSrcDir <- paste0(unique(c(.Library.site, .Library))[1],"/easyGWAS/bin")
      }
    }
  } else if (version$os == "linux" | version$os == "linux-gnu"){
    if(is.null(rSrcDir)){
      if(dir.exists(paste0(paste0(.libPaths(),"/easyGWAS")[paste0(.libPaths(),"/easyGWAS")==find.package("easyGWAS")],"/bin"))){
        rSrcDir <-  paste0(paste0(.libPaths(),"/easyGWAS")[paste0(.libPaths(),"/easyGWAS")==find.package("easyGWAS")],"/bin")
      }
    }
  } else {
    stop('Unsupported architecture or operating system')
  }
  model <- formula
  phenName  <- unlist(strsplit(as.character(model[[2]])," "))[!unlist(strsplit(as.character(model[[2]])," "))%in%c("|")]
  effect <- unlist(strsplit(as.character(model[[3]])," "))[!unlist(strsplit(as.character(model[[3]])," "))%in%c("+")]
  if("1"%in%effect){
    file <- phen[, c(IDname, phenName)]
    colnames(file) <- c(IDname, phenName)
  }else{
    file <- phen[, c(IDname, phenName, effect)]
    colnames(file) <- c(IDname, phenName, effect)
  }
  file <- file[, !duplicated(colnames(file))]
  # get the filename of the executable for this OS
  if (substr(version$platform, 1, 11) == 'x86_64-w64-'){
    fastlmmFileName = file.path(rSrcDir, 'fastlmmC.exe')
  } else if (version$os == "linux" | version$os == "linux-gnu"){
    fastlmmFileName = normalizePath(file.path(rSrcDir, 'fastlmmC'))
    system(paste0("chmod 777 ", fastlmmFileName))
  } else {
    stop('Unsupported architecture or operating system')
  }
  cat('Preparing inputs for FaST-LMM-GWAS...\n')
  # Phenotype file
  phenotypes = data.frame(file)
  # Genotype file
  dataWithIDs = geno
  if(!isTRUE(MarkerRow)){
    IDs <- dataWithIDs[,1]
    dataWithIDs <- data.frame(Marker=colnames(dataWithIDs)[-1], t(dataWithIDs[,-1]))
    colnames(dataWithIDs) <- c("Marker", IDs)
  }
  # SampleIDs has been changed to unique IDs even if they are repeated
  sampleIDs <- phenotypes[,IDname]
  phenoDict = as.matrix(phenotypes[,phenName])
  mode(phenoDict) <- "numeric"
  colnames(phenoDict) <- phenName
  rownames(phenoDict) = phenotypes[,IDname]
  data = as.matrix(dataWithIDs)
  rownames(data) <- rownames(dataWithIDs)
  map <- data.frame(map)
  # create phenotype, fixed effects and fam file
  outPheno = data.frame(Family=1,
                        ID=rownames(phenoDict),
                        phenoDict)
  write.table(outPheno, file = 'phenotype.txt', quote=F, sep = '\t', row.names=F, col.names=F)
  effect <- model.matrix(model, phenotypes)
  FixedEffects = data.frame(Family=1,
                            ID=rownames(phenoDict),
                            effect)
  FixedEffects <- FixedEffects[, !colnames(FixedEffects)%in%phenName]
  write.table(FixedEffects, file = 'FixedEffects.txt', quote=F, sep = '\t', row.names=F, col.names=T)
  PLINK <- as.bed.matrix(t(data),
                         fam=data.frame(famid=1,
                                        id=colnames(data),
                                        father=0,
                                        mother=0,
                                        sex=0,
                                        pheno=-9),
                         bim=data.frame(chr=map[,2],
                                        id=map[,1],
                                        dist=0,
                                        pos=map[,3],
                                        A1="A",
                                        A2="B"))
  # Quality Control Filter
  if(rmMAF) suppressWarnings(PLINK <- select.snps(PLINK, maf>MAF))
  if(rmSNPCall) suppressWarnings(PLINK <- select.snps(PLINK, callrate>SNPcall))
  if(rmINDCall) suppressWarnings(PLINK <- select.inds(PLINK, callrate>INDcall))
  if(rmHZ) suppressWarnings(PLINK <- select.snps(PLINK, hz>HZ))
  if(rmHWE) {
    PLINK <- set.hwe(PLINK)
    suppressWarnings(PLINK <- select.snps(PLINK, hwe>HWE))
  }
  nChr <- length(table(PLINK@snps$chr))
  # create map file and data file
  write.bed.matrix(PLINK, "testmarkers", rds=NULL)
  # create principal components file
  if(nPC>0){
    # Compute Genetic Relationship Matrix
    standardize(PLINK) <- "p"
    K <- GRM(PLINK)
    # Eigen decomposition
    eiK <- eigen(K)
    # deal with small negative eigen values
    eiK$values[ eiK$values < 0 ] <- 0
    # Note: the eigenvectors are normalized, to compute 'true' PCs
    # multiply them by the square root of the associated eigenvalues
    cat('Computing principal components using the genomic relationship matrix...\n')
    PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
    write.table(PC, file = 'pc.txt', quote=F, sep = '\t', row.names=F, col.names=F)
  }
  cat('Running FaST-LMM-GWAS...\n')
  TK = RunMainLoop(formula=formula, fastlmmFileName=fastlmmFileName, nPC=nPC, useG=useG, nChr=nChr)
  # Cleaning
  resultsNames <- c("Pvalue-GWAS.png","QQPlot-GWAS.png","SNPeff-GWAS.png","gwas.txt","log.txt")
  resultGWAS <- paste0("GWAS_", phenName)
  if(!dir.exists(resultGWAS)){
    dir.create(resultGWAS)
  }
  CurrentDir <- getwd()
  FutureDir  <- paste(CurrentDir, resultGWAS, sep="/")
  CopyFile <- function(file){
    FileListSource<- paste(CurrentDir, file, sep="/")
    FileListEnd <- paste(FutureDir, file, sep="/")
    file.copy(from=FileListSource, to=FileListEnd)
  }
  invisible(sapply(resultsNames, CopyFile))
  ToRemove <- c("covariates_pc.txt","covariates.txt","testmarkers.bim","testmarkers.bed","testmarkers.fam",
                "phenotype.txt","FixedEffects.txt","pc.txt","desingMatrix.txt","desingMatrixPC.txt")
  suppressWarnings(invisible(file.remove(ToRemove, resultsNames)))
  cat("\n")
}
