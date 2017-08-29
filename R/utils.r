#' @title Linear Model for Genome-wide Association Studies
#' 
#' @description The lmGWAS (Linear Model for Genome-wide Association Studies) is developed 
#'              as a package for GWAS based on a single SNP analysis .
#' 
#' @param nIDs Number of individuals.
#' @param nSNPs Number of markers.
#' @param h2 Heritability.
#' @param nQTL Number of QTLs.
#' @param CatEff A logical value indicating whether a categorial variable must be used. 
#' @param nCat Number of distinct categories in CatEff.
#' @param CovEff A logical value indicating whether a covariate variable must be used. 
#' 
#' @return A text file with GWAS statistics. Manhattan plot and QQ-Plot. If an interaction 
#'         between the SNP and any other fixed effect is in the model, its pvalues and 
#'         some statistics are also saved in the text file.
#'             
#' @export gpSIM
#' 
gpSIM <- function(nIDs, nSNPs, h2, nQTL, CatEff=TRUE, nCat=5, CovEff=TRUE){
  geno <- matrix(sample(c(0,1,2), nIDs*nSNPs, replace=TRUE), nrow = nSNPs)
  X=scale(t(geno))/sqrt(nSNPs)
  QTL=seq(from=1, to=nSNPs,length=nQTL)
  b <- rep(1,nQTL)
  signal=X[,QTL]%*%b
  signal=signal/sd(signal)*sqrt(h2)
  error=rnorm(n=nIDs)
  error=error/sd(error)*sqrt(1-h2)
  if(isTRUE(CatEff)){
    CatEff <- factor(sample(LETTERS[1:nCat], nIDs, replace=TRUE))
  }else{
    CatEff <- 0
  }
  if(isTRUE(CovEff)){
    CovEff <- round(rnorm(n=nIDs, mean = 6, sd = 2), 2)
  }else{
    CovEff <- 0
  }
  y=round(signal+error+as.numeric(CatEff)+CovEff, 2)
  phen <- data.frame(ID=paste0("ID",1:nIDs), Trait=y, CatEff=CatEff, CovEff=CovEff)
  colnames(geno) <- phen$ID
  geno <- data.frame(Marker=paste0("Marker",1:nSNPs), geno)
  chr=rep(c(1:4), (nSNPs/4))
  nChr <- length(chr)
  if(nChr!=nSNPs){
    if(nChr>nSNPs){
      dif <- abs(nSNPs-nChr)+1
      chr <- chr[-c(1:dif)]
    }else{
      dif <- abs(nSNPs-nChr)
      chr <- c(rep(1, dif), chr)
    }
  }
  map <- data.frame(chr=chr,
                    Marker=paste0("Marker",1:nSNPs),
                    PosA=as.integer((1:nSNPs)*1000),
                    PosB=as.integer((1:nSNPs)*1000))
  return(list(phen=phen, map=map, geno=geno))
}

desingMatrix <- function(formula){
  outFEffe <- read.table("FixedEffects.txt", header = T, sep = "\t")
  desingMatrix = paste('desingMatrix.txt', sep = '')
  write.table(outFEffe, file = desingMatrix, quote=F, sep='\t',row.names=F, col.names=F)
  return(desingMatrix)
}

desingMatrixPC <- function(formula, nPC){
  outFEffe <- read.table("FixedEffects.txt", header = T, sep = "\t")
  pc = read.delim('pc.txt', FALSE, sep = "\t")
  pcCovar = cbind(outFEffe, pc[, 1:nPC])
  desingMatrixPC = paste('desingMatrixPC.txt', sep = '')
  write.table(pcCovar, file = desingMatrixPC, quote=F, sep='\t',row.names=F, col.names=F)
  return(desingMatrixPC)
}

WritePCCovar <- function(nPC){
  CreateCovarFile()
	covar = read.delim("covariates.txt", FALSE)
	pc = read.delim('pc.txt', FALSE, sep = "\t")
  pcCovar = cbind(covar, pc[, 1:nPC])
  pcCovarFileName = paste('covariates_pc.txt', sep = '')
	write.table(pcCovar, file = pcCovarFileName, quote=F, sep='\t',row.names=F, col.names=F)
	return(pcCovarFileName)
}

CreateCovarFile <- function(){
  covar = read.delim('testmarkers.fam', FALSE)
  covarOut ='covariates.txt'
  write.table(covar[, 1:2], file=covarOut, quote=F, sep='\t', row.names=F, col.names=F)
  return(covarOut)
}

RunMainLoop <- function(formula=formula, fastlmmFileName, nPC, useG, nChr){
  #---------------------------------------------------------------------------------------#
  # Running GWAS analyses
  #---------------------------------------------------------------------------------------#
  cat("  performing GWAS...\n")
  if(nPC==0){
    if(all(isTRUE(useG) & is.null(formula))){
      covarFileName <- CreateCovarFile()
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -dfile1Sim testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -simLearnType Full -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", covarFileName, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(all(isTRUE(useG) & !is.null(formula))){
      desingMatrixCOV = desingMatrix(formula)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -dfile1Sim testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -simLearnType Full -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
        } else {
          Sys.setenv(FastLmmUseAnyMklLib=1)
          system(cmd)
      }
    }
    if(all(!isTRUE(useG) & is.null(formula))){
      covarFileName <- CreateCovarFile()
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", covarFileName, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(all(!isTRUE(useG) & !is.null(formula))){
      desingMatrixCOV = desingMatrix(formula)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
  }
  
  if(nPC>0){
    if(all(isTRUE(useG) & is.null(formula))){
      pcCovarFileName <-  WritePCCovar(nPC)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -dfile1Sim testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -simLearnType Full -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", pcCovarFileName, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(all(isTRUE(useG) & !is.null(formula))){
      desingMatrixCOV <- desingMatrixPC(formula, nPC)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -dfile1Sim testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -simLearnType Full -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(all(!isTRUE(useG) & is.null(formula))){
      pcCovarFileName <- WritePCCovar(nPC)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", pcCovarFileName, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(all(!isTRUE(useG) & !is.null(formula))){
      desingMatrixCOV <- desingMatrixPC(formula, nPC)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -dfile1 testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
  }
  gwas <- fread("gwas.txt", header=T, data.table = F)
  gwas[is.na(gwas)] <- 1
  if("SnpWeight"%in%colnames(gwas)){
    gwas$SNPeff <- as.numeric(abs(gwas$SnpWeight))
  }
  if("SNPWeight"%in%colnames(gwas)){
    gwas$SNPeff <- as.numeric(abs(gwas$SNPWeight))
  }
  gwas$Pvalue <- as.numeric(gwas$Pvalue)
  trellis.device(device="png",
                 filename=paste0("QQPlot-GWAS.png"),
                 width=617,
                 height=397)
  qqman::qq(as.numeric(gwas$Pvalue))
  dev.off()
  
  trellis.device(device="png",
                 filename=paste0("Pvalue-GWAS.png"),
                 width=617,
                 height=397)
  qqman::manhattan(gwas, chr="Chromosome", bp="Position", p="Pvalue", snp="SNP", 
                   cex=0.5, cex.axis=0.7, col=c("blue4","orange2"), 
                   suggestiveline=-log10(0.05/nrow(gwas)), genomewideline=F, logp=T)
  dev.off()
  
  trellis.device(device="png",
                 filename=paste0("SNPeff-GWAS.png"),
                 width=617,
                 height=397)
  qqman::manhattan(gwas, chr="Chromosome", bp="Position", p="SNPeff", snp="SNP", 
                   cex=0.5, cex.axis=0.7, col=c("blue4","orange2"), 
                   ylim=c(0, max(gwas$SNPeff)+0.05*max(gwas$SNPeff)),
                   ylab="SNP effect",
                   suggestiveline=F, genomewideline=F, logp=F)
  dev.off()
}
