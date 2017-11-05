desingMatrix <- function(formula){
  outFEffe <- fread("FixedEffects.txt", header = T, sep = "\t", data.table=F)
  desingMatrix = paste('desingMatrix.txt', sep = '')
  write.table(outFEffe, file = desingMatrix, quote=F, sep='\t',row.names=F, col.names=F)
  return(desingMatrix)
}

desingMatrixPC <- function(formula, nPC){
  outFEffe <- fread("FixedEffects.txt", header = T, sep = "\t", data.table=F)
  pc = fread('pc.txt', header = FALSE, sep = "\t", data.table=F)
  pcCovar = cbind(outFEffe, pc[, 1:nPC])
  desingMatrixPC = paste('desingMatrixPC.txt', sep = '')
  write.table(pcCovar, file = desingMatrixPC, quote=F, sep='\t',row.names=F, col.names=F)
  return(desingMatrixPC)
}

#---------------------------------------------------------------------------------------#
# Running GWAS analyses
#---------------------------------------------------------------------------------------#
RunMainLoop <- function(formula=formula, fastlmmFileName, nPC, useG, nChr){
  cat("  performing GWAS...\n")
  if(nPC==0){
    if(isTRUE(useG)){
      desingMatrixCOV = desingMatrix(formula)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -bfile testmarkers -bfileSim testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -simLearnType Full -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(!isTRUE(useG)){
      desingMatrixCOV = desingMatrix(formula)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -bfile testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
  } else {
    if(isTRUE(useG)){
      desingMatrixCOV <- desingMatrixPC(formula, nPC)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -bfile testmarkers -bfileSim testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -simLearnType Full -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
      if (substr(version$platform, 1, 11) == 'x86_64-w64-') {
        shell(cmd) 
      } else {
        Sys.setenv(FastLmmUseAnyMklLib=1)
        system(cmd)
      }
    }
    if(!isTRUE(useG)){
      desingMatrixCOV <- desingMatrixPC(formula, nPC)
      cmd = paste0(fastlmmFileName, " -mpheno 1 -bfile testmarkers -pheno phenotype.txt -SetOutputPrecision 3 -linreg -out gwas.txt -REML -MaxChromosomeValue ", nChr," -covar ", desingMatrixCOV, " 2>log.txt")
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
