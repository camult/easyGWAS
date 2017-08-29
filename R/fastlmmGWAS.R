#' @title Factored Spectrally Transformed Linear Mixed Model GWAS
#' 
#' @description R interface to perform GWAS using Factored Spectrally Transformed Linear Mixed Models (FaST-LMM). 
#'              FaST-LMM, which stands for Factored Spectrally Transformed Linear Mixed Models is a 
#'              program for performing single-SNP and SNP-set genome-wide association studies 
#'              (GWAS) on extremely large data sets. It runs on both Windows and Linux systems, 
#'              and has been tested on data sets with over 120,000 individuals.
#' 
#' @param formula A formula specifying the model.
#' @param genoFileName It is the name of markers file with its extension, i.e., "datafile.txt"
#' @param phenFileName It is the name of phenotypes file with its extension, i.e., "phenfile.txt"
#' @param mapFileName It is the name of map file with its extension, i.e., "mapfile.txt"
#' @param nPC Number of principal components if any.
#' @param rmNonP A logical value indicating whether the non-polimorfic markers must be removed.
#' @param rmMAF A logical value indicating whether the non-polimorfic markers must be removed. 
#' @param useG A logical value indicating whether the kinship must be used. 
#' @param maf A value indicating the minor allele frequency.
#' @param phenName It is the name of the phenotype to be analyzed.
#' @param rSrcDir Optional path to the folder where the FaST-LMM programm are. 
#' @param IDname The name of the trait.
#' @param MarkerRow A logical value indicating whether the markers are in rows. 
#' @param covariate Name of covariate(s) if there is any.
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
#' @import data.table lattice qqman
fastlmmGWAS <- function(formula=NULL, genoFileName, phenFileName, IDname,
                        mapFileName = NULL, nPC=0, useG=FALSE, maf=0.01, covariate=NULL,
                        rmNonP=TRUE, rmMAF=TRUE, rSrcDir=NULL, phenName=NULL, MarkerRow=TRUE){
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
  if(!is.null(formula)){
    form <- formula
    phenName <- as.character(form[[2]])
  }else{
    phenName <- phenName
    form <- formula(paste(phenName," ~ 1"))
  }
  rmNonP <- rmNonP
  rmMAF <- rmMAF
  useG <- useG
  data_file <- normalizePath(genoFileName)
  phenotype_file <- normalizePath(phenFileName)
  map_file <- mapFileName
  if (length(map_file) != 0){
    map_file = normalizePath(map_file)	
  }
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
  phenotypes = data.frame(fread(phenotype_file, head=TRUE))
  # Genotype file
  dataWithIDs = as.matrix(fread(data_file, head=TRUE, showProgress=F))
  if(!isTRUE(MarkerRow)){
    IDs <- dataWithIDs[,1]
    dataWithIDs <- data.frame(Marker=colnames(dataWithIDs)[-1], t(dataWithIDs[,-1]))
    colnames(dataWithIDs) <- c("Marker", IDs)
  }
  # SampleIDs has been changed to unique IDs even if they are repeated
  if(nrow(phenotypes)>length(unique(phenotypes[,IDname]))){
    repeated <- function(phenData, genoData){
      #nRep <- nrow(phenData)/length(unique(phenData[,1]))
      #Rep <- sort(rep(1:nRep, length(unique(phenotypes[,1]))))
      #phenData$Rep <- Rep
      #formula <- update(form, ~ . + Rep)
      repeatedIDs <- phenData[,IDname]
      sampleIDs <- unique(phenData[,IDname])
      markerIDs <- genoData[,1]
      nRepeated <- length(repeatedIDs)/length(sampleIDs)
      fakeIDs   <- paste("fakeIDs", 1:length(repeatedIDs), sep="_")
      auxMarkers <- genoData[,-1]
      repeatMarkers = sapply(repeatedIDs, function(id) c(auxMarkers[,colnames(auxMarkers)%in%id]) )
      colnames(repeatMarkers) <- fakeIDs
      phenData[,IDname] <- fakeIDs
      repeatMarkers <- cbind(markerIDs, repeatMarkers)
      return(list(phenData=phenData, genoData=repeatMarkers))
    }
    repIDs <- repeated(phenData=phenotypes, genoData=dataWithIDs)
    phenotypes <- repIDs$phenData
    dataWithIDs <- repIDs$genoData
    rm(repIDs)
  }
  sampleIDs <- phenotypes[,IDname]
  phenoDict = as.matrix(phenotypes[,phenName])
  mode(phenoDict) <- "numeric"
  colnames(phenoDict) <- phenName
  rownames(phenoDict) = phenotypes[,IDname]
  rownames(dataWithIDs) <- dataWithIDs[, 1]
  data = as.matrix(dataWithIDs[, -1])
  rownames(data) <- rownames(dataWithIDs)
  # create phenotype, fixed effects and fam file
  outPheno = data.frame(Family=1,
                        ID=rownames(phenoDict),
                        phenoDict)
  write.table(outPheno, file = 'phenotype.txt', quote=F, sep = '\t', row.names=F, col.names=F)
  # Changed: 05-15-2017
  AllEff <- unlist(strsplit(as.character(formula[[3]])," "))[!unlist(strsplit(as.character(formula[[3]])," "))%in%c("+")]
  CatEff <- AllEff[!AllEff%in%covariate]
  phenotypes[CatEff] <- lapply(phenotypes[CatEff] , factor)
  #--------------------#
  if(!is.null(formula)){
    effect <- model.matrix(update(formula, ~ . -1), phenotypes)
    FixedEffects = data.frame(Family=1,
                              ID=rownames(phenoDict),
                              effect)
    FixedEffects <- FixedEffects[, !colnames(FixedEffects)%in%phenName]
    write.table(FixedEffects, file = 'FixedEffects.txt', quote=F, sep = '\t', row.names=F, col.names=T)
  }
  outFam = data.frame(Family=1,
                      Individual=rownames(phenoDict),
                      Father=0,
                      Mother=0,
                      Sex=0,
                      Phenotype=-9)
  write.table(outFam, file = 'testmarkers.fam', quote=F, sep = '\t', row.names=F, col.names=F)
  # Quality Control Filter
  # Non-polymorphic Data
  if(rmNonP) NonP <- apply(data, 1, sd, na.rm=T)
  # Minor allele frequency
  if(rmMAF){
    n0 <- apply(data==0,1,sum,na.rm=T)
    n1 <- apply(data==1,1,sum,na.rm=T)
    n2 <- apply(data==2,1,sum,na.rm=T)
    p  <- ((2*n0)+n1)/(2*(n0+n1+n2))
    mafOut <- as.matrix(pmin(p, 1-p))
  }
  # Apply quality control filter if any
  if(all(rmNonP, rmMAF)){
    data <- data[(NonP!=0 & mafOut>=maf), ]
  } 
  if(isTRUE(rmNonP)&!isTRUE(rmMAF)){
    data <- data[(NonP!=0), ]
  }
  if(!isTRUE(rmNonP)&isTRUE(rmMAF)){
    data <- data[(mafOut>=maf), ]
  }
  # Map and genotype files
  markerIDs <- rownames(data)
  nSNPs <- length(markerIDs)
  ReadMapFile <- function(map_file, rSrcDir){
    if (length(map_file) == 0){
      cat("  creating a fake map file\n")
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
      chrDict        <- chr
      names(chrDict) <- rownames(data)
      posDict        <- as.integer((1:nrow(data))*100000)
      names(posDict) <- rownames(data)
    } else {
      cat("  reading in map file:", map_file, '\n')
      map = read.table(map_file, FALSE)
      chrDict = map[,2]
      names(chrDict) = map[,1]
      posDict = map[,3]
      names(posDict) = map[,2]
    }
    return(list(ChrDict = chrDict, PosDict = posDict))
  }
  map = ReadMapFile(map_file, rSrcDir)
  # remove rows with probe ID not found in map file
  rowsToKeep = markerIDs %in% names(map$ChrDict)
  numFound = sum(rowsToKeep)
  numMissing = sum(!rowsToKeep)
  markerIDs = markerIDs[rowsToKeep]
  data = data[rowsToKeep, ]
  mode(data) <- "numeric"
  cat("  number of missing map markers =", numMissing, '\n')
  cat("  number of markers =", numFound, '\n')
  if (numFound == 0) stop("no mapping information found for any markers")
  # get genotype file contents
  dataDict = as.data.frame(apply(data, 1, function(row){
    return(append(c('A', 'B'), row))
  }))
  names(dataDict) = markerIDs
  rownames(dataDict) = append(c('A', 'B'), sampleIDs)
  # get sorted map file contents
  rowsToKeep = names(map$ChrDict) %in% markerIDs
  chr = map$ChrDict[rowsToKeep]
  nChr = length(unique(chr))
  pos = map$PosDict[rowsToKeep]
  tempMap = data.frame(Chr = chr, markerID = names(chr), PosA = pos, PosB = pos)
  sortedMap = suppressWarnings(tempMap[order(as.numeric(as.vector(tempMap$Chr)),
                                             tempMap$Chr, tempMap$PosA, tempMap$markerID), ])
  # create map file and data file
  write.table(sortedMap, file = 'testmarkers.map', quote=F, sep = '\t', row.names=F, col.names=F)
  outData = t(apply(sortedMap, 1, function(row) as.vector(dataDict[, row[2]])))
  write.table(outData, file = 'testmarkers.dat', quote=F, sep = '\t', row.names=T, col.names=F)
  # create principal components file
  if(nPC>0){
    cat('Computing principal components...\n')
    pcData = outData[, 3:ncol(outData)]
    mode(pcData) = 'numeric'
    dataM = apply(pcData, 1, function(row){
      average = mean(row, na.rm = TRUE)
      return(row - average)
    })
    cat("  number of PC used =", nPC, '\n')
    X = svd(dataM) #  singular value decomposition
    scores = t(apply(X$u, 1, function(row) row * sqrt(X$d)))
    write.table(scores[, 1:100], file = 'pc.txt', quote=F, sep = '\t', row.names=F, col.names=F)
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
  ToRemove <- c("covariates_pc.txt","covariates.txt","testmarkers.map","testmarkers.dat","testmarkers.fam",
                "phenotype.txt","FixedEffects.txt","pc.txt","desingMatrix.txt","desingMatrixPC.txt")
  suppressWarnings(invisible(file.remove(ToRemove, resultsNames)))
  cat("\n")
}
