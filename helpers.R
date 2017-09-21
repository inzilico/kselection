# Functions used in make_test_files.R  
# Author: Gennady Khvorykh, http://inZilico.com

library(plyr, quietly = T)

ReadFastPHASE <- function(input, ...){
  # Reads fastPHASE input/output files
  # Args:
  #   Input: full/path/to/filename.inp
  # Returns:
  #   A character vector with sequences 

  # Check input
  if(file.access(input) == -1) stop((sprintf("File %s doesn't exist", input)))
  
  # Load data set
  l <- readLines(input, ...)
  
  # Remove lines with numbers
  l <- l[grep("\\d+", l, invert = T)]
  
  # Remove lines having `#`  
  l <- l[grep("#", l, invert = T)]
  
  # Print info about dataset
  message(sprintf("%s sequences and %s markers are loaded from %s.", 
                  length(l), nchar(l[1]), input))
  l
}

MaskSequence <- function(sequence, positions, symbol = "?"){
  # Substitutes base in a sequence for `mask` symbol
  # Input: 
  #   sequence: character string
  #   positions: vector with positions of characters to be replaced by symbol
  #   symbol: character representing missing value
  # Returns:
  #   Sequence where some bases are replaced by symbol 
  
  # Loop throug all positions to be replaced
  for(x in positions) substr(sequence, start = x, stop = x) <- symbol 
  sequence
}

ApplyMasks <- function(g, masks, pref) {
  # Applies set of masks to sequencies. Saves the result as fastPHASE input. 
  # Args:
  #   g: character vector with sequences
  #   masks: list of masks as binary matrices 
  #   pref: prefix and path for output files
  # Return: 
  #   No values
  
  # Initilize 
  N <- length(g) # Number of sequences
  M <- nchar(g[1]) # Number of markers
  
  # Loop throug all masks 
  for (n in seq_along(masks)) {
    
    # Get indexes of masked genotypes
    ind <- alply(masks[[n]], 1, function(v) which(v == 1))
    
    # Replicate indexes
    ind <- rep(ind, each = 2)
   
    message(sprintf("Applying mask %s...", n))
    
    # Lopp throug all sequences and mask them
    gm <- llply(seq_len(N), .progress = create_progress_bar(name = "text"), 
                function(i, g, ind) MaskSequence(g[i], ind[[i]]), g = g, ind = ind)
    
    # Set output filename
    fn <- sprintf("%s.m%s.inp", pref, n)
    
    # Remove output file if it exists
    if(file.exists(fn)) file.remove(fn)
    
    # Create output file
    write(c(N/2, M), file = fn, ncolumns = 1)
    
    # Write masked sequences to file
    l_ply(gm, write, file = fn, append = T)
    
    message(sprintf("File %s is saved", fn))
  }

}

GetMissingMarkers <- function(data){
  # Determines the positions of missing bases per individual
  # Args:
  #   data: character vector of sequnces
  # Returns:
  #   List with positions of missed genotypes for each individual
  mindex <- sapply(data, function(x) unlist(gregexpr("\\?", x)), USE.NAMES = F)
  mindex[seq(2, length(mindex), 2)]
}

ToBinary <- function(data, jmax){
  # Converts list with the positions into binary matrix
  # Args:
  #  data: list with the positions
  #  jmax: number of markers
  # Returns:
  #  Binary matrix, where 1 means missing and 0 - nonmissing values. 
  sapply(seq_len(jmax), function(j){
    sapply(data, function(i) ifelse(j %in% i, 1, 0))
  })
}

GetMissing <- function(data, M){
  # Determines missing genotypes
  # Args:
  #   data: character vector with sequences
  #   M: number of markers
  # Returns:
  #   Binary matrix, where 1 corresponds to missing genotype.
  
  message("Counting missing genotypes...")
  
  # Get list of missing markers per individual
  tmp <- GetMissingMarkers(data)
  
  # Convert list of missing markers into binary matrix 
  m <- ToBinary(tmp, M)
  
  # Print proportion of missing genotypes
  p <- sum(colSums(m))/(dim(m)[1] * dim(m)[2])
  message("Proportion of missing genotype: ", round(p, 4))
  
  m
}

GenerateMask <- function(m, masked, size){
  # Generate a new mask
  # Args:
  #   m: matrix with zeros
  #   masked: matrix keeping genotypes previously masked
  #   size: vector with the number of genotypes available for masking 
  # Returns:
  #   A new mask as a binary matrix
  
  for (j in seq_len(ncol(m))) {
    ind <- which(masked[, j] == 0)
    add <- sample(ind, size[j])
    m[add, j] <- 1
  }
  m
}

GenerateMaskSet <- function(g, n, p){
  # Generates set of masks
  # Args:
  #   g: character vector with sequences
  #   n: number of masks to be generated    
  #   p: proportion of genotypes to be masked at each loci
  # Returns:
  #   A list of length n containing masks as matrices 
  
  # Initilize variables
  M <- nchar(g[1]) # number of markers
  N <- length(g)/2 # number of individuals
  out <- list()
  
  # Create a binary matrix with originally missing values
  m0 <- GetMissing(g, M)
  
  # Count available genotypes per marker
  size <- round((N - colSums(m0)) * p)
  
  # To keep genotypes that are already masked
  masked <- m0
  
  # Generate diffferent n masks  
  for(i in seq_len(n)) {
    
    # Initiate empty mask
    m <- matrix(0, ncol = M, nrow = N)
    
    message(sprintf("Generating mask %s...", i))
    
    # Generate new mask
    out[[i]] <- GenerateMask(m, masked, size)
    
    # Update genotypes that are alreade masked
    masked <- Reduce('+', out)
    }
  out
}

ViewMaskSet <- function(set){
  # Shows subsets (10 rows, 10 columns) of all masks in a list
  lapply(set, function(x) x[1:10, 1:10])
  
}

CountStat <- function(a0, a1){
  # Counts statistics for imputed alleles 
  # Args: 
  #  a0: original alleles 
  #  a1: imputed alleles
  # Returns:
  #  Vector with allele and genotype errors.
  
  # Check input
  if(length(a0) != length(a1)) 
    stop(sprintf("Vectors are not equal! Original: %s, imputed: %s", 
         length(a0), length(a1)))
  
  out <- a0 == a1
  
  # Count statistics
  nalleles <- length(out)
  ngenotypes <- nalleles/2
  tp <- sum(out) # Number of true positive alleles
  
  alleles <- (nalleles - tp)/nalleles
  
  # Count true positive genotypes
  t1 <- which(out == TRUE)
  t2 <- t1[-1]
  dif <- t2 - t1[-tp]
  
  # How many odd indexes do we have?
  ind <- t1[which(dif == 1)]
  
  # Count accuracy of genotype imputation
  acc <- sum(is.odd(ind))/ngenotypes
  
  # Output statistics
  c(alleles = alleles, genotypes = 1 - acc)

}

seq2mat <- function(vec){
  # Converts sequences into matrix
  
  tmp <- laply(vec, function(x) strsplit(x, split = ""))
  tmp <- ldply(tmp, function(x) x)
  
  tmp <- as.matrix(tmp)
  dimnames(tmp) <- NULL
  tmp
  
}

# Checks whether the number is odd
is.odd <- function(x) x %% 2 != 0

EstimateErrors <- function(origin, masks, imputed, K){
  # Estimates the errors of imputation.
  # Args:
  #   origin: full/path/to/filename.inp, where `filename.inp` is the original fastPHASE file  
  #   masks: full/path/to/masks.RDS
  #   imputed: vector of full/path/to/fname_genotypes.out - files with imputed genotypes.
  #            The files should be in the same order as masks were generated.
  #   K: number of clusters
  # Returns:
  #   Data frame with imputation errors  
  
  # Load original data set
  g0 <- ReadFastPHASE(origin)
  
  # Convert to matrix
  c0 <- seq2mat(g0)
  rm(g0)
  
  # Load masks
  masks <- readRDS(masks)
  
  # Initilize output
  out <- list()
  
  # Loop trough the masks
  for(i in seq_along(masks)){
    
    # Load imputed data for current mask
    g1 <- ReadFastPHASE(imputed[i])
    
    # Convert to matrix
    c1 <- seq2mat(g1)
    rm(g1)
    
    # Select mask
    m1 <- masks[[i]]
    
    # Convert into TRUE/FALSE
    m1 <- m1 == 1
    
    # Replicate rows
    m1 <- m1[rep(seq_len(nrow(m1)), each = 2),]
    
    # Subset imputed alleles by mask
    d0 <- c0[m1]
    d1 <- c1[m1]
    
    # Count statistics
    out[[i]] <- c(CountStat(d0, d1), K = K)
    
  }
  
  ldply(out, function(x) x)
}




