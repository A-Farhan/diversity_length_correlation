### take mean over genomes of predicted min doubling times of each species

## arguments
args <- commandArgs(trailingOnly = T)
dr <- args[1] # input directory that contains species as sub-directories
fin <- args[2] # input file name of grodon results
fout <- args[3] # output table file of species and gr stats

bacs <- list.dirs( dr, recursive = F, full.names = F)
nb <- length(bacs)

mgr <- matrix( data = NA, nrow = nb, ncol = 2, dimnames = list( bacs, c('cub','gr')))
for(i in 1:nb){
  fp <- file.path( dr, bacs[i], fin)
  
  if(!file.exists(fp)){
    print( paste( "File doesn't exist for", bacs[i]))
    next
  }
  
  df <- read.csv( file = fp, row.names = 1, check.names = F)

  ms <- rowMeans( df[ c( 'CUBHE', 'd'), ], na.rm = T)
  mgr[ i, ] <- ms
}

mgr <- as.data.frame( na.omit(mgr))
write.csv( x = mgr, file = fout, quote = F, row.names = T)