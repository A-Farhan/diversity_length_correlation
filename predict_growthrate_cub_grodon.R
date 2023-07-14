## require
suppressPackageStartupMessages( library( gRodon, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages( library( Biostrings, quietly = T,  warn.conflicts = F))

## arguments
args <- commandArgs(trailingOnly = T)
# input directory of CDS FASTA files
ddr <- args[1]
# output file
fout <- file.path( dirname(ddr), 'gr_cub.csv')

# input files
fls <- list.files(ddr)
nf <- length(fls)

sprintf('# input files at path %s = %i', ddr, nf)
# function to predict growth rate
predict_gr <- function(fin){
  genes <- readDNAStringSet(fin)
  highly_expressed <- grepl( "ribosomal protein", names(genes), ignore.case = T)
  ob <- suppressWarnings( predictGrowth(genes, highly_expressed))
  ob <- sapply( ob, round, 3)
  return(ob)
}

# output 
out <- sapply( fls, function(x) predict_gr(fin = file.path( ddr, x)))
write.csv(x = out, file = fout, quote = F)