### calculate mean and variance of nucleotide diversity over genes for each site

## arguments
args <- commandArgs(trailingOnly = T)
fin <- args[1] # input table of site-wise diversity across genes
fout <- args[2] # output table of mean diversities

# load data
div <- read.csv( file = fin, sep='\t', row.names = 1, check.names = F)
# sites
sites <- as.numeric( colnames(div))
# number of sites
nsites <- length(sites)
# number of genes with sites
ngenes <- apply( div, 2, function(x) sum(!is.na(x)))
# mean div for each site
mdivs <- colMeans( div, na.rm = T)
# variance of div for each site
vdivs <- apply( div, 2, var, na.rm=T)
out <- data.frame(site=sites,ngenes=ngenes,mdiv=mdivs,var=vdivs)
write.csv(out,fout,row.names=F,quote=F)
