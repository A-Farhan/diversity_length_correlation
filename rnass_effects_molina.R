### Compute effect size and length of RNA secondary structure selection for multiple bacterial species
## using Molina's method

args <- commandArgs(trailingOnly = T)

indr <- args[1] # input directory
fin <- args[2] # file name for tables of sites' base-pair probabilities 
fout <- args[3] # path of output table of effect size and length based on RNA structure
upr <- as.integer(args[4]) # upstream region
dnr <- as.integer(args[5]) # downstream region

bacs = list.dirs( path = indr, recursive = F)
nb <- length(bacs)

# flanking regions are 40 bases , leaving out the last 10 on both ends  
lf <- (10-upr):-51
rf <- (dnr-50+1):(dnr-10)
nfs <- length( c(lf,rf))

## smoothing function
smoother <- function( x, S=zstat, a=5){
  # smoothing window around the target site
  y <- -a:a
  ep <- exp( -abs(x - y)/a)
  N <- sum( ep, na.rm = T)
  labels <- as.character( x - y)
  stv <- S[ labels]
  num <- sum( stv * ep, na.rm = T)
  out <- num / N
  return(out)
}
###

mt <- matrix( NA, nb, 3)
rownames(mt) <- basename(bacs)
colnames(mt) <- c('ngenes', 'lss', 'maxz')

for( i in 1:nb){
  fpath <- file.path( bacs[i], fin)
  
  if(!file.exists(fpath))
    next
  
  df <- read.csv( file = fpath, check.names = F)
  # remove gene names
  df <- df[ , -1]
  # adjust column names st that first column in 100 bp upstream
  colnames(df) <- as.numeric( colnames(df)) - upr
  # remove rows with NAs
  df <- df[ !apply( df, 1, anyNA), ]
  # number of genes
  G <- nrow(df)
  # average probability for each site
  avgp <- colMeans( df)
  # variance for each site
  varp <- colSums( df * (1-df))/G
  # standard error
  sep <- sqrt(varp/G)
  # average p over flanking regions
  pflank <- ( sum( avgp[ as.character(lf)]) + sum( avgp[ as.character(rf)]))/ (nfs)
  # standard error of flanking regions
  seflank <- sqrt( sum( sep[ as.character(lf)]^2) + sum( sep[ as.character(rf)]^2) ) / (nfs+2)
  # z-statistic
  zstat <-  ( avgp - pflank) / sqrt( sep^2 + seflank^2)
  zstat <- zstat[ as.character(5:(dnr-50))]
  sites <- as.numeric( names(zstat))
  # smoothed statistic
  zsm <- sapply( sites, smoother)
  # the site where smoothed z-statistic drops to 0
  lss <- sites[ zsm < 0][1]
  # maximum value for the smoothed statistic
  maxz <- max(zsm)
  mt[ i, ] <- c( G, lss, maxz)
}

mt <- mt[ !apply( mt, 1, anyNA), ]
write.csv( x = mt, file = fout)