### check correlation b/w effect length & recombination rate using PGLS
setwd('~/github/diversity_length_correlation')
require(ape)
require(phytools)
require(nlme)
require(corrplot)

fin = 'data/species_data.csv'
ftr = 'data/bac64_protein.txt'

df = read.csv( file = fin, row.names = 1)

df <- df[ !is.na(df$seff) & df$seff > 0, ]
df <- df[ !rownames(df) %in% c("Planctomycetes_bacterium"), ]
nb <- nrow(df)

## Distribution of r/m ##
# apply max-min thresholds on r/m
df <- df[ (df$rec_by_mut < 1000) & (df$rec_by_mut > .001), ]
hist( df$rec_by_mut, xlab = 'r/m', main='')
df$rec_by_mut = log(df$rec_by_mut)
hist( df$rec_by_mut, xlab = 'log r/m', main='')
bacs <- rownames(df)
nb <- length(bacs)

### multiple linear regression ###
cf <- df[ , c("leff","rec_by_mut","nstrains")]
# scale variables with mean and SD
cf.scaled <- as.data.frame( scale(cf))
model.formula <- formula( leff ~ rec_by_mut + nstrains)
olsfit <- lm( model.formula, data = cf.scaled)
ols.sm <- summary( olsfit, cor=T)

### phyloGLS ### 
tree.bacs <- read.tree(ftr)
tree.bacs$tip.label <- gsub(pattern = "'",replacement = "",x = tree.bacs$tip.label,fixed = T)
tree.bacs$node.label <- paste0( 'n', seq( nb+1, length.out=nb-1)) 
plot( tree.bacs, cex=0.5, show.node.label=T)
cf.scaled$sp <- rownames(cf.scaled)

## phylogenetic signal in the model
lambda <- seq(0,1,0.01)
logliks <- sapply( lambda, function(x) logLik( gls(
  model.formula, data = cf.scaled,
  correlation = corPagel(value = x, phy = tree.bacs, form=~sp, fixed = T),
  method = 'ML'
)))

plot( lambda, logliks, type='l')
# final gls fit
glsfit <- gls( model.formula, data = cf.scaled,
               correlation = corPagel( value=1, phy=tree.bacs, form=~sp))
gls.sm <- summary( glsfit)

# 95% CI of pagel's lambda
plci <- round( intervals( glsfit, method='var-cov')$corStruct, 2)