setwd('~/github/diversity_length_correlation')
require(ape)
require(phytools)
require(nlme)
require(corrplot)
require(MASS)

fin = 'data/species_data.csv'
ftr = 'data//bac64_protein.txt'

df = read.csv( file = fin, row.names = 1 )

## Distribution of rna ssp len Lss ##
hist(  df$lss, main = '', xlab = 'Lss')
iqr_lss = IQR(df$lss)
# box-cox transform Lss
vlss <- as.vector(df$lss)
b <- boxcox( lm( vlss~1), plotit = F)
rl <- round(b$x[b$y==max(b$y)],1)
lssinv <- 1/vlss
hist( lssinv, main='', xlab='1/Lss')
df$lssinv <- lssinv

## Distribution of rna ssp z-score ##
hist(  df$maxz, main = '', xlab = 'Z-score')
# box-cox transform Lss
vrnz <- as.vector(df$maxz)
b <- boxcox( lm( vrnz~1), plotit = F)
rl <- round(b$x[b$y==max(b$y)],1)
###############

df <- df[ df$seff > 0 & !is.na(df$seff), ]
df <- df[ !rownames(df) %in% c("Planctomycetes_bacterium"), ]
bacs = rownames(df)
nb <- length(bacs)

### correlation b/w sec structure effect length & sizes ###
cf <- df[, c('leff','seff','lssinv','maxz')]
nc <- ncol(cf)
pmat <- sapply( 1:nc, function(i) sapply( 1:nc, function(j) 
  cor.test( x = cf[,i], y = cf[,j], method = 's')$p.value))
colnames(pmat) <- colnames(cf)
rownames(pmat) <- colnames(pmat)
cmat <- round( cor( cf, method = 's'), 2)
rownames(cmat) <- rownames(pmat)
colnames(cmat) <- rownames(cmat)
################

### multiple linear regression ###
# add no. of strains
cf$ns <- df$nstrains
nb <- nrow(cf)
# scale variables with mean and SD
cf.scaled <- as.data.frame( scale(cf))
model.formula <- formula( leff ~ maxz + ns)
olsfit <- lm( model.formula, data = cf.scaled)
ols.sm <- summary( olsfit, cor=T)

# quality check for heteroskedasticity & normal residuals
plot( fitted(olsfit), resid(olsfit))
abline(0,0,lty=2)
hist(resid(olsfit),main='')

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
## alternative model for lambda = 0
fitPagel0 <- gls( model.formula, correlation = 
                    corPagel(value = 0, phy = tree.bacs,form = ~sp, fixed = T), 
                  data = cf.scaled) 
anv <- anova(glsfit,fitPagel0)

## quality check for heteroskedasticity & normal residuals
plot( fitted(glsfit), resid(glsfit))
abline(0,0,lty=2)
hist(resid(glsfit),main='')