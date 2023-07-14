### check correlation b/w effect length & growth rates using PGLS
setwd('~/github/diversity_length_correlation')
require(ape)
require(phytools)
require(nlme)
require(corrplot)
require(MASS)

fin = 'data/species_data.csv'
ftr = 'data/bac64_protein.txt'

df = read.csv( file = fin, row.names = 1)

df <- df[ !is.na(df$seff) & df$seff > 0, ]
df <- df[ !rownames(df) %in% c("Planctomycetes_bacterium"), ]
nb <- nrow(df)

## Distribution of pred minDT ##
X <- as.vector(df$pred_mindt_h)
hist(  X, main = '', xlab = 'pred minDT')
# box-cox transform
b <- boxcox( lm(X~1), plotit = F)
rl <- round(b$x[b$y==max(b$y)],1)
lx <- log(X)
hist( lx, main = '', xlab = 'log pred minDT')
df$lnpdt = lx

## Distribution of CUB-HE ##
X <- as.vector(df$cubhe)
hist(  X, main = '', xlab = 'CUB-HE')
# box-cox transform 
b <- boxcox( lm( X~1), plotit = F)
rl <- round(b$x[b$y==max(b$y)],1)

## Distribution of rRNA count
X <- as.vector(df$rrna)
hist(  X, main = '', xlab = '#rRNA')
# box-cox transform 
b <- boxcox( lm( X~1), plotit = F)
rl <- round(b$x[b$y==max(b$y)],1)

###############

cf <- df[ , c("leff","rrna","cubhe","pred_mindt_h","nstrains")]
nb <- nrow(cf)

### correlation among variables ###
nc <- ncol(cf)
pmat <- sapply( 1:nc, function(i) sapply( 1:nc, function(j) 
  cor.test( x = cf[,i], y = cf[,j], method = 's')$p.value))
colnames(pmat) <- colnames(cf)
rownames(pmat) <- colnames(pmat)
cmat <- round( cor( cf, method = 's'), 2)
rownames(cmat) <- rownames(pmat)
colnames(cmat) <- rownames(cmat)
################

### Collinearity analysis ###
r2s <- c(
  summary( lm( formula = rrna ~ cubhe + pred_mindt_h + nstrains, data = cf))$adj.r.squared,
  summary( lm( formula = cubhe ~ rrna + pred_mindt_h + nstrains, data = cf))$adj.r.squared,
  summary( lm( formula = pred_mindt_h ~ rrna + cubhe + nstrains, data = cf))$adj.r.squared,
  summary( lm( formula = nstrains ~ rrna + cubhe + pred_mindt_h, data = cf))$adj.r.squared)
names(r2s) <- colnames(cf)[c(1,2,3,5)]
vifs <- 1/(1-r2s)

### multiple linear regression ###
# scale variables with mean and SD
cf.scaled <- as.data.frame( scale(cf))
model.formula <- formula( leff ~ rrna + cubhe + nstrains)
olsfit <- lm( model.formula, data = cf.scaled)
ols.sm <- summary( olsfit, cor=T)

# quality check for heteroskedasticity & normal residuals
par.og <- par( mfrow = c(1,2), oma = c(0,0,3,0))
plot( fitted(olsfit), resid(olsfit))
abline(0,0,lty=2)
hist(resid(olsfit),main='')
mtext(text = 'Quality checks for OLS',side = 3,outer = T)
mtext(text = paste('N = ',nb),side = 3,adj = 1,outer = T)
par(par.og)

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
mtext( text = bquote( paste( italic(N), ' = ', .(nb))), side = 3, 
       adj = 0.9, line = 2, cex = 1.2)
# final gls fit
glsfit <- gls( model.formula, data = cf.scaled,
               correlation = corPagel( value=1, phy=tree.bacs, form=~sp))
gls.sm <- summary( glsfit)

# 95% CI of pagel's lambda
plci <- round( intervals( glsfit, method='var-cov')$corStruct, 2)
ln <- bquote( paste( lambda,'=',.(plci[2]),'(',.(plci[1]),', ',.(plci[3]),')'))
mtext( text=ln, side=3, line=0, adj=1)

## quality check for heteroskedasticity & normal residuals
par.og <- par( mfrow = c(1,2), oma = c(0,0,3,0))
plot( fitted(glsfit), resid(glsfit))
abline(0,0,lty=2)
hist(resid(glsfit),main='')
mtext(text = 'Quality checks for GLS',side = 3,outer = T)
par(par.og)
