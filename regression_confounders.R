### explore correlation between effect length & confounders
setwd('~/github/diversity_length_correlation')
require(ape)
require(nlme)
require(corrplot)
require(MASS)

fin = 'data/species_data.csv'
ftr = 'data/bac64_protein.txt'
fplot = 'figures/distributions_confounders.png'

df = read.csv( file = fin, row.names = 1)
df <- df[ !is.na(df$seff) & df$seff > 0, ]
nb <- nrow(df)

cf <- df[ , c("leff","nstrains","ngenes","genelen","div","gc")]
# no. of variables
K <- ncol(cf)

### overall spearman correlation ###
scor.r <- round( cor(cf,method='s'), 2)
scor.p <- cor.mtest( cf, method='s')$p
scor.p_adj <- round( p.adjust( p = scor.p[1,-1], method = 'holm'), 3)

### Distributions of variables ###
png(file.path(fplot))
par.og <- par(mfrow=c(3,2), mar = c(2,2,2,2), oma = c(4,4,2,2))
sapply( 1:K, function(i) hist(x = cf[,i], main = colnames(cf)[i]))
mtext(text = paste('N = ',nrow(cf)),side = 3,outer=T)
par(par.og)
dump <- dev.off()

# remove gene-length outlier
cf <- cf[ cf$genelen < 1500, ]
nb <- nrow(cf)
# box-cox transform diversity
dv <- as.vector( cf$div)
b <- boxcox(lm(dv~1))
rl <- round(b$x[b$y==max(b$y)],1)
ldv <- log(dv)
hist(ldv, main='log div',xlab='')
cf$div <- ldv

### Collinearity analysis ###
r2s <- c(
  summary( lm( formula = nstrains ~ ngenes + genelen + div + gc, data = cf))$adj.r.squared,
  summary( lm( formula = ngenes ~ nstrains + genelen + div + gc, data = cf))$adj.r.squared,
  summary( lm( formula = genelen ~ nstrains + ngenes + div + gc, data = cf))$adj.r.squared,
  summary( lm( formula = div ~ nstrains + ngenes + genelen + gc, data = cf))$adj.r.squared,
  summary( lm( formula = gc ~ nstrains + ngenes + genelen + div, data = cf))$adj.r.squared)
names(r2s) <- colnames(cf)[-1]
vifs <- 1/(1-r2s)

### multiple linear regression ###
# scale variables with mean and SD
cf.scaled <- as.data.frame( scale(cf))
model.formula <- formula( leff ~ nstrains + ngenes + genelen + div)
olsfit <- lm( model.formula, data = cf.scaled)
ols.sm <- summary( olsfit, cor=T)

### phyloGLS ###
tree.bacs <- read.tree(ftr)
tree.bacs$tip.label <- gsub(pattern = "'",replacement = "",x = tree.bacs$tip.label,fixed = T)
cf.scaled$sp <- rownames(cf.scaled)
cf.scaled <- cf.scaled[ tree.bacs$tip.label, ]

## final gls fit
glsfit <- gls( model.formula, data = cf.scaled,
               correlation = corPagel( value = 1, phy = tree.bacs, form=~sp))
gls.sm <- summary(glsfit)

# 95% CI of pagel's lambda
plci <- round( intervals( glsfit, method='var-cov')$corStruct, 2)
## alternative model for lambda = 0
fitPagel0 <- gls( model.formula, correlation = 
                    corPagel(value = 0, phy = tree.bacs,form = ~sp, fixed = T), 
                  data = cf.scaled) 
anv <- anova(glsfit,fitPagel0)