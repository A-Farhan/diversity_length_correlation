setwd('~/github/diversity_length_correlation')

scinot <- function(num,s=1){
  # scientific notation E form
  sci = formatC( num, format = "e", digits = s)
  # separate components
  parts = unlist( strsplit( x = sci, split = 'e', fixed = T))
  # decimal part
  d = parts[1]
  # exponent part
  p = as.numeric( parts[2])
  # generate x10 form and return components
  sci = paste0( d, ' x 10')
  out = c( sci, p)
  return(out)
}

## require
require(ape)
require(phytools)
require(nlme)
require(MASS)

## file paths
fin = 'data/species_data.csv'
frs = 'data/rnass_profiles.csv'
ftr = 'data/bac64_protein.txt'
fout = 'figures/leff_correlates.png'

df <- read.csv(file = fin, row.names = 1, stringsAsFactors = F)

## initialize plotting device
png( fout, width = 900, height = 600, res = 110)
par( oma = c(0,0,0,0), mar = c(4.2,4.5,2,2), mfrow = c(2,3))

# keep those with positive effect size
# across species leff distribution
ac <- df$leff[df$seff>0]
ac <- ac[!is.na(ac)]
nac <- length(ac)
acden <- density(ac)
# within ecoli leff distribution
wn <- rnorm(n = nac, mean = df["Escherichia_coli","leff"], 
            sd = df["Escherichia_coli","SEleff"])
wnden <- density(wn)

cols = c( rgb(red = 255,green = 165,blue = 0,alpha = 200,maxColorValue = 255),
          rgb(red = 135,green = 206,blue = 235,alpha = 200,maxColorValue = 255))
plot( acden, main = '', ylim = c(0,0.10), xlab = expression(italic(L)[e]), cex.lab=1.2)
polygon(acden,col=cols[1])
lines(wnden)
polygon(wnden,col=cols[2])
legend('topleft', legend = c('Within species','Across species'), fill = rev(cols), bty='n', cex=1.2)

## plot 2: rna ss profiles
rs <- read.csv( file = frs, row.names = 1, check.names = F)
nb <- nrow(rs)
sites <- as.numeric( colnames(rs))
sites <- sites[ (sites >= -50) & (sites <= 150)]
rs <- rs[ , as.character(sites)]

mx <- max(rs)
mn <- min(rs)

color = rgb(red = 255,green = 165,blue = 0,alpha = 100,maxColorValue = 255)
plot( sites, rs[1, ], type='l', ylim = c(mn,mx),col=color,
      xlab = 'Site', ylab = 'Z-score', cex.lab=1.2)
for(i in 2:nb) lines( sites, rs[ i, ], col=color)
abline( h = 0, lty = 2)
abline( v = 0, lty = 2)
mtext(text = bquote( paste( italic(N), ' = ', .(nb))), side = 3, adj = 1, cex=0.7)

## plot 3: rec rate
# exclude Planctomycetes
df <- df[ rownames(df)!="Planctomycetes_bacterium", ]
# data for positive Se
df <- df[ df$seff > 0 & !is.na(df$leff), ]
nb <- nrow(df)

# apply max-min thresholds on r/m
cf <- df[ , c("leff","rec_by_mut")]
cf <- cf[ (cf$rec_by_mut < 1000) & (cf$rec_by_mut > .001), ]
cf$rec_by_mut = log(cf$rec_by_mut)
cf$sp <- rownames(cf)
bacs <- rownames(cf)
nb <- length(bacs)

model.formula <- formula( leff ~ rec_by_mut)
olsfit <- lm( model.formula, data = cf)
ols.sm <- summary( olsfit, cor=T)

plot( leff ~ rec_by_mut, data = cf, xlab='log r/m', ylab = bquote(paste(italic(L[e]))),
      cex.lab=1.2)
abline(olsfit)
mtext(text = bquote( paste( italic(N), ' = ', .(nb))), side = 3, adj = 1, cex = 0.7)
po <- round( ols.sm$coefficients[2,4], 4)
mtext(text = bquote( paste( italic(P), ' = ', .(po))), side = 3, adj = 0, cex = 0.7)

## plot 4: rRNA count

cf <- df[ , c("leff","rrna")]
# boxcox transformation of rRNA count : sqrt
cf$sqrt = sqrt(cf$rrna)
nb <- nrow(cf)

model.formula <- formula( leff ~ sqrt)
olsfit <- lm( model.formula, data = cf)
ols.sm <- summary( olsfit, cor=T)

plot( leff ~ sqrt, data = cf, xlab = expression(sqrt('#rRNA')), ylab = bquote( paste( italic(L[e]))), cex.lab=1.2)
abline( olsfit)
po = round( ols.sm$coefficients[2,4], 5)
sp = scinot(po)
mtext(text = bquote( paste( italic(P), ' = ', .(sp[1])^.(sp[2]))), 
      side=3, adj=0, cex=0.7)
mtext(text = bquote( paste( italic(N), ' = ', .(nb))), side=3, adj=1, cex=0.7)

## plot 5: CUB
# box-cox trasnformation of CUB : sqrt
cf <- df[ , c("leff","cubhe")]
cf$sqrtcub = sqrt(cf$cubhe)
nb <- nrow(cf)

model.formula <- formula( leff ~ sqrtcub)
olsfit <- lm( model.formula, data = cf)
ols.sm <- summary( olsfit, cor=T)

plot( leff ~ sqrtcub, data = cf, xlab=expression(sqrt(CUB[HE])), ylab=bquote(paste(italic(L[e]))), cex.lab=1.2)
abline(olsfit)
mtext(text = bquote(paste( italic(N), ' = ', .(nb))), side = 3, adj = 1, cex=0.7)
p = round( ols.sm$coefficients[2,4], 7)
sp = scinot(p)
mtext( text = bquote( paste( italic(P), ' = ', .(sp[1])^.(sp[2]))), 
       side = 3, adj = 0, cex=0.7)

## plot 6: DTs

cf <- df[ , c("leff","pred_mindt_h","tnorm_expdt_h")]
# log pred dt
cf$lnpdt = log(cf$pred_mindt_h)
# log exp dt
cf$lnedt = log(cf$tnorm_expdt_h)
nb <- nrow(cf)
# regressions w/ pred dt
olspred <- lm( leff ~ lnpdt, data = cf)
olspred.sm <- summary( olspred, cor=T)
pop <- round( olspred.sm$coefficients[2,4], 4)
# regressions w/ exp dt
subcf <- cf[ !apply( cf, 1, anyNA), ]
olsexp <- lm( leff ~ lnedt, data = subcf)
olsexp.sm <- summary( olsexp, cor=T)
poe <- round( olsexp.sm$coefficients[2,4], 4)

plot( leff ~ lnpdt, data = cf, xlab = 'log min Doubling Time', ylab = bquote(italic(L[e])), pch=20, col='orange',cex.lab=1.2)
points( leff ~ lnedt, data = cf, pch=20, col='skyblue')
abline(olspred,col='orange')
abline(olsexp,col='skyblue')
legend('topright', legend = c(expression(DT[pred]),expression(DT[exp])),
       col = c('orange','skyblue'), pch=20, bty='n')
sp = scinot(pop)
se = scinot(poe)
ln <- bquote( paste( italic(P[pred]), ' = ', .(sp[1])^.(sp[2]), 
                     ', ', italic(P[exp]), ' = ', .(se[1])^.(se[2])))
mtext(text = ln, side = 4, adj = 0, line = 1, cex=0.7)
ln <- bquote( paste( italic(N[exp]), ' = ', .(nrow(subcf))))
mtext(text = ln, side = 3, adj = 1, cex=0.7)

mtext(text = c('A','B','C'), side = 3, line = -2, outer = T, 
      cex = 1.2, at = c(0.02,0.37,0.7))
mtext(text = c('D','E','F'), side = 3, line = -22, outer = T, 
      cex = 1.2, at = c(0.02,0.37,0.7))

dump <- dev.off()