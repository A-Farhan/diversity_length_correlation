setwd('~/github/diversity_length_correlation')

fin = 'data/simulatedSp_summary.csv'
fout = 'figures/figure_sim.png'

# all data
df = read.csv( file = fin, row.names = 1)

## plotting device
reso = 300
imagelen = 2.5
psize = 4
png( fout, width = imagelen*2, height = imagelen*2, 
     units = "in", res = reso, pointsize = psize)
# plotting parameters
cxlab = 2.5
cxxs = 2.2
ctxt = 2.2

cols = c( rgb(red = 255,green = 165,blue = 0,alpha = 200,maxColorValue = 255))

par( oma = c(0,0.2,0.1,0), mar = c(7,7,3,2) + 0.2, mfrow = c(2,2),
     mgp = c(5,1.5,0))


# panel A : example avg-site-div
fasd = 'data/example_sim_avg_sitediv.csv'

dasd = read.csv( file = fasd)
dasd$zsite = dasd$site-1
plot( div ~ zsite, data = dasd,
      xlab = 'Site', ylab = 'Site diversity',
      col = cols[1],
      cex.lab = cxlab, cex.axis = cxxs, cex = 1.5)
# asr parameters
asrps <- df[ 'rspecies10', ]
# asr curve corresponding to these parameters
predy <- asrps[1,'dmax'] + (asrps[1,'dmin']-asrps[1,'dmax'])*exp(-dasd$zsite*asrps[1,'c'])
lines( dasd$zsite, predy, lty = 3)
abline( v = asrps[1,'leff'], lty = 3, col = 'grey60')
mtext(text = bquote( paste( italic(L)[e], ' = ', .(asrps[1,'leff']))),
      side = 3, at = asrps[1,'leff'], cex = ctxt, line = 0.5)

# panel B : correlation b/w S and Se

## scale the data with mean and SD
df.scaled <- as.data.frame( apply( df, 2, scale))

plot( seff ~ S, data = df,
      xlab = bquote( paste('Scaled selection coefficient, -',italic(S))),
      ylab = bquote( paste('Effect size, ',italic(S)[e])),
      col = cols[1],
      cex.lab = cxlab, cex.axis = cxxs, cex = 1.5)
cob <- cor.test(x = df$S, y = df$seff, method = 's')
pv <- round( cob$p.value, 7)
r <- round( cob$estimate, 2)
nb = sum(!is.na(df.scaled$seff))
ln <- bquote( paste( rho[Spearman], ' = ', .(r)))
mtext( text = ln, side = 3, adj = 1, cex = ctxt, line = 0.5)
mtext(text = bquote( paste( italic(N), ' = ', .(nb))), side = 3, adj = 0, 
      cex = ctxt, line = 0.5)

## panel C : correlation b/w sellen and Le
plot( leff ~ sellen, data = df,
      xlab = bquote( paste('Length of selection region, ',italic(L)[s])),
      ylab = bquote( paste( 'Effect length, ',italic(L)[e])),
      col = cols[1],
      cex.lab = cxlab, cex.axis = cxxs, cex = 1.5)
ob <- lm( leff ~ sellen, data = df)
abline(ob,lty=2)
sob <- summary(ob)
ln <- bquote( paste( italic(R)^2, ' = ',
                     .( round(sob$adj.r.squared, 2))))
mtext( text = ln, side = 3, adj = 1, cex = ctxt, line = 0.5)

## panel D : correlation b/w leff and divlencor
plot( scor ~ leff, data = df,
      xlab = bquote( paste( 'Effect length, ',italic(L)[e])),
      ylab = bquote( paste( 'Spearman correlation, ',rho)),
      col = cols[1],
      cex.lab = cxlab, cex.axis = cxxs, cex = 1.5)
cob <- cor.test(x = df$leff, y = df$scor, method = 's')
pv <- round( cob$p.value, 7)
r <- round( cob$estimate, 2)
ln <- bquote( paste( rho[Spearman], ' = ', .(r)))
mtext( text = ln, side = 3, adj = 0, cex = ctxt, line = 0.5)


ob <- lm( scor ~ leff, data = df)
sob <- summary(ob)
abline(ob,lty=2)
ln <- bquote( paste( italic(R)^2, ' = ',
                     .( round(sob$adj.r.squared, 2))))
mtext( text = ln, side = 3, adj = 1, cex = ctxt, line = 0.5)
sob <- summary(ob)
#######################
dump <- dev.off()