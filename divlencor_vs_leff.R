### Plot estimate of spearman correlation between gene-length and diversity against Effect Length

setwd('~/github/diversity_length_correlation')
fdl = 'data/glendiv_cors.csv'
fef = 'data/species_data.csv'
fout = 'figures/divlencor_vs_leff.png'

dl = read.csv(file = fdl, row.names = 1)
ef = read.csv(file = fef, row.names = 1, stringsAsFactors = F)

dl = dl[ dl$scor > 0, , drop=F]
ef = ef[ ef$seff > 0, ]


bacs <- intersect( rownames(dl), rownames(ef))
nb <- length(bacs)
dl <- dl[ bacs, , drop=F]
ef <- ef[ bacs, ]

df <- data.frame( ef=ef$leff, cor=dl$scor)
ob <- lm( cor ~ ef, data = df)
sob <- summary(ob)
r <- round( sob$adj.r.squared, 2) 

## plotting
png(fout)
plot( df, xlab = expression(italic(L)[e]), ylab = expression(rho))
abline(ob)
mtext(text = bquote( paste( italic(N), ' = ', .(nb))), side = 3, adj = 1)
mtext(text = bquote( paste( italic(R)^2, ' = ', .(r))), side = 3, adj = 0)
dump <- dev.off()