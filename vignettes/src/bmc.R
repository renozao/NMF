# Scripts runs to produce the figures in the BMC paper
# 
# Author: renaud
###############################################################################

# install and load NMF package
lib.dir <- 'lib'
dir.create(lib.dir, showWarnings=FALSE)
install.packages('NMF_0.1.tar.gz', repos=NULL, lib=lib.dir)
library(NMF, lib=lib.dir)

# define a seed
.seed <- 123456

# load Golub data
data(esGolub)
#esGolub <- syntheticNMF(500, 3, 20, noise=TRUE)

# estimate rank for Golub dataset
rank.nrun <- 50
rank.range <- seq(2,6)
res.estimate <- nmfEstimateRank(esGolub, rank.range, method='brunet'
				, nrun=rank.nrun, conf.interval=TRUE, seed=.seed)
save(res.estimate, file='res.estimate.rda')

# Full run of Brunet algorithm
nmf.nrun <- 200
res.brunet <- nmf(esGolub, 3, 'brunet', nrun=nmf.nrun, seed=.seed, .options='tv')
save(res.brunet, file='res.brunet.rda')

# Comparison of methods
res.comp <- nmf(esGolub, 3, list('brunet', 'lee', 'ns', 'lnmf'), seed='nndsvd', .options='tv')
save(res.comp, file='res.comp.rda')

# save all in one file
save(res.estimate, res.brunet, res.comp, file='res.bmc.rda')

if( FALSE ){
# generate plots
png('consensus.png')
metaHeatmap(res.brunet, class=esGolub$Cell)
dev.off()

png('metagenes.png')
metaHeatmap(fit(res.brunet), class=esGolub$Cell)
dev.off()
}