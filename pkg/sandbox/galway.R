

# install needed packages
package.dir <- 'C:/Documents and Settings/User/My Documents/My Programs/R/'
dep <- c('roxygen_0.1.zip', 'RColorBrewer_1.0-2.zip', 'R.utils_1.1.7.zip', 'R.oo_1.4.8.zip', 'R.methodsS3_1.0.3.zip')
install.packages(paste(package.dir, dep, sep='/'), repos=NULL)
dep.net <- c('roxygen', 'RColorBrewer', 'R.utils')
install.packages(dep.net, dep=TRUE)


disk.letter <- 'G'
setwd(paste(disk.letter, ':/Documents/UCT/projects/NMF/Rnmf/tmp', sep=''))
source('../package/R/package.R', chdir=TRUE)

# matrix
n <- 50; r <- 3; m <- 20
V <- syntheticNMF(n, r, m, noise=TRUE)
res <- nmf(V,r)

# ExpressionSet
load('../package/data/esGolub.rda')
res <- nmf(esGolub,r)

# test suite
source('../package/R/devtools.R', chdir=TRUE)
runTestNMF('test.src')

nmf.options(track.interval=40)
res <- nmf(V, r, list('bru', 'lee', 'ns', 'off', 'snmf/r'), track=TRUE)


disk.letter <- 'D'
setwd(paste(disk.letter, ':/Documents/UCT/projects/NMF/Rnmf/tmp', sep=''))
library('NMF', lib='../lib')
source('../package/R/devtools.R', chdir=TRUE)
runTestNMF('test.pkg')

data(esGolub)
res <- nmf(esGolub, 3, list('br', 'R_br', 'le', 'R_le', 'off', 'R_of', 'ln', 'R_ln'), seed=123456)
res <- nmf(esGolub, 3, 'lee', seed=123456)
reso <- nmf(esGolub, 3, 'R_lee', seed=123456)
