# Basic NMF seeding methods
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include registry-seed.R
NULL

## Register base seeding methods
# None: do nothing and return object unchanged
setNMFSeed('none', function(object, x, ...){object}, overwrite=TRUE)
# Random: use function rnmf
setNMFSeed('random', rnmf, overwrite=TRUE)
