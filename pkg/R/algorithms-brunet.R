# Original Matlab algorithm from Brunet et al. (2004)
# 
# Author: Renaud Gaujoux
# Created: 23 Nov 2012
###############################################################################

setNMFMethod('.M#brunet',
	, objective= 'KL'
	, algorithm = 'brunet_wrap'
	, mfiles=c('brunet.m', 'brunet_wrap.m')
)
