# Original Matlab algorithm from Brunet et al. (2004)
# 
# Author: Renaud Gaujoux
# Created: 23 Nov 2012
###############################################################################

setNMFMethod('.M#brunet',
	, objective= 'KL'
	, algorithm = 'nmf_wrap'
	, mfiles=mfiles('brunet.m', 'nmf_wrap')
)
