## This Makefile use the Makefile provided by the package pkgmaker
## Copyright (C) 2012 Renaud Gaujoux

ifneq (${R_HOME},)
RSCRIPT=$(R_HOME)/bin/Rscript
else 
RSCRIPT=`echo R`script
endif

MK=$(shell $(RSCRIPT) --vanilla -e "pkgmaker::vignetteMakefile('NMF')")

include $(MK)
