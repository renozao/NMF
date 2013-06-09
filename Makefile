## This Makefile use the package Makefile provided by the package pkgmaker
## Copyright (C) 2013 Renaud Gaujoux

RSCRIPT=Rscript

MK=$(shell $(RSCRIPT) --vanilla -e "pkgmaker:::packageMakefile()")
include $(MK)

