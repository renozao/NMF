## This Makefile was inspired from the RcppGSL package
## Copyright (C) 2011 Romain François and Edd Eddelbuettel
## It was modifed by Renaud Gaujoux to make it more generic and to generate the 
## fake vignettes on the fly.
## Copyright (C) 2011 Renaud Gaujoux

## There is an old bug in texidvi that makes it not swallow the ~
## marker used to denote whitespace. This is actually due to fixing
## another bug whereby you could not run texidvi on directory names
## containing a tilde (as we happen to do for Debian builds of R
## alpha/beta/rc releases). The 'tilde' bug will go away as it
## reportedly has been squashed upstream but I am still bitten by it
## on Ubuntu so for now Dirk will insist on pdflatex and this helps.
AUTHOR_USER=renaud
MAKE_R_PACKAGE=NMF
MAKEPDF=1

##---------------------------------------------------------------------
## Everything below this should be generic and work for any package provided that
## they have the following directory inst/doc setting:
## - inst/vignettes/src: contains the Rnw files for normal vignettes
## - inst/vignettes/unitTests: contains an R file <PKGNAME>-unitTests.R and an .Rnw file
## that run the unit tests and generate the unit test vignette respectively
##---------------------------------------------------------------------

SRC_DIR=src
whoami=$(shell whoami)
RNW_SRCS = $(notdir $(wildcard src/*.Rnw))
PDF_OBJS=$(RNW_SRCS:.Rnw=.pdf)

# add unit tests if necessary
ifneq ("$(wildcard ../tests)","")
# TODO: reactive unit tests
#PDF_OBJS:=$(R_PACKAGE)-unitTests.pdf $(PDF_OBJS)
endif

ifneq (,$(findstring $(AUTHOR_USER),$(whoami)))
LOCAL_MODE=1
MAKEPDF=1
endif

ifeq (${R_HOME},)
NOT_CHECKING=1
TMP_INSTALL_DIR := $(shell mktemp -d)
export R_LIBS=$(shell pwd)/../../../lib
export MAKE_R_PACKAGE
endif

# Define command for temporary installation (used when make is directly called,
# i.e. when not in check/build/INSTALL)
ifdef NOT_CHECKING
define do_install
	# Installing the package in temporary library directory $(TMP_INSTALL_DIR)
	-$(RPROG) CMD INSTALL -l "$(TMP_INSTALL_DIR)" ../. >> Rinstall.log 2>> Rinstall.err
	@if test ! -d "$(TMP_INSTALL_DIR)/$(R_PACKAGE)"; then \
	echo "ERROR: Temporary installation failed: see Rinstall.log"; \
	echo "# Removing temporary library directory $(TMP_INSTALL_DIR)"; \
	rm -rf $(TMP_INSTALL_DIR); \
	exit 1; else \
	echo "# Package successfully installed"; fi;
	@echo "# Adding temporary library '$(TMP_INSTALL_DIR)' to R_LIBS"
	$(eval R_LIBS := $(TMP_INSTALL_DIR):$(R_LIBS))
endef
else
define do_install
endef	
endif

#ifdef LOCAL_MODE
#define update_inst_doc
#	# Updating PDF files in inst/doc
#	mv -f *.pdf ../inst/doc
#endef
#else
#define update_inst_doc
#endef	
#endif

all: init $(PDF_OBJS) do_clean
	@echo "# All vignettes in 'vignettes' are up to date"

init:
	# Generating vignettes for package $(R_PACKAGE)
	# Detected vignettes: $(PDF_OBJS)
ifdef LOCAL_MODE
	# Mode: Local Development
else
	# Mode: Production
endif

clean:
	rm -fr *.tex *.bbl *.blg *.aux *.out *.log *.err *-blx.bib unitTests-results

clean-all: clean
	rm -fr *.pdf ../inst/doc/*.pdf *.Rnw

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript

.SECONDARY:

do_clean:
	# Removing temporary install directory '$(TMP_INSTALL_DIR)'
	@-rm -rf $(TMP_INSTALL_DIR);

# Generate .pdf from .tex
%.pdf: ${SRC_DIR}/%.Rnw
	$(eval VIGNETTE_BASENAME := $(shell basename $@ .pdf))
	$(do_install)
	# Generating vignette $@ from ${SRC_DIR}/$*.Rnw
	# Using R_LIBS: $(R_LIBS)
	$(RSCRIPT) --vanilla -e "pkgmaker::rnw('${SRC_DIR}/$*.Rnw', '$*.tex');"
	# Generating pdf $@ from $*.tex
ifdef MAKEPDF
ifdef LOCAL_MODE
	# Using pdflatex
	@pdflatex $(VIGNETTE_BASENAME) >> $(VIGNETTE_BASENAME)-pdflatex.log
	-bibtex $(VIGNETTE_BASENAME)
	@pdflatex $(VIGNETTE_BASENAME) >> $(VIGNETTE_BASENAME)-pdflatex.log
	@pdflatex $(VIGNETTE_BASENAME) >> $(VIGNETTE_BASENAME)-pdflatex.log
	# Compact vignettes
	$(RSCRIPT) -e "tools::compactPDF('$(VIGNETTE_BASENAME).pdf')"
else
	# Using tools::texi2dvi
	$(RSCRIPT) -e "tools::texi2dvi( '$(VIGNETTE_BASENAME).tex', pdf = TRUE, clean = FALSE )"
	-bibtex $(VIGNETTE_BASENAME)
	$(RSCRIPT) -e "tools::texi2dvi( '$(VIGNETTE_BASENAME).tex', pdf = TRUE, clean = TRUE )"
endif
endif	
	# Remove temporary LaTeX files (but keep the .tex)
	rm -fr $(VIGNETTE_BASENAME).toc $(VIGNETTE_BASENAME).log \
	$(VIGNETTE_BASENAME).bbl $(VIGNETTE_BASENAME).blg \
	$(VIGNETTE_BASENAME).aux $(VIGNETTE_BASENAME).out $(VIGNETTE_BASENAME)-blx.bib	
	# Update fake vignette file ./$*.Rnw
	$(RSCRIPT) --vanilla -e "pkgmaker::makeFakeVignette('${SRC_DIR}/$*.Rnw', '$*.Rnw')"
	$(update_inst_doc)

ifneq ("$(wildcard ../tests)","")
%-unitTests.pdf: ../tests/doUnit.R
	$(do_install)
	$(eval VIGNETTE_BASENAME := $(shell basename $@ .pdf))
	# Generating vignette for unit tests: $@
	# Using R_LIBS: $(R_LIBS)
	# Cleanup unit test result directory and sweave source
	if -e $(R_PACKAGE)-unitTests.Rnw; then rm $(R_PACKAGE)-unitTests.Rnw; fi
	rm -fr unitTests-results/*
	# Run unit test file: 'unitTests/doUnit.R'
	$(RSCRIPT) --vanilla -e "source('unitTests/doUnit.R'); doUnit('$(R_PACKAGE)', '$(R_PACKAGE)-unitTests.Rnw');" >> unitTests.log
	# Sweave generated file: $(R_PACKAGE)-unitTests.Rnw
	$(RPROG) CMD Sweave $(R_PACKAGE)-unitTests.Rnw
	# Run texi2dvi on $(R_PACKAGE)-unitTests.tex
	$(RSCRIPT) --vanilla -e "tools::texi2dvi( '$(R_PACKAGE)-unitTests.tex', pdf = TRUE, clean = TRUE )"
	# Cleanup latex file $(R_PACKAGE)-unitTests.tex
	rm -fr $(R_PACKAGE)-unitTests.tex
	$(update_inst_doc)
	
endif
