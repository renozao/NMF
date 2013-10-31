# Makefile for R package NMF
#

MAKEFLAGS += --print-directory
#MAKEFLAGS += --silent

NAME=NMF
# load version number from the DESCRIPTION file
VERSION:=$(shell while IFS=: read var1 var2; do \
	if [ "$${var1}" = "Version" ]; then \
	echo $${var2}; \
	break; \
	fi; \
	done < "DESCRIPTION";)

RBIN=R
PKG_GZ=$(NAME)_$(VERSION).tar.gz
SRC_DIR=${CURDIR}
TMP_DIR=$(SRC_DIR)/../build/$(NAME)_$(VERSION).make.${RBIN}
LIB_PATH=$(TMP_DIR)/../lib
DOC_DIR=doc
PKG_DIR=$(TMP_DIR)/$(NAME)
RWIN_DIR=/home/renaud/Rwin32
DOC_PATH=./doc/$(NAME)

all: build install

version:
	@echo
	#################################
	@echo "# MAKE for package $(NAME)"	 
	#################################
	@echo "# Version: '$(VERSION)'"
	@echo "# using code directory '$(SRC_DIR)'"
	@echo "# using build directory '$(TMP_DIR)'"
	@echo 

init: clean
	@if [ ! -d $(TMP_DIR) ]; \
	then \
	echo -n "# Create build directory... ";	mkdir -p $(TMP_DIR) && echo 'OK' ; \
	fi
	@echo -n "# Create package source directory... ";
	@mkdir $(PKG_DIR) && echo 'OK' ;

clean: version
	@if [ -d $(TMP_DIR) ]; \
	then \
	echo -n "# Clean build directory... "; \
	rm -rf $(TMP_DIR)/* \
	&& echo 'OK'; \
	fi

install: version	
	@echo "#################################"
	@echo "##       Install package       ##"
	@echo "#################################"	
	@echo "# Use install directory '$(LIB_PATH)'"
	
	@cd $(TMP_DIR); \
	if [ ! -d $(LIB_PATH) ]; \
	then \
	echo -n "# Creating install directory... "; \
	mkdir $(LIB_PATH) \
	&& echo 'OK'; \
	fi; \
	$(RBIN) CMD INSTALL -l $(LIB_PATH) $(PKG_GZ)

build: build_only

build_only: 
	@echo "#################################"
	@echo "##       Build package         ##"
	@echo "#################################"
	
	@cd $(TMP_DIR) \
	&& $(RBIN) CMD build --no-vignettes $(SRC_DIR) #\
	#&& $(RBIN) CMD build --no-vignettes $(NAME) #\
	#&& $(RBIN) CMD build --no-vignettes --binary $(NAME)

check: update
	
	@echo "#################################"
	@echo "##       Check package         ##"
	@echo "#################################"
	@echo
	
	@cd $(TMP_DIR) \
	&& $(RBIN) CMD check --no-vignettes $(NAME)

checkAll: update
	
	@echo "#################################"
	@echo "##       Check package         ##"
	@echo "#################################"
	@echo
	
	@cd $(TMP_DIR) \
	&& $(RBIN) CMD check $(NAME)

quick_check: update
	
	@echo "#################################"
	@echo "##       Check package         ##"
	@echo "#################################"
	@echo
	
	@cd $(TMP_DIR) \
	&& $(RBIN) CMD check --no-vignettes --no-examples $(NAME)

vignette: update

	@echo "#################################"
	@echo "##       Check Vignette        ##"
	@echo "#################################"
	@echo
	
	@cd $(TMP_DIR) \
	&& $(RBIN) CMD check --no-examples --no-manual --no-codoc --no-tests $(NAME)
	
check10: update
	
	@echo "#################################"
	@echo "##       Check package         ##"
	@echo "#################################"
	@echo
	
	@cd $(TMP_DIR) \
	&& R10 CMD check --no-vignettes --no-example $(NAME)

checkAll10: update
	
	@echo "#################################"
	@echo "##       Check package         ##"
	@echo "#################################"
	@echo
	
	@cd $(TMP_DIR) \
	&& R10 CMD check $(NAME)

wincheck:
	mkdir -p $(TMP_DIR)/win/library
	#&& make PKGDIR=$(PKG_DIR) RLIB=$(TMP_DIR)/win/library pkg-$(NAME) \

	cd "$(RWIN_DIR)/src/gnuwin32" \
	&& make PKGDIR=$(PKG_DIR) RLIB=$(TMP_DIR)/win/library lazyload-$(NAME) \
	&& cd "$(RWIN_DIR)/library" \
	zip -r9X $(NAME).zip $(NAME);

html: 
	# HTML build not implemented yet
	exit 1

	@echo
	#################################
	# Build HTML help
	#################################
	@echo
		
	mkdir -p $(DOC_DIR)/$(NAME)
	cp -r $(TMP_DIR)/$(NAME).Rcheck/ $(DOC_DIR)
	cd $(DOC_DIR)/$(NAME)/html && mv 00Index.html index.html 
	
update: update_no_vignette Vignettefiles
update_no_vignette: skeleton Rdfiles Datafiles extrafiles srcfiles #Rfiles
		
skeleton: init
	@cd $(PKG_DIR); \
	echo -n "# Create package skeleton... "; \
	mkdir R man; \
	cd $(TMP_DIR); \
	echo "files <- packageDescription('', '$(SRC_DIR)', fields='Collate'); files <- strsplit(files, '[\n ]'); files <- files[[1]]; file.copy(file.path('$(SRC_DIR)/R/', files), '$(PKG_DIR)/R');"  > make.doc.R ; \
	$(RBIN) CMD BATCH --no-save make.doc.R make.doc.Rout \
	&& rm make.doc.R && rm make.doc.Rout \
	&& echo 'OK'		

skeleton2:
	
	@echo -n "# Create package skeleton [PLAIN]... ";
	@rm -rf $(TMP_DIR).plain; \
	mkdir -p $(TMP_DIR).plain; \
	 
	@cd $(TMP_DIR).plain; \
	echo "source('$(SRC_DIR)/R/package.R', chdir=TRUE); make.doc('.')" > make.doc.R ; \
	$(RBIN) CMD BATCH --no-save make.doc.R > make.doc.Rout \
	&& rm make.doc.R && rm make.doc.Rout \
	&& echo 'OK'
		
Rdfiles: 
	
	@echo -n "# Update Rd files... ";	
	@rm -f $(PKG_DIR)/man/* \
	&& cp ./Rd/*.Rd $(PKG_DIR)/man \
	&& echo 'OK'
	
	@echo -n "# Update intallation files: DESCRIPTION, NAMESPACE, NEWS, etc ... "
	@cp -f ./DESCRIPTION $(PKG_DIR) \
	&& cp -f ./NAMESPACE $(PKG_DIR) \
	&& cp -f ./NEWS $(PKG_DIR) \
	&& echo 'OK'
	

Datafiles:
	@echo -n "# Copy RData files... ";		 
	@rm -rf $(PKG_DIR)/data \
	&& mkdir -p $(PKG_DIR)/data \
	&& cp ./data/*.rda $(PKG_DIR)/data \
	&& echo 'OK'
	
Vignettefiles:
	@echo -n "# Copy Vignettes files... ";		 
	@rm -rf $(PKG_DIR)/inst/doc \
	&& cp -r ./inst/doc $(PKG_DIR)/inst \
	&& echo 'OK'
	
extrafiles:
	@echo -n "# Copy extra files/directories ... ";		 
	@mkdir -p $(PKG_DIR)/inst/ \
	&& cp -rf ./inst/* $(PKG_DIR)/inst/ \
	&& rm -rf `find $(PKG_DIR)/inst/ -type d -name .svn`\
	&& echo 'OK'

srcfiles:
	@echo -n "# Copy source files/directories ... ";		 
	@mkdir -p $(PKG_DIR)/src/ \
	&& cp -rf ./src/* $(PKG_DIR)/src/ \
	&& rm -rf `find $(PKG_DIR)/src/ -type d -name .svn`\
	&& rm -rf `find $(PKG_DIR)/src/ -type f -name *.o`\
	&& rm -rf `find $(PKG_DIR)/src/ -type f -name *.so`\
	&& echo 'OK'

	
Rfiles:
	@echo -n "# Copy source R files... ";		 
	@rm -rf $(PKG_DIR)/inst/extra \
	&& mkdir -p $(PKG_DIR)/inst/extra \
	&& cp ./NMF-src.R $(PKG_DIR)/inst/extra \
	&& echo 'OK'

utest: all run_utest

run_utest:
	@echo
	#################################
	# Running Unit Tests
	#################################
	@echo
	mkdir -p ../utests \
	&& cd ../utests \
	&& mkdir -p _$(VERSION) && cd _$(VERSION) \
	&& $(RBIN) --vanilla -e "library($(NAME), lib='$(LIB_PATH)'); pkgmaker::source_files('$(SRC_DIR)/tests');"
	#&& $(RBIN) --vanilla -e "library($(NAME), lib='$(LIB_PATH)'); source('$(SRC_DIR)/dev/devtools.R', chdir=TRUE); runTestNMF('_$(VERSION)');"
	##&& $(RBIN) --vanilla -e "library($(NAME), lib='$(LIB_PATH)'); source('$(SRC_DIR)/dev/devtools.R', chdir=TRUE); runTestNMF(file='algorithms', testFun='test.port');"

#  R --vanilla -e "require(highlight); driver <- HighlightWeaveLatex(boxes = TRUE, bg = 'white' ); Sweave( '../../NMFged/inst/doc/markers.Rnw', driver = driver ); "
	