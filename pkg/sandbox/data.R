###% Data generation for the NMF package
###%
###% The script generates sample data to be used with the NMF package:
###% - esGolub: an ExpressionSet object that contains the Golub dataset and all the covariates
###% - CellType: an ExpressionSet object that contains the cell-type specific blood dataset from Palmer
###%
###% @author Renaud Gaujoux
###% @creation 17.04.2009

outdir <- 'data'

#####################
# GOLUB dataset: create an ExpressionSet with ALL.AML data from Brunet
#####################
data.dir <- '../data/golub'
express <- read.delim(file.path(data.dir,'ALL_AML_data.txt'), header=FALSE)
genes <- read.delim(file.path(data.dir,'ALL_AML_genes.txt'), header=TRUE, row.names=1)
pheno <- read.delim(file.path(data.dir,'ALL_AML_samples.txt'), header=FALSE, row.names=1)	
# add columns: Name, ALL.AML, Cell
pheno$Sample	<- rownames(pheno)	
pheno$ALL.AML <- as.factor(sub("([a-z]+)_[a-z0-9]+_?.*","\\1",pheno$Sample, ignore.case=TRUE))
pheno$Cell <- as.factor(sub("[a-z]+_[a-z0-9]+_?(.*)","\\1",pheno$Sample, ignore.case=TRUE))	
rownames(express) <- rownames(genes)
colnames(express) <- pheno$Sample
phenoMetaData <- data.frame(labelDescription=c("Sample name from the file ALL_AML_data.txt","ALL/AML status","Cell type"))
genesMetaData <- data.frame(labelDescription=c("Short description of the gene"))

# create experimental data
expData <- new("MIAME", name = "Renaud Gaujoux",
    lab = "UCT NBN", contact = "renaud@cbio.uct.ac.za",
    title = "Golub dataset", 
    abstract = "NMF Package",
    url = "cbio.uct.ac.za",
    other = list(notes = "")
)

esGolub <- new('ExpressionSet', exprs=as.matrix(express)
							, phenoData=new('AnnotatedDataFrame', data=pheno, varMetadata=phenoMetaData)
							, featureData=new('AnnotatedDataFrame', genes, varMetadata=genesMetaData)
							, experimentData = expData)
save(esGolub, file='esGolub.rda')
#####################


#####################
# Cell-type
#####################
data.dir <- '../data/celltype'
# load expression data
expr <- read.table(file.path(data.dir, 'Log2 ratio data.txt'), sep='\t', header=TRUE)
rownames(expr) <- expr[, 'LUID']
feature <- expr[, c('LUID', 'NAME', 'GWEIGHT', 'Symbol')]
expr <- expr[, -which(colnames(expr) %in% colnames(feature))]

# load phenotypic data
pheno <- read.table(file.path(data.dir, 'phenodata.txt'), sep='\t', header=TRUE, row.names='ID')

# create experimental data
expData <- new("MIAME", name = "Renaud Gaujoux",
    lab = "UCT NBN", contact = "renaud@cbio.uct.ac.za",
    title = "Cell-type specific gene expression profiles of leukocytes in human peripheral blood [Palmer 2006]", 
    abstract = "Apply NMF to try to distinguish separate different cell-types",
    url = "cbio.uct.ac.za",
    other = list(notes = "The data has comes from the supplementary data of the paper.")
)

# create the expression set
CellType <- new("ExpressionSet", exprs=as.matrix(expr)
	, phenoData = new("AnnotatedDataFrame", data=pheno)
	, featureData = new("AnnotatedDataFrame", data=feature)
	, experimentData = expData)


# remove genes with missing values
NAs <- esApply(CellType, 1, function(x) any(is.na(x)))
CellType = subset(CellType, subset= !NAs)

if( FALSE ){
	desc.genes <- function(eSet)
	{
		cat('Number of probes:', nrow(eSet),"\n")
		cat('Number of genes:', length(unique(fData(eSet)$Symbol)), "\n")
		
	}

	# flag genes with expression that varies some more than a given fold over the mean across the samples
	library(genefilter)
	fold.from.mean <- function(n=3, fold=3){ function(x){ f <- abs(x/mean(x)); (sum(f  >= fold) >= n)  } }

	ff <- filterfun(fold.from.mean(3,3))
	ff <- filterfun(cv(3))
	genes.ok <- genefilter(eSet, ff) 
	sum(genes.ok)

	eSet = subset(eSet, subset= genes.ok)
}

# save data in a file
save( CellType, file=file.path(outdir, 'CellType.rda'))

############################

