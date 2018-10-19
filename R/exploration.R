# Functions easily test and compare ranks and algorithms.
# Functiosn to extract informations of the final NMF analysis.
# 
# Author: Cintia C Palu
# Creation: 25 May 2018
# Integrated to NMF: 19 Oct 2018
###############################################################################

#' @param data table to be analysed
#'
#' @param ann information on the data to be displayed on the heatmaps
#' @param r range of ranks to be tested
#' @param prefix a string used to name the generated plots
#' @param seed a number to initialise the algorithm
#' @param .opt 
#'
#' @export
#' @rdname nmfExplore
#' @aliases nmfExplore

nmfExplore <- function(data, ann, r, prefix, seed = 123456, .opt = "vP"){
  #Function to support the decision on algorithms and rank value
  
  folder <- gsub("[^[:alnum:] ]", '.', prefix)
  i <- 1
  while(file.exists(folder)){
    folder <- paste0(gsub("[^[:alnum:] ]", '.', prefix), i)
    i <- i+1
  }
  dir.create(folder)
  
  postfix <- paste(r[1], r[length(r)], sep = 'to')
  
  nrange <- nmf(data, r, method = nmfAlgorithm(), nrun = 50, .opt = .opt, seed = seed)
  save(nrange, file = paste0(folder, '/', prefix, '_nrange.rda'))
  
  for(i in names(nrange)){
    png(filename = paste0(folder, '/', prefix, '_', gsub('\\/', '.', i), '_', postfix, "_consensusmap.png"), 
        width = 1600, height = 1200, units = "px")
    par(mfrow = c(4, 3))
    consensusmap(nrange[[i]], annCol = ann)
    dev.off()
  }
  par(mfrow = c(1, 1))
  
  N.rand <- randomize(data)
  # estimate quality measures from the shuffled data
  # only run NMF using the algorithms that didn't have errors in nrange
  
  nrand <- nmf(N.rand, r, method = names(nrange), nrun = 50, .opt = .opt, seed = 123456)
  save(nrand, file = paste0(folder, '/', prefix, '_nrand.rda'))
  assign(paste0(prefix, "_nrand"), nrand, .GlobalEnv)
  
  for(i in names(nrange)){#[names(nrange)%in%names(nrand)]){
    png(filename = paste0(folder, '/', prefix, '_', gsub('\\/', '.', i), '_', postfix, "_ranksurvey.png"), 
        width = 1600, height = 1200, units = "px")
    print(plot(nrange[[i]], nrand[[i]]))#, method = i)
    dev.off()
  }
  
  return(nrange)
}

#' @param meta 
#'
#' @param fileID 
#' @param data 
#'
#' @export
#' @rdname exportAnalysis
#' @aliases exportAnalysis
 
exportAnalysis <- function(meta, fileID, data){
  g <- apply(meta, 2, function (x) which(x == max(x)))
     
  data <- cbind(g, data)
  data <- data[names(sort(g)), ]
  write.table(data, file = paste(fileID, "metagene.tsv", sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
  return(g)
}

#' @param nmf.result 
#'
#' @param method 
#' @param original.data 
#'
#' @export
#' @rdname getFeatures
#' @aliases getFeatures

getFeatures = function(nmf.result, method = 'max', original.data){
  groups <- (extractFeatures(nmf.result, method))
  cat(paste('\n\nExtract Groups -', toupper(method), '\n'))
  if( any(!is.na(groups))){
    for (g in 1:length(groups)){
      cat(paste('\nGroup', g, '\n'))
      if(is.na(groups[g])){
        cat('No top most contributing feature associated\n')
      }else{
        print(rownames(original.data)[groups[[g]]])
      }
    }
  }else{
    cat('No metagene has a top most contributing feature associated\n')
  }
  return(groups)
}

#' @param original.data the original data matrix
#'
#' @param nmf.result nmf object from which we will extract information
#' @param prefix the name of the folder to save the results
#' @param postfix a string informing the algorithm and rank used, it will also be used in the file names
#'
#' @export
#' @rdname exploreNmf
#' @aliases exploreNmf

exploreNmf <- function(original.data, nmf.result, prefix, postfix){
  # Extracting information of the main features that contribute to defining
  # each metagene, as well as identifying to which cluster the samples belong to.
  
 
  setwd(dir = paste0('./', prefix))
  rankSum <- summary(nmf.result)
  write.table(rankSum, file = paste(prefix, postfix, "rankSum.tsv", sep = '.'), sep = "\t", row.names = TRUE, 
              col.names = NA, dec = ".")
  
  ###################
  ###################
  ## Data analysis ##
  ###################
  ###################
     cat('\nStarting analysis\n')
  # Getting the matrix
  w <- basis(nmf.result)
  h <- coef(nmf.result)
  
  #fit(nmf.result)
  
  V.hat <- fitted(nmf.result)
  write.table(V.hat, file = paste(prefix, postfix, "V.HAT.tsv", sep = '.'), row.names = TRUE, col.names = NA, sep = "\t")
  
  if(is.null(rownames(original.data))){
    cat('\nNaming rows in the original dataset\n')
    rownames(original.data) = paste0('R', 1:dim(original.data)[1])
  }
  
  groups_kim <- getFeatures(nmf.result, "kim", original.data)
  groups_max <- getFeatures(nmf.result, "max", original.data)
  
  ####################
  #    Exporting     #
  # relevant Genes   #
  # for each cluster #
  ####################
  cat('\nExporting metagenes information\n')
  groups_h <- exportAnalysis(h, file = paste(prefix, postfix, ".H", sep = '.'), t(original.data))
  groups_w <- exportAnalysis(t(w), file = paste(prefix, postfix, "W", sep = '.'), original.data)
  write.table(w, file = paste(prefix, postfix, "W.tsv", sep = '.'), row.names = TRUE, col.names = NA, sep = "\t")
  write.table(h, file = paste(prefix, postfix, "H.tsv", sep = '.'), row.names = TRUE, col.names = NA, sep = "\t")
  
  cat('\nSaving variables\n')
  save(list = c('w', 'h', 'V.hat', 'groups_w', 'groups_h', 'groups_max', 'groups_kim'), 
       file = paste(prefix, postfix, "exploreNMF.rda", sep = '.'))
  setwd(dir = '../')
  cat('\nExploratory Analysis finished. ')
  cat('Environment objects saved in the file:\n"')
  cat(paste(prefix, postfix, 'exploreNMF.rda"\n', sep = '.'))
}

#' @param original.data 
#'
#' @param ann 
#' @param r 
#' @param prefix 
#' @param .opt 
#' @param alg 
#' @param maxIter 
#' @param nrun 
#'
#' @export
#' @rdname runNmf
#' @aliases runNmf

runNmf <- function (original.data, ann, r, prefix, .opt = "vP", alg, maxIter = 30000, nrun = 1000){
  
  nmf.result <- nmf(original.data, r, method = alg, .opt = .opt, maxIter = maxIter, nrun = nrun)
  
  folder <- gsub("[^[:alnum:] ]", '.', prefix)
  if(!file.exists(folder)){
    dir.create(folder)
  }
  postfix <- paste0(alg, '.r', r)
  
  save(nmf.result, file = paste0(folder, '/', folder, '.', postfix, '.rda'))
  exploreNmf(original.data = original.data, nmf.result = nmf.result, prefix = folder, postfix = postfix)
  
  return(nmf.result)
}


