library(WGCNA)
library(CoExpNets)




#' Modules enrichment
#'
#' This function determine the module enrichment of a GCN
#' @param net A GCN created with \code{\link{createGCN}}. It can be an R object or the full path of a file
#' containing it. If a filename is provided this function will first look for a file with the same name but
#' ending in ".enrichment.csv". If the file exists the function will return its content.
#' @param markers.path Folder containing user-defined lists of genes to be used as marker genes to determine
#' modules enrichment. They must have in the first row the name of the enrichment and then, a gene ID on every
#' row. Gene IDs must be expresed using the same format as in the network specified in the first parameter.
#' @param return.processed When True the functions returns the -log10 of the p-value obtained replacing the -inf
#' values by the highest value obtained. When False the function returns the p-values. This parameter is False by default.
#' @param significanceThreshold When this value is lower than 1 this function returns the name of the celltype
#' enrichment of each module which is exclusively enriched by a single celltype with a p-value lower than the
#' value received by this parameter. In this case, the value of return.processed is ignored.
#' @param outputCorrectedPvalues When true some corrections are made to the standard fisher test obtained. Firstly a Bonferroni
#' correction, then a regression correction to reduce bias of using too big enrichment sets. This paremeter is True by default.
#' @param debug Show some messages about the test made
#' @return a matrix with a column for each module in the network and a row for each enrichment type. Values in the matrix
#' reflects the p-value resulting from the test made to determine if the module is enriched or not. Alternatively, when
#' significanceThreshold < 1 a named vector with the enrichment name for each module significantly and exclusivelly
#' enriched by a celltype.
#' @export
getModulesEnrichment <- function(net,
                                 markers.path = ".",
                                 return.processed=F,
                                 significanceThreshold = 1,
                                 outputCorrectedPvalues = T,
                                 debug = F) {
  
  my.fisher.test = function( a, b, shared, total ) {
    m = matrix(c(shared, b-shared, a-shared, total-a-b+shared), ncol=2,nrow=2)
    f = fisher.test(m, alternative="greater")
    return(f$p.value)
  }
  
  if ( typeof(net) == "character" ) {
    enrichment.filename = paste0( net, ".enrichment.csv")
    if ( file.exists( enrichment.filename ) ) {
      enrichment.by.module = read.csv( enrichment.filename, stringsAsFactors=F)
      if ( colnames(enrichment.by.module)[1] == "X" ) {
        rownames(enrichment.by.module) = enrichment.by.module[,1]
        enrichment.by.module = enrichment.by.module[,-1]
      }
      return (enrichment.by.module)
    }
    net = readRDS(net)
  }
  
  modules = unique(net$moduleColors)
  markers.set = getMarkerGenes( markers.path )
  print(markers.set)
  # We remove duplicated and genes not included in the network we are going to analyze to reduce bias
  for ( i in 1:length(markers.set) ) {
    if ( debug ) {
      print( paste0( "Enrichment file", names(markers.set)[i], " contains ", 
                     length(markers.set[[i]]), " genes"))
    }
    markers.set[[i]] = unique(markers.set[[i]])
    markers.set[[i]] = markers.set[[i]][is.element(markers.set[[i]], names(net$moduleColors))]
    if ( debug ) {
      print( paste0( "After removing duplicated and genes not included in the network we have: ", 
                     length(markers.set[[i]]), " genes"))
    }
  }
  m = NULL
  for ( module in modules ) {
    genes.module = names(net$moduleColors[net$moduleColors == module] )
    n = NULL
    for ( enrichment.name in names(markers.set) ) {
      markers = markers.set[[enrichment.name]]
      shared = intersect( genes.module, markers )
      size.module = table(net$moduleColors)[module]
      size.markers = length(markers)
      size.shared = length(shared)
      ft = my.fisher.test( size.module, size.markers, size.shared, length(net$moduleColors) )
      if ( debug ) {
        corrected.p.value = length(modules) * length(markers.set) * ft
        corrected.p.value = min(1, corrected.p.value)
        print( paste0( "Module ", module, " has ", size.module, " genes, ", size.shared, 
                       " of which are also included in the ", size.markers, 
                       " genes from enrichment ", enrichment.name, " p-value: ", ft, 
                       " corrected p-value: ", corrected.p.value ) )
      }
      n = c(n, ft)
    }
    m = cbind(m,n)
  }
  
  # Change 0 by minimum possible value in R
  m = apply(m, c(1,2), FUN=function(x) {max(.Machine$double.xmin, x)})
  
  if ( outputCorrectedPvalues ) {
    # Bonferroni correction
    m = matrix(p.adjust( m,  method = "bonferroni", n = length(modules) * length(markers.set) ), 
               nrow = length(markers.set))
    #m = m * length( modules ) * length( markers.set )
    #m = apply(m, c(1,2), FUN=function(x) {min(1,x)})
    
    # Regresion correction
    x = unlist(lapply( markers.set, length))
    names(x) = names(markers.set)
    y = apply( m, 1, 
               FUN = function(x) {sum(-log10(x))/length(x) } )
    mylm = lm( y ~ x )
    base = predict( mylm, data.frame(x=unlist(lapply( markers.set, length))))
    m = ( 10^-(-log10(m) - base) )
    
    # Filter out too big values
    m = ifelse( m >= 1, 1, m)
  }
  print(m)
  if(sum(m < significanceThreshold) == 0){
    warning("This is weird: getModulesEnrichment detects no enrichment. Perhaps genes are not symbols?")
  }
  
  rownames(m) = names(markers.set)
  colnames(m) = modules
  
  # Filter out non significantly and exclusively enriched modules
  if ( significanceThreshold < 1 ) {
    if ( outputCorrectedPvalues ) {
      # significanceThreshold = significanceThreshold / ncol(m) * length(markers.set)
    }
    m = m[, as.vector(apply( m, 2, FUN = function(x) { sum( x < significanceThreshold ) == 1  } ) ), drop=F]
    return ( apply( m, 2, FUN = function(x) { rownames(m)[which(x <= significanceThreshold)] } )	)
  }
  
  # Process values transforming them in the range 0..Max
  if ( return.processed ) {
    m = -log10(m)
    if ( sum( is.infinite( m ) ) > 0 ) {
      m[ is.infinite( m ) ] = max( m[ !is.infinite( m ) ] )
    }
  }
  
  return( m )
}


#' Cleaning the primary enrichment signal
#'
#' This function generates a new expression profile by removing the primary signal of the groups of genes whose only and primary
#' enrichment is indicated in the following parameter.
#' @param expr.data The expression data to be cleaned. It can be a full path file name or a data frame with
#' genes in columns and samples in rows.
#' @param target.enrichment A string with the kind of enrichment that should be removed. This name (or list of
#' names) should be one from the list of enrichment marker genes' files provided.
#' @param net A GCN created with \code{\link{createGCN}} from the same expression data. It can be a full path
#' of a file containing it or an R object. This network reflects the primary signal and it is used to determine
#' what modules need to be cleaned in order to find the secondary signal. If no network is provided it will
#' be created. Therefore, providing one is also a way to reduce the time spent by this function.
#' @param significanceThreshold Max enrichment p-value for a module to be considred as significantly enriched.
#' @param modules.to.clean List of modules to be filtered out. This parameter is NULL by default, when a list of
#' module names is provided target and significanceThreshold parameters are ignored.
#' @param markers.path Folder containing user-defined lists of genes to be used as marker genes to determine
#' modules' enrichment. This is done using  \code{\link{WGCNA::userListEnrichment}} function, so they must be in a compatible
#' format. Gene IDs must be expresed using the same format as in the expression data specified in the first parameter.
#' @return  the expression data filtered that, hopefully, show the secondary role of some genes in the same format
#' of the input expression data. The function returns NULL when no data has been filtered out.
#' @export
removePrimaryEffect <- function( expr.data, target.enrichment, net = NULL,
                                 significanceThreshold = 0.05,
                                 modules.to.clean = NULL, markers.path = NULL,
                                 how = "pca",
                                 nComp = 20,
                                 corThreshold = 0.8) {
  
  if ( typeof(expr.data) == "character" ) {
    expr.data = readRDS(expr.data)
  }
  net.filename = "."
  if ( is.null( net ) ) {
    net = getDownstreamNetwork(tissue="GPCNA",
                               expr.data=expr.data,
                               job.path = NULL) 
    #net = createGCN( expr.data )
  } else {
    if ( typeof( net ) == "character" ) {
      net.filename = net
      net = readRDS( net.filename )
    }
  }
  #JB
  modules.to.clean = modules.to.clean[modules.to.clean != "grey"]
  #JB
  if ( ( is.null(target.enrichment) || length(target.enrichment) == 0 ) && ( is.null(modules.to.clean) || length(modules.to.clean) == 0 ) ) {
    stop( "target.enrichment must be a character list containing one or mor enrichment names when modules.to clean is NULL or empty" )
    return (NULL)
  }
  if ( is.null( modules.to.clean ) ) {
    enriched.modules = getModulesEnrichment( net = net, markers.path = markers.path,
                                             significanceThreshold = significanceThreshold )
    
    modules.to.clean = names(enriched.modules[ enriched.modules %in% target.enrichment])
    print("Remove primary effect")
    print(modules.to.clean)
  }
  if ( length(modules.to.clean) == 0) return (NULL)
  cleaned.data = expr.data[,!(net$moduleColors %in% modules.to.clean)]
  
  rmvGenes = 0
  for ( module in modules.to.clean ) {
    ee = net$MEs[paste0("ME", module)]
    expr = expr.data[,net$moduleColors == module]
    cors = cor(ee,expr)
    keepGenes = abs(cors) < corThreshold
    expr = expr[,keepGenes]
    lrmvGenes = sum(net$moduleColors == module) - sum(keepGenes)
    rmvGenes = rmvGenes + lrmvGenes
    cat("Removing", lrmvGenes,"genes\n" )
    cat("Cleaning module ", module,"\n")
    if(how == "reg"){
      res = apply( expr, 2, function(y) { lm( y ~ ., data=ee )$residuals })
    }else{#We apply PCA based method
      mu = colMeans(expr)
      axes = prcomp(expr,scale=T,rank=nComp)
      expr = axes$x[,2:nComp] %*% t(axes$rotation[,2:nComp])
      res = scale(expr, center = -mu, scale = FALSE)
    }
    cleaned.data = cbind(cleaned.data, res)
  }
  cat("Total of genes removed",rmvGenes,"\n")
  return ( cleaned.data )
}

#' Resume changes in enrichment for each gene
#'
#' This function analyces two networks determining the enrichment of each gene in each network. Both
#' networks must be made from the same set of genes. The list of genes studied can be provided or all of them
#' will be analyzed.
#' @param primary.net A GCN created with \code{\link{createGCN}}
#' @param secondary.net A GCN created with \code{\link{createGCN}}
#' @param significanceThreshold Max enrichment p-value for a module to be considred as significantly enriched.
#' @param genes The list of gene's IDs to be used
#' @param markers.path Folder containing user-defined lists of genes to be used as marker genes to determine
#' modules' enrichment. This is done using WGCNA::userListEnrichment function, so they must be in a compatible
#' format. Gene IDs must be expresed using the same format as in the networks.
#' @return A data frame showing a row for each gene and the module and enrichment of that gene in both networks.
#' When the module is not enriched a - is shonw. When is enriched by more than a type, they are shown as a list.
#' @export
enrichmentEvolution <- function( primary.net, secondary.net, significanceThreshold = 0.05, 
                                 genes = NULL, markers.path = "." ) {
  remove.na = function( x ) { x[!is.na(x)] }
  
  primary.net.filename = "."
  if ( typeof(primary.net) == "character" ) {
    primary.net.filename = primary.net
    primary.net = readRDS(primary.net)
  }
  secondary.net.filename = "."
  if ( typeof(secondary.net) == "character" ) {
    secondary.net.filename = secondary.net
    secondary.net = readRDS(secondary.net)
  }
  
  # primary.em = getModulesEnrichment( net = primary.net, markers.path = markers.path, significanceThreshold = significanceThreshold )
  # secondary.em = getModulesEnrichment( net = secondary.net, markers.path = markers.path, significanceThreshold = significanceThreshold )
  #
  # if ( is.null(genes) ) {
  #   genes = names(primary.net$moduleColors)
  # }
  # tabla = data.frame( gene = genes,
  #                     primary.module = as.character(primary.net$moduleColors[genes]),
  #                     primary.enrichment = as.character(primary.em[primary.net$moduleColors[genes]]),
  #                     secondary.module = as.character(secondary.net$moduleColors[genes]),
  #                     secondary.enrichment = as.character(secondary.em[secondary.net$moduleColors[genes]]),
  #                     stringsAsFactors=F,
  #                     row.names=NULL )
  # tabla$primary.enrichment[is.na(tabla$primary.enrichment)] = unlist(apply(primary.enrichment.by.module[,tabla$primary.module[is.na(tabla$primary.enrichment)],drop=F ], 2, FUN=function(x) { if (sum(x<=significanceThreshold) == 0) return("-") else return (paste(primary.enrichment.names[which(x <= significanceThreshold)], collapse=", ")) }))
  # tabla$secondary.enrichment[is.na(tabla$secondary.enrichment)] = unlist(apply(secondary.enrichment.by.module[,remove.na(tabla$secondary.module[is.na(tabla$secondary.enrichment)]),drop=F ], 2, FUN=function(x) { if (sum(x<=significanceThreshold) == 0) return("-") else return (paste(secondary.enrichment.names[which(x <= significanceThreshold)], collapse=", ")) }))
  # tabla$secondary.enrichment[is.na(tabla$secondary.enrichment)] = "-"
  # return(tabla)
  
  
  primary.enrichment.by.module = getModulesEnrichment( net = primary.net, markers.path = markers.path)
  secondary.enrichment.by.module = getModulesEnrichment( net = secondary.net, markers.path = markers.path)
  
  # Nos quedamos por un lado con los nombres de las columnas ( tipos de celula )
  primary.enrichment.names = rownames(primary.enrichment.by.module)
  # Nos quedamos sólo con las columnas  que tienen modulos significativa y exclusivamente enriquecidos por un único tipo de célula
  primary.em = primary.enrichment.by.module[,as.vector(apply(primary.enrichment.by.module, 2, 
                                                             FUN = function(x) { 
                                                               sum( x <= significanceThreshold)} ) == 1), 
                                            drop = F]
  # la lista de módulos etiquetada por el tipo de célula en el que están significativa y exclusivamente enriquecidos
  primary.em = apply( primary.em, 2, FUN = function(x) { primary.enrichment.names[which(x <= significanceThreshold)] } )
  
  # Nos quedamos por un lado con los nombres de las columnas ( tipos de celula )
  secondary.enrichment.names = rownames(secondary.enrichment.by.module)
  # Nos quedamos sólo con las columnas  que tienen modulos significativa y exclusivamente enriquecidos por un único tipo de célula
  secondary.em = secondary.enrichment.by.module[,as.vector(apply(secondary.enrichment.by.module, 2, FUN = function(x) { sum( x <= significanceThreshold)} ) == 1), drop=F]
  # la lista de módulos etiquetada por el tipo de célula en el que están significativa y exclusivamente enriquecidos
  secondary.em = apply( secondary.em, 2, FUN = function(x) { secondary.enrichment.names[which(x <= significanceThreshold)] } )
  
  if ( is.null(genes) ) {
    genes = names(primary.net$moduleColors)
  }
  tabla = data.frame( gene = genes,
                      primary.module = as.character(primary.net$moduleColors[genes]),
                      primary.enrichment = as.character(primary.em[primary.net$moduleColors[genes]]),
                      secondary.module = as.character(secondary.net$moduleColors[genes]),
                      secondary.enrichment = as.character(secondary.em[secondary.net$moduleColors[genes]]),
                      stringsAsFactors=F,
                      row.names=NULL )
  tabla$primary.enrichment[is.na(tabla$primary.enrichment)] = unlist(apply(primary.enrichment.by.module[,tabla$primary.module[is.na(tabla$primary.enrichment)],drop=F ], 2, FUN=function(x) { if (sum(x<=significanceThreshold) == 0) return("-") else return (paste(primary.enrichment.names[which(x <= significanceThreshold)], collapse=", ")) }))
  tabla$secondary.enrichment[is.na(tabla$secondary.enrichment)] = unlist(apply(secondary.enrichment.by.module[,remove.na(tabla$secondary.module[is.na(tabla$secondary.enrichment)]),drop=F ], 2, FUN=function(x) { if (sum(x<=significanceThreshold) == 0) return("-") else return (paste(secondary.enrichment.names[which(x <= significanceThreshold)], collapse=", ")) }))
  tabla$secondary.enrichment[is.na(tabla$secondary.enrichment)] = "-"
  return(tabla)
  
}

#' Resume changes in the variance of genes after removing some enrichment signal
#'
#' This function compare the variance of expression data before and after removing some enrichment signal using
#' \code{\link{removePrimaryEffect}}.
#' @param original.expr.data The original expression data. It can be a full path file name or a data frame with
#' genes in columns and samples in rows.
#' @param cleaned.expr.data The expression data obtained by removing some enrichment signal with \code{\link{removePrimaryEffect}}.
#' @param genes The list of gene's IDs to be used
#' @return  A data frame showing a row for each gene and the variance of that gene in both data sets.
#' @export
varianceEvolution <- function( original.expr.data, cleaned.expr.data, genes = NULL ) {
  if ( typeof(original.expr.data) == "character" ) {
    original.expr.data = readRDS( original.expr.data )
  }
  if ( typeof(cleaned.expr.data) == "character" ) {
    cleaned.expr.data = readRDS( cleaned.expr.data )
  }
  
  original.var = apply(original.expr.data, 2, var)
  cleaned.var = apply(cleaned.expr.data, 2, var)
  
  if ( is.null(genes) ) {
    genes = colnames( original.expr.data)
  }
  return (data.frame( gene = genes,
                      original.variance = original.var[genes],
                      cleaned.variance = cleaned.var[genes],
                      stringsAsFactors = F,
                      row.names=NULL ) )
}

#' Resume Read marker gene files and return a list of marker genes
#'
#' This function read marker gene files and return a list of marker genes.
#' @param markers.path Folder containing user-defined lists of genes to be used as marker genes to determine
#' modules' enrichment. This is done using WGCNA::userListEnrichment function, so they must be in a compatible
#' format. Gene IDs must be expresed using the same format as in the expression data specified in the first parameter.
#' @return  A list of marker gene lists, one for each file read.
#' @export
getMarkerGenes <- function( markers.path = system.file("markers", "",
                                                       package = "GPCNA") ) {
  files = list.files(path=markers.path, full.names = T)
  files = files[grep(pattern=".txt$", files)]
  markernames = gsub(".txt","", basename(files))
  if ( length(markernames) == 0 ) {
    stop( paste0( "There are no txt files to be used as marker genes in path ", markers.path) )
  }
  markers = list()
  for ( marker in markernames ) {
    #markers[[marker]] = CoExpNets::fromEnsembl2GeneName(read.csv( paste0( markers.path, "/", marker, ".txt" ), 
    #                              header=T, stringsAsFactors = F)[,1])
    markers[[marker]] = read.csv( paste0( markers.path, "/", marker, ".txt" ), 
                                  header=T, stringsAsFactors = F)[,1]
  }
  return(markers)
}

getCellSubtypes = function(markers.path=system.file("markers", "",
                                                    package = "GPCNA")){
  markers.set = getMarkerGenes( markers.path )
  list( MA = names(markers.set)[grep(toupper("^Astrocyte_"), 
                                     toupper(names(markers.set)))],
        OLG = names(markers.set)[grep(toupper("^oligodendrocyte_"), 
                                      toupper(names(markers.set)))],
        N = names(markers.set)[grep(toupper("^Neuron_"), 
                                    toupper(names(markers.set)))],
        MG = names(markers.set)[grep(toupper("^microglia_"), 
                                     toupper(names(markers.set)))])
}

#' Title Process an enrichment matrix
#' 
#' This function takes an enrichment matrix and generates a vectors of enrichments per module.
#' Determines the module enrichment of a GCN taking into account subsets of markers. 
#' Each subset is called meta-enrichment.
#' 
#'
#' @param enrichment.by.module The enrichment matrix. Usually, is the matrix 
#' generated by \code{\link{getModulesEnrichment}}.
#' @param celltypes The cell types we want to test
#' @param cellsubtypes Named list of user-defined lists of genes. Each entry of this list is a 
#' vector of a subset of the filenames used by \code{\link{getModulesEnrichment}}. 
#' All the genes included in the files will be considered part of the meta-enrichment whose
#' name is the name of the list.
#' @param significanceThreshold When this value is lower than 1 this function returns the name of the celltype
#' enrichment of each module which is exclusively enriched by a single celltype with a p-value lower than the
#' value received by this parameter.
#' @param singleSignal If set to TRUE only modules with enrichment for a single cell type are considered 
#' @param as.is If set to TRUE, a list with ALL signals enrichment and p-values are returned
#' @return a named vector of meta-enrichment names. Each entry is named after the 
#' module it has been found. Only those modules including
#' only one or more significantly p-values from the same meta-enrichment 
#' will receive the meta-enrichment name. The rest will have the '-'
#' character.
#' @export
#'
#' @examples
aggModEnrichmentConf = function( enrichment.by.module,
                                 celltypes = c("N","MG","OLG","MA"),
                                 cellsubtypes = getCellSubtypes(), 
                                 significanceThreshold = 0.05,
                                 singleSignal = T,
                                 as.is = T) {
  
  enrichment.names = rownames( enrichment.by.module )
  #print(enrichment.names)
  toret = apply( enrichment.by.module, 2,
                 #Work this function for each module´s signals
                 FUN = function(x) {
                   signalmask = x <= significanceThreshold
                   if(sum(signalmask)){
                     detectedcts = enrichment.names[signalmask]
                     #First get the key
                     #print(detectedcts)
                     keys = unique(unlist(lapply(detectedcts,function(x){
                       strsplit(x,"_")[[1]][[1]]
                     })))
                     #print(keys)
                     ctids = unlist(lapply(keys,function(x){
                       names(cellsubtypes)[unlist(lapply(cellsubtypes,function(y){
                         return(length(grep(x,y)) > 0)
                       }))]
                     }))
                     #print(ctids)
                     if(singleSignal & length(ctids) == 1){
                       return(list(pval=min(x[signalmask]),enr=ctids))
                     }else if(singleSignal){
                       return(list(pval=1,enr="-"))
                     }
                     if(!singleSignal){
                       pvals = NULL
                       for(key in keys){
                         #Separate the marker sets for this key
                         mks4singlekey = detectedcts[grep(key,detectedcts)]
                         pvals = c(pvals,min(x[names(x) %in% mks4singlekey]))
                       }
                       return(list(pval=pvals,enr=ctids))
                     }
                   }
                   return(NULL)
                 })
  if(as.is)
    return(toret)
  mods = NULL
  signals = NULL
  for(i in 1:length(toret)){
    mod = names(toret)[i]
    signal = "-"
    if(!is.null(toret[[i]]$enr)){
      signal = toret[[i]]$enr[which.min(toret[[i]]$pval)]
    }
    mods = c(mods,mod)
    signals = c(signals,signal)
  }
  names(signals) = mods
  return(signals)
}

#' Title Generate a function catalog from a gene expression profiling
#'
#' This is the main method of GMSCA It needs a gene expression matrix, and markers to test for enrichment
#' of the modules. Gene names are expected to be at column names in the matrix. Markers and gene names 
#' must be in the same naming space.
#' 
#' 
#' 
#' @param expr.data Gene expression profile, genes at columns, samples at rows. Column names are used as 
#' gene names
#' @param job.path GMSCA generates a primary network and as many secondary networks (and expression
#' matrices) as cell types considered. Everything is stored in this folder
#' @param significanceThreshold 0.05 by default, it is the final p-value threshold to consider a module enriched
#' for gene markers
#' @param markers.path Where to find the gene marker sets 
#' @param celltypes Names of cell types to be used internally to tag files and cell types
#' @param cellsubtypes A list with marker set files for each cell type
#' @param tissue A string to be used to name files
#'
#' @return
#' @export
#'
#' @examples
createPSGCN = function(expr.data,
                       job.path,
                       significanceThreshold = 0.05, 
                       markers.path = system.file("markers", "",
                                                  package = "GPCNA"),
                       celltypes = c("N","MG","OLG","MA"), 
                       cellsubtypes = getCellSubtypes(), 
                       tissue = NULL,
                       ...){
  
  stopifnot(!is.null(job.path))
  stopifnot(!is.null(tissue))
  stopifnot(typeof(tissue) == "character")
  stopifnot(dir.exists(job.path))
  stopifnot(length(grep("^ENSG",colnames(expr.data))) == length(colnames(expr.data)))
  
  
  pnetname = paste0(job.path,"/primaryNet",tissue,".rds")
  pdataname = paste0(job.path,"/primaryExprData",tissue,".rds")
  
  if(typeof(expr.data) == "character"){
    expr.dataf = expr.data
    expr.data = readRDS(expr.dataf)
  }
  
  net = CoExpNets::getDownstreamNetwork(expr.data=expr.data,
                             job.path = NULL,
                             ...)  
  saveRDS(net,pnetname)
  saveRDS(expr.data,pdataname)
  
  createSGCN(pdataname,pnetname,...)
}


createSGCN = function( original.data.filename, 
                       primary.net.filename,
                       significanceThreshold = 0.05, 
                       markers.path = system.file("markers", "",
                                                  package = "GPCNA"),
                       celltypes = c("N","MG","OLG","MA"), 
                       cellsubtypes = getCellSubtypes(), 
                       singleSignal = F,
                       tissue = NULL,
                       ...) {
  
  stopifnot(file.exists(original.data.filename))
  stopifnot(file.exists(primary.net.filename))
  
  if ( is.null( tissue ) ) {
    tissue = unlist(strsplit(original.data.filename, "[.]"))[1]
  }
  original.data = readRDS( original.data.filename )
  primary.net = readRDS( primary.net.filename )
  primary.enrichment.by.module = getModulesEnrichment( net = primary.net, 
                                                       markers.path = markers.path,
                                                       significanceThreshold = 1)
  
  if(sum(primary.enrichment.by.module) == 0){
    warning("No secundary network can be created. No cell enrichment detected")
    return(NULL)
  }
  primary.enrichment.names = rownames(primary.enrichment.by.module)
  
  primary.module.enrichment = aggModEnrichmentConf( enrichment.by.module = primary.enrichment.by.module,
                                                    celltypes = celltypes,
                                                    cellsubtypes = cellsubtypes,
                                                    significanceThreshold = significanceThreshold,
                                                    singleSignal = singleSignal, 
                                                    as.is = F)
  
  cts = intersect(celltypes,primary.module.enrichment)
  for(ct in cts){
    modules.to.clean = names(primary.module.enrichment)[which(primary.module.enrichment == ct)]
    print("Modules to clean");
    print(modules.to.clean)
    cleaned.data = NULL
    if ( length(modules.to.clean) > 0 ) {
      cleaned.data = removePrimaryEffect( expr.data = original.data.filename, 
                                          target.enrichment = NULL, 
                                          net = primary.net.filename, 
                                          significanceThreshold = significanceThreshold, 
                                          modules.to.clean = modules.to.clean, 
                                          markers.path = markers.path )
    }
    if ( is.null(cleaned.data) ) {
      print("Nothing done creating secondary net")
    } else {
      secondary.data.filename = paste0( original.data.filename, ".No.", ct, ".rds" )
      secondary.net.filename = paste0( primary.net.filename, ".No.", ct, ".rds" )
      saveRDS( cleaned.data, secondary.data.filename )
      print( paste0( "Generating secondary co-expression network from removing ", ct ) )
      secondary.net = CoExpNets::getDownstreamNetwork(expr.data=cleaned.data,
                                           job.path = NULL,
                                           ...) 
      saveRDS( secondary.net, secondary.net.filename )
    }  
  }
  getPredictions(primary.net.filename=primary.net.filename,
                 original.data.filename=original.data.filename,
                 markers.path = markers.path,
                 significanceThreshold = significanceThreshold,
                 celltypes = celltypes)
  
  
}

getPredictions = function(primary.net.filename,
                          original.data.filename,
                          markers.path = system.file("markers", "",
                                                     package = "GPCNA"),
                          significanceThreshold = 0.05, 
                          celltypes=c("N","MG","OLG","MA"), 
                          cellsubtypes = getCellSubtypes()){
  supertabla = data.frame()
  original.data = readRDS(original.data.filename)
  primary.net = readRDS(primary.net.filename)
  for ( src.ct in c("-",celltypes) ) {
    
    print( paste0( "Procesing: ", src.ct ) )
    if(src.ct == "-"){
      for(localct in celltypes){
        print( paste0( "Subprocesing: ", localct ) )
        secondary.net.filename = paste0( primary.net.filename, ".No.", localct, ".rds" )
        secondary.data.filename = paste0( original.data.filename, ".No.", localct, ".rds" )
        secondary.net = readRDS( secondary.net.filename )
        cleaned.data = readRDS( secondary.data.filename )
        variance.changes = varianceEvolution( original.data, cleaned.data, genes = NULL )
        
        secondary.enrichment.by.module = getModulesEnrichment( net = secondary.net,
                                                               markers.path = markers.path,
                                                               outputCorrectedPvalues = T)
        secondary.enrichment.names = rownames(secondary.enrichment.by.module)
        
        primary.module.enrichment = summaryModuleEnrichment(primary.net.filename = primary.net.filename,
                                                            significanceThreshold = significanceThreshold,
                                                            markers.path = markers.path,
                                                            celltypes = celltypes,
                                                            cellsubtypes = cellsubtypes)
        
        modules.to.clean = names(primary.module.enrichment)[which(primary.module.enrichment == src.ct)]
        print( "modules.to.clean: " )
        print( modules.to.clean )
        processed.genes = names(primary.net$moduleColors[primary.net$moduleColors %in% modules.to.clean])
        if ( length( processed.genes ) == 0 ) next
        secondary.module.enrichment = aggModEnrichmentConf( enrichment.by.module = secondary.enrichment.by.module,
                                                            celltypes = celltypes,
                                                            cellsubtypes = cellsubtypes,
                                                            significanceThreshold = significanceThreshold,
                                                            singleSignal = F, F)
        
        processed.genes.celltype.changes = data.frame( gene = processed.genes,
                                                       primary.module = as.character(primary.net$moduleColors[processed.genes]),
                                                       primary.enrichment = src.ct,
                                                       secondary.module = as.character(secondary.net$moduleColors[processed.genes]),
                                                       secondary.enrichment = as.character(secondary.module.enrichment[secondary.net$moduleColors[processed.genes]]),
                                                       sec.net = secondary.net.filename,
                                                       stringsAsFactors=F,
                                                       row.names=NULL )
        
        
        processed.genes.variance.changes = variance.changes[variance.changes$gene %in% processed.genes, ]
        markers = getMarkerGenes( markers.path )
        crv = 1 - processed.genes.variance.changes$cleaned.variance / processed.genes.variance.changes$original.variance
        
        tabla = cbind( processed.genes.celltype.changes, crv, processed.genes.variance.changes[, 2:3])
        tabla$third = ifelse( crv < .33, 1, ifelse( crv >= 0.33 & crv <= 0.66, 2, 3) )
        # type 1 es que ha cambiado de tipo de célula, type 2 es que se mantiene, type 3 es que se desactiva
        tabla$type = ifelse( tabla$primary.enrichment != tabla$secondary.enrichment & tabla$secondary.enrichment != "-", 1,
                             ifelse( tabla$primary.enrichment == tabla$secondary.enrichment & tabla$secondary.enrichment != "-", 2,
                                     ifelse(tabla$primary.enrichment != tabla$secondary.enrichment & tabla$secondary.enrichment == "-", 3, 0)) )
        target = c(src.ct)
        if ( !is.null( cellsubtypes ) ) {
          if ( exists( src.ct, cellsubtypes ) ) 
            target = as.vector(unlist(cellsubtypes[src.ct]))
        }
        tabla$marker = tabla$gene %in% unlist(markers[target])
        supertabla = rbind( supertabla, tabla )
      }
      
    }else{
      secondary.net.filename = paste0( primary.net.filename, ".No.", src.ct, ".rds" )
      secondary.data.filename = paste0( original.data.filename, ".No.", src.ct, ".rds" )
      if(file.exists(secondary.data.filename) & file.exists(secondary.net.filename)){
        secondary.net = readRDS( secondary.net.filename )
        cleaned.data = readRDS( secondary.data.filename )
        
        secondary.enrichment.by.module = getModulesEnrichment( net = secondary.net,
                                                               markers.path = markers.path,
                                                               outputCorrectedPvalues = T)
        secondary.enrichment.names = rownames(secondary.enrichment.by.module)
        
        #    celltype.changes = enrichmentEvolution( primary.net = primary.net, secondary.net = secondary.net, genes = NULL,  significanceThreshold = significanceThreshold, markers.path = markers.path )
        variance.changes = varianceEvolution( original.data, cleaned.data, genes = NULL )
        
        primary.module.enrichment = summaryModuleEnrichment(primary.net.filename = primary.net.filename,
                                                            significanceThreshold = significanceThreshold,
                                                            markers.path = markers.path,
                                                            celltypes = celltypes,
                                                            cellsubtypes = cellsubtypes)
        
        modules.to.clean = names(primary.module.enrichment)[which(primary.module.enrichment == src.ct)]
        print( "modules.to.clean: " )
        print( modules.to.clean )
        processed.genes = names(primary.net$moduleColors[primary.net$moduleColors %in% modules.to.clean])
        if ( length( processed.genes ) == 0 ) next
        secondary.module.enrichment = aggModEnrichmentConf( enrichment.by.module = secondary.enrichment.by.module,
                                                            celltypes = celltypes,
                                                            cellsubtypes = cellsubtypes,
                                                            significanceThreshold = significanceThreshold,
                                                            singleSignal = F, F)
        processed.genes.celltype.changes = data.frame( gene = processed.genes,
                                                       primary.module = as.character(primary.net$moduleColors[processed.genes]),
                                                       primary.enrichment = src.ct,
                                                       secondary.module = as.character(secondary.net$moduleColors[processed.genes]),
                                                       secondary.enrichment = as.character(secondary.module.enrichment[secondary.net$moduleColors[processed.genes]]),
                                                       sec.net = secondary.net.filename,
                                                       stringsAsFactors=F,
                                                       row.names=NULL )
        
        
        processed.genes.variance.changes = variance.changes[variance.changes$gene %in% processed.genes, ]
        
        markers = getMarkerGenes( markers.path )
        
        crv = 1 - processed.genes.variance.changes$cleaned.variance / processed.genes.variance.changes$original.variance
        
        tabla = cbind( processed.genes.celltype.changes, crv, processed.genes.variance.changes[, 2:3])
        tabla$third = ifelse( crv < .33, 1, ifelse( crv >= 0.33 & crv <= 0.66, 2, 3) )
        # type 1 es que ha cambiado de tipo de célula, type 2 es que se mantiene, type 3 es que se desactiva
        tabla$type = ifelse( tabla$primary.enrichment != tabla$secondary.enrichment & tabla$secondary.enrichment != "-", 1,
                             ifelse( tabla$primary.enrichment == tabla$secondary.enrichment & tabla$secondary.enrichment != "-", 2,
                                     ifelse(tabla$primary.enrichment != tabla$secondary.enrichment & tabla$secondary.enrichment == "-", 3, 0)) )
        target = c(src.ct)
        if ( !is.null( cellsubtypes ) ) {
          if ( exists( src.ct, cellsubtypes ) ) {
            target = as.vector(unlist(cellsubtypes[src.ct]))
          }
        }
        tabla$marker = tabla$gene %in% unlist(markers[target])
        supertabla = rbind( supertabla, tabla )
      }
      
    }
  }
  
  supertabla = supertabla[!(supertabla$primary.enrichment == "-" & supertabla$secondary.enrichment == "-"),]
  return(supertabla)
}

summaryModuleEnrichment = function(primary.net.filename,
                                   significanceThreshold = 0.05, 
                                   markers.path = system.file("markers", "",
                                                              package = "GPCNA"),
                                   celltypes = c("N","MG","OLG","MA"),
                                   singleSignal = F,
                                   as.is = F,
                                   cellsubtypes = getCellSubtypes()){
  
  primary.net = readRDS(primary.net.filename)
  primary.enrichment.by.module = getModulesEnrichment( net = primary.net, 
                                                       markers.path = markers.path,
                                                       significanceThreshold = 1)
  
  if(sum(primary.enrichment.by.module) == 0){
    warning("No secundary network can be created. No cell enrichment detected")
    return(NULL)
  }
  aggModEnrichmentConf( enrichment.by.module = primary.enrichment.by.module,
                        celltypes = celltypes,
                        cellsubtypes = cellsubtypes,
                        significanceThreshold = significanceThreshold,
                        singleSignal = singleSignal,
                        as.is = as.is)
}
