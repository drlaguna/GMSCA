library(WGCNA)
library(flashClust)
library(dynamicTreeCut)

#' Best Coorrelated Module Selection
#'
#' The function selects the best module for a gene on the basis of the highest
#' correlation of that gene's expression with the eigengenes
#' @param gene numerical vector with the gene expression
#' @param centroids Matrix with many rows as \code{gene} components and a column for each module eigengene
#' @param signed is the network signed?
#' @return The index of the eigengene within the matrix passed as argument
getBestModuleCor <- function(gene,centroids,signed=TRUE){
  return(which.max(corDistance(centroids,gene,signed)))
}

createCentroidMatrix <- function(eigengenes){
  my.matrix <- NULL
  for(eigengene in eigengenes){
    my.matrix <- cbind(my.matrix,eigengene)
  }
  return(my.matrix)
}

getNewCentroids <- function(expr.data,partition.in.colors,centroid.labels,mgg){
  if(sum(partition.in.colors == "grey") < mgg)
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)$eigengenes
  else
    eg.vectors = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)$eigengenes

  names(eg.vectors) <- substring(names(eg.vectors),3)
  eg.vectors <- eg.vectors[,centroid.labels]
  return(eg.vectors)
}

getExchangedGenes <- function(old.partition,new.partition){
  stopifnot(length(old.partition) == length(new.partition))
  return(old.partition[old.partition != new.partition])
}

#' Distance based on correlation
#'
#' Function that generates the correlation for distance, in line
#' with a signed network
#' @param a numerical vector with gene expression for gene a
#' @param b numerical vector with gene expression for gene a
#' @param signed is the network signed?
#' @return The signed correlation if signed is TRUE
corDistance = function(a,b,signed=TRUE){
  if(signed)
    return(0.5 * (1 + stats::cor(a,b)))
  return(abs(stats::cor(a,b)))
}

#' Apply K-means algoritm on a GCN
#'
#' This function takes a hopefully correctly created WGCN network and, based on gene expression it
#' modifies the module color assignment by using k-means heuristic.
#' The description of the parameters show how the function works.
#' @param net must be a RDS object storing the network. Networks can be created using \code{\link{createGCN}} function.
#' @param expr.data can be a full path file name or a data frame with genes in columns and samples in
#' rows, giving the expression data used to construct the WGCNA network. The function expects genes in the
#' column order to correspond to the order of nameColors of net.file param.
#' @param  beta The soft threshold parameter as it was used with WGCNA
#' @param n.iterations A max number of iterations to run the k-means
#' @param net.type it is used in the same sense as WGCNA
#' @return The network, post-processed with a k-means heuristic including two additional properties
#' \enumerate{
#'   \item iterations The number of iterations that the algorithm actually made.
#'   \item exchanged.genes The number of exchanged genes in the last iteration of the algorithm.
#' }#'
applyKMeans <- function(net,
                        expr.data,
                        beta=6,
                        n.iterations=20,
                        net.type="signed" ){


  min.genes.for.grey = 20

  #ALGORITHM SPECIFICATION
  #Step 1. Let D be the expression data in which dij in D represents the expression value for
  #sample i and gene j, being s samples and g genes in total.
  #Step 2. Construct the partition by the WGCNA process, let P_D={m_1, m_2, ..., m_n} be
  #that partition where m_k is the k-th module.
  #Step 3. Get the eigengenes for each module within the partition, E={e_1, e_2, ..., e_n}
  #Step 4. Set up the k-means clustering
  #Step 4.1. Set k to n
  #Step 4.2. Set the centroids C to the eigengenes E, thus C to E
  #Step 5. Run the algorithm and monitor its evolution
  #Step 5.1 Set iterations to 0
  #Step 5.2 Create a new partition P', given C with n modules such that, for each gene, 1 <=
  #		j <= g, g_j belongs to the module c_t in C such that a distance meassure d(g_j,c_t) is
  #		minimum.
  #		Step 5.3 Calculate eigengenes of P', giving a new E'
  #		Step 5.4 Evaluate the progress. If progress done, set iterations to iterations + 1 and
  #		C to E' and go to step 5.2
  #Step 5.5 Finish
  #

  #Gather the current partition we start from
  partition.in.colors <- net$moduleColors

  if(sum(partition.in.colors == "grey") < min.genes.for.grey)
    eigengenes = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=TRUE)
  else
    eigengenes = moduleEigengenes(expr.data,partition.in.colors, excludeGrey=F)

  #This variable is fixed and used as a reference to indicate the
  #modules used (they are colours but the position within the vector is
  #also relevant)
  centroid.labels <- substring(names(eigengenes$eigengenes),3)

  k <- length(eigengenes$eigengenes)
  #Centroids must be a matrix with as much colums as centroids,
  #as much rows as samples
  centroids <- createCentroidMatrix(eigengenes$eigengenes)

  #Step 5
  #For storage of all the partitions created during the iterations
  partitions <- list()
  #A partition will be a list of as much elements as genes and for the
  #i-th position it stores the index of the module the ith gene belongs
  #to, and the color can be found in "centroid.labels"
  new.partition <- match(partition.in.colors, centroid.labels)
  names(new.partition) <- centroid.labels[new.partition]
  partitions[[1]] <- new.partition


  #Launch the iterations
  meg = 20 # Minimum exchanged genes
  exchanged.genes = meg + 1
  iteration = 1

  while(exchanged.genes > meg & iteration <= n.iterations){

    new.partition <- apply(expr.data,MARGIN=2,getBestModuleCor,centroids=centroids, signed=(net.type == "signed"))

    partitions[[iteration + 1]] <- new.partition
    #Get the control values for the new partition
    exchanged.gene.count <- length(getExchangedGenes(partitions[[iteration]],
                                                     partitions[[iteration + 1]]))

    new.partition.in.colors <- centroid.labels[unlist(new.partition)]
    centroids <- getNewCentroids(expr.data, new.partition.in.colors, centroid.labels, min.genes.for.grey)

    exchanged.genes = exchanged.gene.count
    iteration = iteration + 1
  }

  # Extract network from partitions
  index = length(partitions)

  colors = unique(names(partitions[[1]]))
  col.number = unique(partitions[[1]])


  new.net <- NULL
  new.net$beta = beta
  new.net$type = net.type
  if(index == 1){
    new.net$moduleColors = names(partitions[[index]])
  }else{
    new.net$moduleColors = colors[match(partitions[[index]],col.number)]
  }
  new.net$moduleLabels = partitions[[index]]
  names(new.net$moduleColors) = colnames(expr.data)
  names(new.net$moduleLabels) = colnames(expr.data)

  #If there are some grey genes as NA, add them again
  new.net$moduleColors[is.na(new.net$moduleColors)] = "grey"

  if(sum(new.net$moduleColors == "grey") >= meg)
    new.net$MEs <- moduleEigengenes(expr.data,new.net$moduleColors,softPower=beta, excludeGrey=F)$eigengenes
  else
    new.net$MEs <- moduleEigengenes(expr.data,new.net$moduleColors,softPower=beta, excludeGrey=T)$eigengenes

  new.net$iterations = iteration-1
  new.net$exchanged.genes = exchanged.genes
  return(new.net)
}

#' Gene Coexpression Network Creation
#'
#' This function creates a WGCN network using WGCNA with the additonal step described in https://doi.org/10.1186/s12918-017-0420-6.
#' @param expr.data can be a full path file name or a data frame with genes in columns and samples in
#' rows, giving the expression data used to construct the WGCNA network.
#' @param net.type it is used in the same sense as WGCNA.
#' @param min.cluster.size The minimum number of genes used during the dynamic tree cut phase en WGCNA pipeline.
#' @param n.iterations A max number of iterations to run the k-means algorithm.
#' @param save.adj.filename The name of the file in which to save the adjacency matrix generated by WGCNA.
#' If you do not want to save this matrix does just leave this parameter with its default value.
#' @param save.tom.filename The name of the file in which to save the TOM matrix generated by WGCNA.
#' If you do not want to save this matrix does just leave this parameter with its default value.
#' @return The network, post-processed with a k-means heuristic. The object returned has the following properties:
#' \enumerate{
#'   \item moduleColors a named vector whose elements are colors as assigned by WGCNA. The names() function on this
#'   propertie should return the gene names for each vector's element. Notice that the order of genes in moduleColors
#'   is the same as they appear in columns in expr.data parameter's object.
#'   \item moduleLabels contains a list with the gene names in the same order as in moduleColors.
#'   \item MEs is a matrix with a column for each module and a row for each sample. It contains the module eigengenes.
#'   \item MM contains a list with the module membership of each gene in the same order as in moduleColors.
#'   \item beta The soft threshold parameter used to create the adjacency matrix in WGCNA.
#'   \item type it is used in the same sense as WGCNA.
#'   \item iterations The number of iterations that the algorithm actually made.
#'   \item exchanged.genes The number of exchanged genes in the last iteration of the algorithm.
#' }#'
#' @export
createGCN <- function( expr.data,
                       net.type="signed",
                       min.cluster.size=30,
                       n.iterations = 30,
                       save.adj.filename=NULL,
                       save.tom.filename=NULL ){

  max.k.cutoff = 150
  r.sq.cutoff = 0.8
  net <- NULL

  if(typeof(expr.data) == "character") {
    expr.data = readRDS(expr.data)
  }

  stopifnot( net.type == "signed" | net.type == "unsigned" )

  #We assume gene names are at columns
  gene.names <- colnames(expr.data)
  sample.names <- rownames(expr.data)

  # We use pearson correlation and signed networks
  b.study = pickSoftThreshold(expr.data, powerVector=c(1:30), moreNetworkConcepts=TRUE, networkType=net.type, corOptions=list(use='p') )

  # Choosing beta from SFT.R.sq cut-off and max.k cut-off
  beta = Inf
  tomin = b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff & as.numeric(b.study$fitIndices$slope) < 0 & as.numeric(b.study$fitIndices$max.k) > max.k.cutoff,"Power"]
  if ( length(tomin) != 0 ) {
    beta = min(tomin)
  }

  if(beta == Inf){
    #OK, lets drop-off the maxk.cutoff
    beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff & as.numeric(b.study$fitIndices$slope) < 0,"Power"])
  }

  if(beta == -Inf & is.na(b.study$powerEstimate)){
    stop("There is something wrong, min beta is ", beta, " and suggested is", b.study$powerEstimate,"\n")
  }

  if(!is.na(b.study$powerEstimate)){
    if (beta == -Inf){
      cat("We'll use WGCNA's suggested beta\n")
      beta = b.study$powerEstimate
    } else if (beta - b.study$powerEstimate > 5){
      beta = trunc(0.5*(beta + b.study$powerEstimate))
      cat("We'll use average between WGCNA's suggested beta and ours.\n")
    }
  }

  if (beta == Inf) {
    beta = 21
    cat("Warning, the final beta is",beta,"and SFT is compromised\n")
  }

  cat("The final beta value to use is:",beta,"\n")


  #Create the adjacency matrix of genes coming from the expression data, with the beta
  #passwd as argument
  print("Creating adjacency matrix")
  adj = adjacency(expr.data, power = beta, type = net.type, corOptions = "use = 'p'")
  print("Created")
  if ( !is.null( save.adj.filename ) ) {
    saveRDS( adj, save.adj.filename )
  }
  # Topological Overlap Matrix (TOM)
  # Turn adjacency into topological overlap
  print("Creating TOM")
  if(net.type == "signed") {
    TOM = TOMsimilarity( adj )
  } else {
    TOM = TOMsimilarity( adj, TOMType="signed")
  }
  if ( !is.null( save.tom.filename ) ) {
    saveRDS( TOM, save.tom.filename )
  }
  print("Created")
  #Now we can delete adjacency
  dissTOM = 1-TOM

  print("Deleting adjacency matrix")
  rm(adj)
  gc()



  # Clustering using TOM
  # Call the hierarchical clustering function that makes the clustering
  # dendrogram to grow until the leaves are genes

  geneTree = flashClust(as.dist(dissTOM), method = "average")

  # Dynamic Tree Cut
  # We like large modules, so we set the minimum module size relatively high
  # Module identification using dynamic tree cut

  n.mods = 0
  deep.split = 2
  while(n.mods < 10 & deep.split < 5){
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = deep.split,
                                pamRespectsDendro = FALSE, minClusterSize = min.cluster.size)
    n.mods = length(table(dynamicMods))
    deep.split = deep.split + 1
  }

  print("Deleting TOM and dissTOM")
  rm( dissTOM )
  rm( TOM )
  gc()

  # Convert numeric lables into colors
  #This will print the same, but using as label for modules the corresponding colors
  dynamicColors = labels2colors(dynamicMods)

  # Calculate eigengenes
  MEList = moduleEigengenes(expr.data, colors = dynamicColors, excludeGrey=TRUE)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss), method = "average")

  MEDissThres = 0.1 #### MERGING THRESHOLD
  # Call an automatic merging function
  merge = mergeCloseModules(expr.data, dynamicColors, cutHeight = MEDissThres, verbose = 0, unassdColor="grey",getNewUnassdME = FALSE)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules
  mergedMEs = merge$newMEs

  #Prepare for creating the network objecto to return
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(400))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs

  net$MEs <- MEs
  rownames(net$MEs) <- sample.names
  net$moduleLabels <- moduleLabels
  net$moduleColors <- moduleColors
  net$beta = beta
  net$type = net.type
  names(net$moduleColors) <- gene.names
  net$geneTree <- geneTree
  print( "Computing Module Membership" )
  outnet = applyKMeans( net=net, expr.data=expr.data, beta=beta, n.iterations, net.type = net.type)
  gene.names = names(outnet$moduleColors)
  if ( sum(duplicated(gene.names)) > 0 ) {
    stop("There is some gene names duplicated, please remove this befor trying to create a Coexpression Network\n")
  }
  mm = NULL
  for ( gene in gene.names ) {
    module = outnet$moduleColors[ gene.names == gene ]
    mm[[gene]] = cor(outnet$MEs[paste0("ME",module)], expr.data[,gene])
  }
  outnet$MM = mm;
  return(outnet)
}


#' Modules enrichment
#'
#' This function determine the module enrichment of a GCN using WGCNA::userListEnrichment
#' @param net A GCN created with \code{\link{createGCN}}. It can be an R object or the full path of a file
#' containing it. If a filename is provided this function will first look for a file with the same name but
#' ending in ".enrichment.csv". If the file exists the function will return its content.
#' @param markers.path Folder containing user-defined lists of genes to be used as marker genes to determine
#' modules enrichment. This is done using WGCNA::userListEnrichment function so they must be in a compatible
#' format. Gene IDs must be expresed using the same format as in the network specified in the first parameter.
#' @param return.processed When True the functions returns the -log10 of the p-value obtained replacing the -inf
#' values by the highest value obtained. When False the function returns the p-values.
#' @return a matrix with a column for each module in the network and a row for each enrichment type. Values in the matrix
#' reflects the p-value resulting from the test made to determine if the module is enriched or not.
#' @export
getModulesEnrichment = function(net,
                                markers.path = ".",
                                return.processed=F){

  if ( typeof(net) == "character" ) {
    enrichment.filename = paste0( net, ".enrichment.csv")
    if ( file.exists( enrichment.filename ) ) {
      enrichment.by.module = read.csv( enrichment.filename, stringsAsFactors=F)
      return (enrichment.by.module)
    }
    net = readRDS(net)
  }
  modules = unique(net$moduleColors)

  files = list.files(path=markers.path,full.names = T)
  files = files[grep(pattern=".txt$",files)]
  markernames = gsub(".txt","",basename(files))

  ctypes = markernames
  ctypedata = matrix(ncol=length(modules),nrow=length(files))
  ctypedata[,] = 1
  rownames(ctypedata) = ctypes
  colnames(ctypedata) = modules

  all.gene.names = names(net$moduleColors)

  tmp.enr.f = tempfile()
  en = userListEnrichment( all.gene.names, net$moduleColors, files, nameOut = tmp.enr.f)
  if(file.exists(tmp.enr.f)){
    enrichment = read.csv(tmp.enr.f, stringsAsFactors=F)
    file.remove(tmp.enr.f)
    for(i in 1:nrow(enrichment)){
      module = enrichment$InputCategories[i]
      for(j in 1:length(files)){
        if(grepl(ctypes[j],enrichment$UserDefinedCategories[i])) {
          category = ctypes[j]
          ctypedata[category,module] = enrichment$CorrectedPvalues[i]
          break
        }
      }
    }
  }
  if ( return.processed ) {
    ctypedata = ctypedata[,apply(ctypedata,2,function(x){ any(x < 1)}),drop=FALSE]
    ctypedata = -log10(ctypedata)
    ctypedata[is.infinite(ctypedata)] = max(ctypedata[!is.infinite(ctypedata)])
    ctypedata = ctypedata[order(rownames(ctypedata)),,drop=FALSE]
  }
  return (ctypedata)
}

#' Cleaning the primary signal
#'
#' This function remove the primary signal of the groups of genes whose primary enrichment is indicated in
#' the following parameter.
#' @param expr.data The expression data to be cleaned. It can be a full path file name or a data frame with
#' genes in columns and samples in rows.
#' @param target.enrichment A string with the kind of enrichment that should be removed. This name (or list of
#' names) should be one from the list of enrichment marker genes' files provided.
#' @param net A GCN created with \code{\link{createGCN}} from the same expression data. It can be a full path
#' of a file containing it or an R object. This network reflects the primary signal and it is used to determine
#' what modules need to be cleaned in order to find the secondary signal. If no network is provided it will
#' be created. Therefore, providing one is also a way to reduce the time spent by this function.
#' @param markers.path Folder containing user-defined lists of genes to be used as marker genes to determine
#' modules' enrichment. This is done using WGCNA::userListEnrichment function, so they must be in a compatible
#' format. Gene IDs must be expresed using the same format as in the expression data specified in the first parameter.
#' @return  the expression data filtered that, hopefully, show the secondary role of some genes in the same format
#' of the input expression data.
#' @export
removePrimaryEffect = function( expr.data, target.enrichment, net = NULL, markers.path = NULL ) {
  if ( typeof(expr.data) == "character" ) {
    expr.data = readRDS(expr.data)
  }
  net.filename = "."
  if ( is.null( net ) ) {
    net = createGCN( expr.data )
  } else {
    if ( typeof( net ) == "character" ) {
      net.filename = net
      net = readRDS( net.filename )
    }
  }
  enrichment.filename = paste0( net.filename, ".celltype.csv")
  if ( file.exists( enrichment.filename ) ) {
    enrichment.by.module = read.csv( enrichment.filename, stringsAsFactors=F)
  } else {
    if ( is.null( markers.path ) ) markers.path = dirname( net.filename )
    enrichment.by.module = getModulesEnrichment( net = net, markers.path = markers.path)
  }
  enrichment.names = enrichment.by.module[,1]
  enrichment.by.module = enrichment.by.module[,-1][,as.vector(apply(enrichment.by.module[,-1], 2, FUN = function(x) { sum( x <= 0.05)} ) == 1)]
  enriched.modules = apply( enrichment.by.module, 2, FUN = function(x) { enrichment.names[which(x <= 0.05)] } )
  modules.to.clean = enriched.modules[ enriched.modules != target.enrichment]
  cleaned.data = expr.data[,!(net$moduleColors %in% names(modules.to.clean))]

  for ( module in names(modules.to.clean) ) {
    ee = net$MEs[paste0("ME", module)]
    expr = expr.data[,net$moduleColors == module]
    print(paste0("Cleaning module ", module) )
    res = apply( expr, 2, function(y) { lm( y ~ ., data=ee )$residuals })
    cleaned.data = cbind(cleaned.data, res)
  }
  net$type = "signed"
  if ( !is.null(net$type)) net.type = net$type
  cleaned.net = createGCN( expr.data = cleaned.data, net.type = net.type )
  return (list( expr.data = cleaned.data, gcn = cleaned.net ))
}


#' Showing changes in enrichment between networks
#'
#' This function compares the enrichment of a set of genes in two networks. Usually, the primary and the secondary network.
#' @param primary.net A GCN created with \code{\link{createGCN}}. It can be a full path of a file containing it or an R object.
#' @param secondary.net A GCN created with \code{\link{createGCN}} from a expression data containing the same genes as the
#' previous network. It can be a full path of a file containing it or an R object.
#' @param genes Subset of genes to be compared. The default value is null and that means all genes will be used in the comparison.
#' @return A dataframe with a row for each gene and the module and enrichment in both networks.
#' @export
enrichmenEvolution = function( primary.net, secondary.net, genes = NULL ) {
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
  if ( is.null( genes ) ) {
    genes = names( primary.net$moduleColors );
  }

  primary.enrichment.by.module = getModulesEnrichment( net = primary.net, markers.path = markers.path)
  secondary.enrichment.by.module = getModulesEnrichment( net = secondary.net, markers.path = markers.path)

    # Nos quedamos por un lado con los nombres de las columnas ( tipos de celula )
  enrichment.names = enrichment.by.module[,1]
  # Nos quedamos sólo con las columnas  que tienen modulos significativa y exclusivamente enriquecidos por un único tipo de célula
  em = enrichment.by.module[,-1][,as.vector(apply(enrichment.by.module[,-1], 2, FUN = function(x) { sum( x <= 0.05)} ) == 1)]
  # la lista de módulos etiquetada por el tipo de célula en el que están significativa y exclusivamente enriquecidos
  em = apply( em, 2, FUN = function(x) { enrichment.names[which(x <= 0.05)] } )

  tabla = data.frame( gene = genes, module = as.character(net$moduleColors[genes]), celltype = as.character(em[net$moduleColors[genes]]), stringsAsFactors=F )
  tabla$celltype[is.na(tabla$celltype)] = unlist(apply(enrichment.by.module[,tabla$module[is.na(tabla$celltype)],drop=F ], 2, FUN=function(x) { if (sum(x<=0.05) == 0) return("-") else return (paste(enrichment.names[which(x <= 0.05)],"(",x[x<=0.05],")", collapse=", ")) }))
  return(tabla)
}

