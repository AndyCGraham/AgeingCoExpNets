
#Remove Genes with expression below a threshold in any group
getExp <- function(dataset, traits, threshold, group = "Group"){
  #Get vector of gene names passing filter in any group
  Genes = unlist(lapply(unique(traits[,colnames(traits) == group]), function(Group){
    #Find the group
    Filtered = dataset[,traits[,colnames(traits) == group] == Group]
    #Get genes expressed above threshold in this group
    return(rownames(Filtered[rowMeans(Filtered) > threshold, ]))
  }
  ))
  Clean = dataset[rownames(dataset) %in% Genes,]
  return(Clean)
}

#Remove Genes with expression below a threshold in any group
getExpSC <- function(dataset, traits, threshold, clusters){
  #Get vector of gene names passing filter in any cluster
  Genes = unlist(lapply(unique(traits[,colnames(traits) == clusters]), function(Cluster){
    #Find the cells in this cluster
    Filtered = dataset[,traits[,colnames(traits) == clusters] == Cluster]
    #Get genes expressed in an above threshold percentage of cells in this cluster
    return(rownames(Filtered[rowSums(Filtered > 0)  > (threshold * ncol(Filtered)), ]))
  }
  ))
  Clean = dataset[rownames(dataset) %in% Genes,]
  return(Clean)
}

CoExpPrepare = function(dataset, traits, varThresh, expThresh, batches = NULL, plots = FALSE, modelMatrix){
  femaleSamples = colnames(dataset);
  traitRows = match(femaleSamples, rownames(traits));
  traits <- as.data.frame(traits[traitRows, ])
  
  #exclude non-variant genes
  getVar = apply(dataset, 1, var)
  dataset <- dataset[getVar > varThresh & !is.na(getVar), ] 
  
  #Exclude lowly expressed genes
  dataset <- getExp(as.matrix(dataset), traits, expThresh)
  
  dataset = log2(dataset+1)
  
  if(!is.null(batches)){
    
    #Make plot output folder if not already present
    dir.create(plots)
    #Get covariate names from model matrix to see effect of batch correction on them
    Covs = colnames(modelMatrix)[-1]
    
    #We will use only batch for getting the MDS plot coloured per sample
    traits[,colnames(traits)==batches] = as.factor(traits[,colnames(traits)==batches])

    #MDS takes a while to compute so we will use only a maximum of 100 samples to represent 
    mask = if(ncol(dataset) > 100){
      sample(1:ncol(dataset),100) 
      } else { 
        1:ncol(dataset)
      }
    #Get one color for each batch category
    colors = rainbow(length(unique(traits[,colnames(traits)==batches])))
    finalcolors = colors[as.numeric(traits[,batches])[mask]]
    
    #Plot samples with IDs and colour with their corresponding batch
    #samples will be disposed in 2D such that the distance between samples reflects differences
    #in log2 expression
    png(filename = paste0(plots, "Samples by ",batches, ".png"))
    limma::plotMDS(dataset[,mask],col=finalcolors,
                   main=paste0("MDS using ",batches)) 
    legend("topright",fill=colors,
           legend=levels(traits[,batches]))
    dev.off()
    
    femaleSamples = colnames(dataset);
    traitRows = match(femaleSamples, rownames(traits));
    traits <- as.data.frame(traits[traitRows, ])
    
    png(filename = paste0(plots, "Prince plot samples by ",batches, ".png"))
    pcres = prince(as.matrix(dataset), as.data.frame(traits[,colnames(traits) %in% c(batches, Covs)]),top=20, imputeknn = TRUE)
    CoExpNets::princePlot(prince=pcres,main="All samples, batch uncorrected")
    dev.off()
    
    dataset = limma::removeBatchEffect(dataset, batch = traits[,colnames(traits)==batches], design = modelMatrix)
    
    png(filename = paste0(plots, "Prince plot corrected samples by ",batches, ".png"))
    pcres = prince(as.matrix(dataset), as.data.frame(traits[,colnames(traits) %in% c(batches, Covs)]),top=20, imputeknn = TRUE)
    CoExpNets::princePlot(prince=pcres,main="All samples, batch corrected")
    dev.off()
    
    
    #MDS takes a while to compute so we will use only 100 samples to represent instead of the total
    mask = sample(1:ncol(dataset),ncol(dataset))
    #Get one color for each batch category
    finalcolors = colors[as.numeric(traits[,batches])[mask]]
    #Plot samples with IDs and colour with their corresponding batch
    #samples will be disposed in 2D such that the distance between samples reflects differences
    #in log2 expression
    png(filename = paste0(plots, "Corrected samples by ",batches, ".png"))
    limma::plotMDS(dataset[,mask],col=finalcolors,
                   main=paste0("MDS using ",batches, " after correction"))
    legend("topright",fill=colors,
           legend=levels(traits[,batches]))
    dev.off()
  }
  
  #Make sure all samples are matching 
  femaleSamples = colnames(dataset);
  traitRows = match(femaleSamples, rownames(traits));
  traits <- as.data.frame(traits[traitRows, ])
  
  #Transpose
  PreparedData <- t(dataset)
  return(PreparedData)
}

##Prepare single cell data for coexpnets
CoExpPrepareSC = function(dataset, traits, varThresh, expThresh, clusters){
  
  #exclude non-variant genes
  getVar = apply(dataset, 1, var)
  dataset <- dataset[getVar > varThresh & !is.na(getVar), ]
  
  #Exclude lowly expressed genes
  dataset <- getExpSC(as.matrix(dataset), traits, expThresh, clusters = clusters)
  
  #Log and transpose
  dataset = t(log2(dataset+1))
  return(dataset)
  
}

##Get co-expression metrics and genes for a module
mms_TOM <- function(dataset, module, tissue, Beta){
  genes = CoExpNets::getGenesFromModule(which.one="new",
                                        tissue=tissue,
                                        module=module) #get genes from module
  mms = CoExpNets::getMM(which.one="new",tissue=tissue,
                         genes= genes,expr.data.file = dataset) #Get MM of genes
  mmsO <- mms[order(mms$mm, decreasing = TRUE),] #order by MM
  moduleExpr <- dataset[,colnames(dataset) %in% genes]
  createTOM = function(expr.data.file,
                       beta=Beta,
                       save.as=NULL,
                       net.type="signed",
                       debug=F){
    
    stopifnot(beta > 0 & beta < 40)
    stopifnot(net.type == "signed" | net.type == "unsigned")
    
    if(typeof(expr.data.file) == "character"){
      print(paste0("Creating matrix ",save.as," from expression data ",expr.data.file))
      expr.data <- readRDS(expr.data.file)
    }else{
      expr.data <- expr.data.file
    }
    cat("Creating TOM for",ncol(expr.data),"genes and",nrow(expr.data),"samples, beta",
        beta,"and type",net.type,"\n")
    if(debug)
      expr.data = expr.data[,1:1000]
    adjacency = adjacency(expr.data, power = beta, type = net.type)
    print("Adjacency matrix created")
    # Topological Overlap Matrix (TOM)
    # Turn adjacency into topological overlap
    print("Creating TOM")
    TOM = TOMsimilarity(adjacency)
    colnames(TOM) = colnames(expr.data)
    rownames(TOM) = colnames(TOM)
    
    if(!is.null(save.as)){
      cat("Saving TOM at",save.as,"\n")
      saveRDS(TOM,save.as)
    }
    else return(TOM)
  }
  TOM <- as.data.frame(CoExpNets::createTOM(expr.data.file = moduleExpr, beta = Beta))
  TOM$from <- rownames(TOM)
  TOM <- gather(TOM, from)
  genes <- genes[genes %in% TOM[,1]]
  TOM$to <- rep(genes)
  TOM$Zero <- rep(0)
  TOM$MOO39 <- rep("M0039")
  TOM <- TOM[,c(1, 3, 4, 5, 2)]
  MMS_tom <- list(mmsO, TOM)
  return(MMS_tom)
}


corWithCatTraitsSC = function(tissue,which.one,covlist,covs=NULL,retPVals=T){
  if(is.null(covs))
    covs = getCovariates(tissue=tissue,which.one=which.one)
  
  if(!is.null(covlist))
    covs = covs[,covlist,drop=F]
  
  for(i in 1:ncol(covs)){
    if(typeof(covs[,i]) ==  "character")
      covs[,i] = as.factor(covs[,i])
  }
  factor.mask = unlist(lapply(covs,is.factor))
  cat("We will work with",sum(factor.mask),"factors\n")
  stopifnot(sum(factor.mask) > 0)
  
  MEs = getNetworkEigengenes(tissue=tissue,which.one=which.one)
  
  
  fcm = matrix(nrow=ncol(MEs),ncol=sum(factor.mask))
  index = 1
  for(i in which(factor.mask)){
    #cat("Factor",colnames(trait.data)[i],"\n")
    #cat("Levels",levels(trait.data[,i]))
    #print(trait.data[,i])
    for(j in 1:ncol(MEs)){
      #print(paste0(i,j)
      if(length(unique(covs[,i])) > 1){
        form = eg ~ cov
        data.in = data.frame(MEs[,j],covs[,i])
        colnames(data.in) = c("eg","cov")
        fcm[j,index] = anova(aov(form,data.in))$`Pr(>F)`[1]
      }else
        fcm[j,index] = 1
      
    }
    fcm[,index] = p.adjust(fcm[,index],method="BH")
    index = index + 1
  }
  
  
  
  if(sum(!factor.mask) > 0){
    moduleTraitCor = cor(MEs,covs[,!factor.mask,drop=FALSE], use="p")
    #Generate the p-values for significance of a given matrix of correlations, for all modules,
    #between traits data and eigengenes, both from samples
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nrow(MEs))
    moduleTraitPvalue = cbind(moduleTraitPvalue,fcm)
    colnames(moduleTraitPvalue) = c(colnames(covs)[!factor.mask],
                                    colnames(covs)[factor.mask])
  }else{
    moduleTraitPvalue = fcm
    colnames(moduleTraitPvalue) = colnames(covs)[factor.mask]
  }
  if(retPVals)
    toReturn = moduleTraitPvalue
  else
    toReturn = -log10(moduleTraitPvalue)
  rownames(toReturn) = gsub("ME","",names(MEs))
  moduleTraitPvalue = -log10(moduleTraitPvalue)
  moduleTraitPvalue[moduleTraitPvalue > 200] = 150
  WGCNA::labeledHeatmap(Matrix=moduleTraitPvalue,
                        
                        xLabels=colnames(moduleTraitPvalue),
                        yLabels=gsub("ME","",names(MEs)),
                        ySymbols=names(MEs),
                        colorLabels=TRUE,
                        colors=rev(heat.colors(50)),
                        cex.text=0.5,
                        zlim = c(0,100),
                        main="Module-trait relationships")
  return(toReturn)
}

preservationOneWayNew <- function(network,
                                  expr.data.files=NULL,
                                  tissues=c("snig","putm"),
                                  permutations=200,
                                  maxModuleSize=5000,
                                  maxGoldModuleSize=400,
                                  randomSeed=1){
  
  cat("Entering preservation\n")
  n1.shortname = tissues[1]
  n2.shortname = tissues[2]
  
  if(typeof(network) == "character"){
    print(paste0("Reading network ",network))
    network1 <- readRDS(network)
  }else{
    network1 = network
  }
  
  
  print( tissue1 <- tissues[1] )
  print( tissue2 <- tissues[2] )
  print(paste0(tissue1," vs. ", tissue2))
  
  if( tissues[1] == tissues[2] ) stop("Can't do a preservation against self")
  options(stringsAsFactors = FALSE)
  
  print(paste0("Reading expression data for tissue ",tissues[1]))
  if(typeof(expr.data.files[[1]]) == "character"){
    expression.data1 <- readRDS(expr.data.files[[1]])
    cat(expr.data.files[[1]],"\n")
    
  }else{
    expression.data1 = expr.data.files[[1]]
  }
  print(expression.data1[1:5,1:5])
  
  print(paste0("Reading expression data for tissue ",tissues[2]))
  if(typeof(expr.data.files[[2]]) == "character"){
    expression.data2 <- readRDS(expr.data.files[[2]])
    cat(expr.data.files[[2]],"\n")
  }else
    expression.data2 = expr.data.files[[2]]
  print(expression.data2[1:5,1:5])
  
  
  ## Prepare the data
  
  #First we check for good genes and samples
  intersect.g = intersect(colnames(expression.data1),colnames(expression.data2))
  expression.data1 = expression.data1[,match(intersect.g,colnames(expression.data1))]
  expression.data2 = expression.data2[,match(intersect.g,colnames(expression.data2))]
  
  network1 = network1$moduleColors[match(intersect.g,names(network1$moduleColors))]
  network2 = network1
  
  cat("We'll use",ncol(expression.data1),"genes for the preservation analysis\n")
  
  multiExpr <- list()
  multiExpr [[1]] <- list(data = expression.data1)
  multiExpr [[2]] <- list(data = expression.data2)
  
  names(multiExpr) <- c(tissues[1], tissues[2])
  
  checkSets(multiExpr, checkStructure = FALSE, useSets = NULL)
  multiColor <- list( network1) #, network2 )
  names(multiColor) <- c(tissues[1]) #, tissues[2])
  ## Run the preservation statistics and save
  enableWGCNAThreads()
  print( WGCNAnThreads() )
  system.time( {
    mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                             nPermutations = permutations,
                             networkType="signed",
                             maxGoldModuleSize=maxGoldModuleSize,
                             ## no. of permutations
                             randomSeed=randomSeed,
                             verbose=3,
                             maxModuleSize=maxModuleSize)
  })
  
  getPreservationStatisticsOneWay(tissues=tissues,
                                  presRes=mp)
  
}


getPreservationStatisticsOneWay <- function(tissues,presRes){
  
  Zsummary <- NULL
  MedianRank <- NULL
  
  tissue1 = paste0("ref.",tissues[1])
  tissue2 = paste0("inColumnsAlsoPresentIn.",tissues[2])
  Z.tmp.list <- list(NULL)
  MR.tmp.list <- list(NULL)
  cat(tissue1, "\t", tissue2, "\n")
  fn	<- paste(tissue1, "vs", tissue2,sep="")
  mp = presRes
  return(list(Z=as.data.frame(mp$preservation$Z[[tissue1]][[tissue2]],stringsAsFactors=F),
              logp=as.data.frame(mp$preservation$log.pBonf[[tissue1]][[tissue2]],stringsAsFactors=F)))
  
}

getZsummaryPress = function(tissues,presRes,module,statistic="Zsummary.pres"){
  pData = getPreservationStatisticsOneWay(tissues=tissues,presRes=presRes)$Z
  mask = rownames(pData) == module
  return(pData[mask,statistic])
}

getMeanZsummaryPress = function(tissue,tissues,
                                package="CoExp10UKBEC",
                                folder="micro19K",
                                module,
                                statistic="Zsummary.pres"){
  path = system.file(package = package)
  if(path == ""){
    cat("Can't find package ",package,"\n")
    return(NULL)
    
  }
  path = paste0(path,"/",folder)
  means = 0
  for(tother in tissues){
    #THAL.mic.net.19K.rds_vs_WHMT.mic.net.19K.rds_preserv.rds
    if(tissue < tother)
      fpres = paste0(path,"/",tissue,".mic.net.19K.rds_vs_",tother,".mic.net.19K.rds_preserv.rds")
    else
      fpres = paste0(path,"/",tother,".mic.net.19K.rds_vs_",tissue,".mic.net.19K.rds_preserv.rds")
    
    tlist = c(tissue,tother)
    cat("Reading from",fpres,"\n")
    partial = getZsummaryPress(tissues=tlist,
                               presRes=readRDS(fpres),
                               module=module,
                               statistic=statistic)
    print(partial)
    means = means + partial
  }
  return(means/length(tissues))
}

RankcorWithNumTraits = function(tissue,which.one,covlist,covs=NULL){
  MEs = getNetworkEigengenes(tissue=tissue,which.one=which.one)
  if(is.null(covs))
    covs = getCovariates(tissue=tissue,which.one=which.one)
  covs = covs[,covlist]
  moduleTraitCor = cor(MEs, covs, use = "pairwise.complete.obs", method = "spearman")
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(MEs))
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  if(length(covlist) == 1){
    textMatrix = cbind(rep("--",length(textMatrix)),textMatrix)
    moduleTraitCor = cbind(rep(0,nrow(moduleTraitCor)),moduleTraitCor)
    covlist = c("Dummy",covlist)
  }else
    dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  b = WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                            xLabels = covlist,
                            yLabels = names(MEs),
                            ySymbols = names(MEs),
                            colorLabels = FALSE,
                            colors = blueWhiteRed(50),
                            textMatrix = textMatrix,
                            setStdMargins = FALSE,
                            cex.text = 0.5,
                            zlim = c(-1,1),
                            main = paste0("Module-trait relationships"))
  print(b)
  return(moduleTraitCor)
}



NewcorWithNumTraits = function(tissue,which.one,covlist,covs=NULL){
  MEs = getNetworkEigengenes(tissue=tissue,which.one=which.one)
  if(is.null(covs))
    covs = getCovariates(tissue=tissue,which.one=which.one)
  covs = covs[,covlist]
  moduleTraitCor = cor(MEs, covs, use = "pairwise.complete.obs", method = "pearson")
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(MEs))
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  if(length(covlist) == 1){
    textMatrix = cbind(rep("--",length(textMatrix)),textMatrix)
    moduleTraitCor = cbind(rep(0,nrow(moduleTraitCor)),moduleTraitCor)
    covlist = c("Dummy",covlist)
  }else
    dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  b = WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                            xLabels = covlist,
                            yLabels = names(MEs),
                            ySymbols = names(MEs),
                            colorLabels = FALSE,
                            colors = blueWhiteRed(50),
                            textMatrix = textMatrix,
                            setStdMargins = FALSE,
                            cex.text = 0.5,
                            zlim = c(-1,1),
                            main = paste0("Module-trait relationships"))
  print(b)
  return(moduleTraitCor)
}


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Assess_exp = function(dataset, traits, Genes, title, var){
  require("ggpubr")
  
  #Make sure traits and metadata mach
  datTraits = as.data.frame(traits[match(rownames(dataset), rownames(traits)),])
  
  #Get gene expression of module genes
  datExpr0 = as.data.frame(t(dataset))
  Cyan <- as.data.frame(datExpr0[rownames(datExpr0) %in% Genes,])
  rownames(Cyan) <- Genes[Genes %in% rownames(datExpr0)]
  Ages <- as.data.frame(datTraits[rownames(datTraits) %in% colnames(datExpr0),])
  Ages$x <- rownames(datTraits) 
  
  #Get mean module gene expression
  Cyan <- as.data.frame(apply(Cyan, 2,as.numeric))
  Cyan <- as.data.frame(t(Cyan))  
  colnames(Cyan) <- Genes[Genes %in% rownames(datExpr0)]
  Cyan$mean <- rowMeans(Cyan)
  Cyan$x <- rownames(Cyan)
  Cyan <- merge(Cyan, Ages, by = "x")
  Cyan$mean <- as.numeric(Cyan$mean)
  Cyan = Cyan[!is.na(Cyan[,colnames(Cyan) == var]),] 
  Cyan$var = Cyan[,colnames(Cyan) == var]
  
  #Run ANOVA or T-test
  if(length(unique(Cyan[,colnames(Cyan) == var])) > 2){
    res.aov = aov(mean ~ var, Cyan)
    DTres = DunnettTest(Cyan$mean, Cyan$var)
    print(paste0(title, " ANOVA results:"))
    print(summary(res.aov))
    print(DTres)
  } else {
    res = t.test(mean ~ var, Cyan)
    print(paste0(title, " t-test results:"))
    print(res)
  }
  
  #Get scores for each group to plot in prism
  SUM = sapply(unique(Cyan[,colnames(Cyan) == var]), function(x){
    return(Cyan[Cyan$var == x,]$mean)
  })
  names(SUM) = unique(Cyan[,colnames(Cyan) == var])
  
  Cyan_sum <- data_summary(Cyan, varname="mean", 
                           groupnames=var)
  
  boxplt = ggboxplot(Cyan, x = var, y = "mean",
                color = var, palette = rainbow(n=length(unique(Cyan[,colnames(Cyan) == var]))),
                ylab = "Normalised Counts", xlab = var)
  print(boxplt)
  return(SUM)
}

GetModuleExpsforMonocle = function(module){
  turquoise = Assignments[Assignments$`net$moduleColors` == c(module),]
  turquoise = merge(turquoise, datExpr0, by = "Genes")
  turquoise[nrow(turquoise)+1,] = rep(NA)
  turquoise = turquoise[,-c(1,2)]
  turquoise[is.na(turquoise)] = 0
  b = colMeans(turquoise)
  turquoise[nrow(turquoise)+1,] = b
  point.data[,paste(module, "Network", sep = "")] = as.numeric(b)
  return(point.data)
}
