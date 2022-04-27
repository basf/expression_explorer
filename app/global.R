library(ggplot2)
library(ggrepel)
library(shiny)
library(reshape2)
library(ape)
library(DESeq2)
library(cyjShiny)
library(VennDiagram)
library(grid)
library(gridExtra)
library(shinyBS)
library(DT)
library(stringr)
library(readr)
library(trewjb)
library(shinyWidgets)

validate = shiny::validate

###read GO names
goDescriptions = read.table('go_names.tab',header=T,sep='\t',quote='',row.names=1)

# change baseFolder based on host
#baseFolder = ifelse(Sys.getenv("RSTUDIO") == "1", "./data", "/data")
# Based on how you mount this folder, you might want to change this
baseFolder='./data'

ds <- list.files(path=baseFolder, recursive = F, full.names = TRUE)
ds <- sapply(strsplit(ds,'/'), function(x) x[length(x)])
ds <- ds[grep(".tab", ds, invert = TRUE)]

info <- NULL
for (f in ds) {#Check for info of this dataset
  if (file.exists(paste(baseFolder,f,'info',sep='/'))) {
    i <- read.table(paste(baseFolder,f,'info',sep='/'),sep='\t')
    info <- c(info,as.character(i[1,1]))
  } else {
    info <- c(info,'No additional information available')
  }
}
allData <- data.frame(Datasets=ds,Info=info)  

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x,na.rm=T)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

getPositions <- function(click,data) {
  
  panel = click$panelvar1
  col <- click$y
  lvls <- levels(data$timepoint)
  row <- lvls[round(click$x)]
  
  return(which(data$value < col+0.1*col & data$value > col-0.1*col & data$timepoint==row & data$sample==panel))
}

getPoint <- function(click,data) {
  dd <- data[getPositions(click,data),]
  if (nrow(dd) > 1) {#Subselect to only the closest
    dif <- abs(dd$value-click$y)
    return(dd[which.min(dif),])
  }
  dd
}

### This function does a lot of heavy lifting when it comes to loading datasets!
# Probably too much lifting to be honest

loadDataset <- function(folder) {
  #Check for normalized counts
  folder = paste(baseFolder,folder,sep='/')
  meta <- read.table(paste(folder,'meta.tab',sep='/'),header=T,row.names=1,sep='\t',stringsAsFactors=F)
  
  if (!file.exists(paste(folder,'norm_counts.tab',sep='/'))) {
    if (!file.exists(paste(folder,'counts.tab',sep='/'))) {
      return("Could not find raw counts file. Contact your admin!")
    }
    data <- read.table(paste(folder,'counts.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
    mm <- data.frame(row.names=colnames(data),condition=colnames(data),type="paired-end")
    dds <- DESeqDataSetFromMatrix(countData=data,
                                  colData=mm,
                                  design=~condition)
    dds <- estimateSizeFactors(dds)
    norm <- counts(dds,normalized=TRUE)
    all_diffs = NULL
    for (i in 1:(ncol(norm)-1)) {
      for (j in (i+1):ncol(norm)) {
        main_diff = median(norm[,i]/norm[,j],na.rm=T)
        all_diffs = c(all_diffs,main_diff)
      }
    }
    write.table(norm,paste(folder,'norm_counts.tab',sep='/'),quote=F,sep='\t')
  }
  
  fpkm <- NULL
  genes <- NULL
  
  if (!file.exists(paste(folder,'fpkm_counts.tab',sep='/')) && file.exists(paste(folder,'genes.tab',sep='/'))) {
	  if (!file.exists(paste(folder,'counts.tab',sep='/'))) {
		return("Could not find raw counts file. Contact your admin!")
	  }
	  
		data <- read.table(paste(folder,'counts.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
		mm <- data.frame(row.names=colnames(data),condition=colnames(data),type="paired-end")
		dds <- DESeqDataSetFromMatrix(countData=data,
                                  colData=mm,
                                  design=~condition)
		dds <- estimateSizeFactors(dds)
	  
		genes <- read.table(paste(folder,'genes.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
		genes = genes[rownames(data),,F]
		if ('Length' %in% colnames(genes)) {
			genes$Start = 1
			genes$End = genes$Start + genes$Length
			genes$Strand = '*'
			#genes$chr = 'chrom'
		} 
		genesGR <- GRanges(genes)
		rowRanges(dds) <- genesGR
		res_fpkm <- fpkm(dds)
		write.table(res_fpkm,paste(folder,'fpkm_counts.tab',sep='/'),quote=F,sep='\t')
  } 
  
  
  
  if (file.exists(paste(folder,'fpkm_counts.tab',sep='/'))) {
		fpkm <- read.table(paste(folder,'fpkm_counts.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
		fpkm <- fpkm[order(rowSums(fpkm),decreasing=T),]
		genes <- read.table(paste(folder,'genes.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
  }
  
  span <- NULL
  if (file.exists(paste(folder,'span.tab',sep='/'))) {
    span <- read.table(paste(folder,'span.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
  }
  
  mutations <- NULL
  if (file.exists(paste(folder,'mutations.tab',sep='/'))) {
    mutations <- read.table(paste(folder,'mutations.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
  }
  
  prot <- NULL
  
  if (file.exists(paste(folder,'counts_prot.tab',sep='/'))) {
    prot <- read.table(paste(folder,'counts_prot.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
    prot <- prot[order(rowSums(prot),decreasing=T),]
  }
  
  norm <- read.table(paste(folder,'norm_counts.tab',sep='/'),header=T,row.names=1,sep='\t',check.names=F,quote='')
  norm <- norm[order(rowSums(norm),decreasing=T),]

  if(!file.exists(paste(folder,'edgeData.tab',sep='/'))){
    ######### Compute all correlations
    final <- sqrt(fpkm)
    final <- as.matrix(final[which(rowSums(final)>10),])
    c <- cor(t(final),method='pearson')
    
    #Make list of all paths in this
    is = (c >= 0.95 | c <= -0.95)
    keep = which(rowSums(is)>1)
    is <- is[keep,keep]
    is[lower.tri(is,diag=T)] <- FALSE
    c <- c[keep,keep]
  
    edgeData <- data.frame(reshape2::melt(c), stringsAsFactors=FALSE)
    k <- data.frame(reshape2::melt(is), stringsAsFactors=FALSE)
    print(identical(edgeData$Var1, k$Var1))
    edgeData$keep <- k$value # 5,631,129
    edgeData <- subset(edgeData, keep == TRUE) # 42785
    
    colnames(edgeData) <- c("source", "target", "score", "keep")
  
    write.table(edgeData[,1:3],paste(folder,'edgeData.tab',sep='/'),quote=F,sep='\t')
    
    id <- rownames(fpkm)
    name <- id
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    write.table(nodeData,paste(folder,'nodeData.tab',sep='/'),quote=F,sep='\t') 
  }

  edgeData <- read.table(paste(folder,'edgeData.tab',sep='/'),header=T,sep='\t',quote='')
  nodeData <- read.table(paste(folder,'nodeData.tab',sep='/'),header=T,sep='\t',quote='')
  
  if (file.exists(paste(folder,'interraction_network.tab',sep='/'))) {
    ref_interraction = read.table(paste(folder,'interraction_network.tab',sep='/'),header=T,sep='\t',quote='')
    ref_interraction$score = 0
    dd = ref_interraction[,c(1,3,4,2)]
    colnames(dd) = c('source','target','score','type')
    edgeData$type = 'correlation'
    #Add this to the edge data!
    edgeData = rbind(edgeData,dd)
  }
  
  ann <- read.table(paste(folder,'ann.tab',sep='/'),header=T,row.names=1,sep='\t',quote='',comment='',stringsAsFactors=F)
  if (!"synonym" %in% colnames(ann)) {
	ann$synonym = rownames(ann)
  }

  #### Read info about where the RNA-seq data is located. This will allow us to show the jbrowser
  jbrowse="";
  jbrowse_server="";
  if(file.exists(paste(folder,'jbrowse.tab',sep='/'))){
    jbrowse = read.table(paste(folder,'jbrowse.tab',sep='/'),sep='\t')
    jbrowse = jbrowse$V1
    jbrowse_server = jbrowse$V2
  }

  #### Read pre-saved gene lists
  geneSets = NULL
  if(file.exists(paste(folder,'genesets.tab',sep='/'))){
    geneSets = read.table(paste(folder,'genesets.tab',sep='/'),sep='\t',header=T)
  }
  
  sampleSets = NULL
  if(file.exists(paste(folder,'samplesets.tab',sep='/'))){
    sampleSets = read.table(paste(folder,'samplesets.tab',sep='/'),sep='\t',header=T)
  }
  
  diffLong <- list.files(folder, pattern = '*_diffExp.csv', full.names = TRUE)
  diffShort <- list.files(folder, pattern = '*_diffExp.csv', full.names = FALSE)
  
  diff <- list()
  diffPos = 1
  if(length(diffLong) > 0){
    for(i in 1:length(diffLong)){
    sets1 <- sapply(strsplit(diffShort[i], '_vs_'), '[[', 1)
    sets2 <- sapply(strsplit(diffShort[i], '_vs_'), '[[', 2)
    ss = unlist(strsplit(sets2, '-'))
    sets2 = paste(ss[1:(length(ss)-1)],collapse='-')
    dgs <- gsub('_diffExp.csv','',ss[length(ss)])
  
    diff[[diffPos]] <- list(path = diffLong[i], set1 = sets1, set2 = sets2, dgs = dgs)
    diffPos = diffPos + 1
    }} else{
    diff <- NULL
  }
  diffExpr = NULL
  if(file.exists(paste(folder,'diff_expr.tab',sep='/'))){
    diffExpr = read.table(paste(folder,'diff_expr.tab',sep='/'),sep='\t')
    for(i in 1:nrow(diffExpr)) {
      sets1 <- sapply(strsplit(diffExpr[i,2], '_vs_'), '[[', 1)
      sets2 <- sapply(strsplit(diffExpr[i,2], '_vs_'), '[[', 2)
      ss = unlist(strsplit(sets2, '-'))
      sets2 = paste(ss[1:(length(ss)-1)],collapse='-')
      dgs <- gsub('_diffExp.csv','',ss[length(ss)])
      
      diff[[diffPos]] <- list(path = paste(folder,diffExpr[i,1],sep='/'), set1 = sets1, set2 = sets2, dgs = dgs)
      diffPos = diffPos + 1
    }
  }

  res = list()
  res$meta <- meta
  res$norm <- norm
  res$ann <- ann
  #res$expression <- expression
  res$expression <- NULL
  res$edgeData <- edgeData
  res$nodeData <- nodeData
  res$c <- c
  res$fpkm <- fpkm
  res$diff <- diff
  res$folder <- folder
  res$jbrowse <- jbrowse
  res$jbrowse_server <- jbrowse_server
  res$prot <- prot
  res$span <- span
  res$mutations <- mutations
  res$geneSets <- geneSets
  res$sampleSets <- sampleSets
  res$geneCoordinates <- genes
  return(res)
  
}

makeDESeqObject <- function(data,set1,set2) {
  
  set1 <- set1[which(set1 %in% colnames(data))]
  set2 <- set2[which(set2 %in% colnames(data))]
  if (length(set1)<1 | length(set2)<1) {
    return(NULL)
  }
  final <- floor(data[,c(set1,set2)])
  ids = c(rep('set1',length(set1)),rep('set2',length(set2)))
  mm <- data.frame(row.names=colnames(final),condition=ids,type="paired-end")
  
  dds <- DESeqDataSetFromMatrix(countData=final,
                                colData=mm,
                                design=~condition)
  
  #Filter out bogus transcripts
  dds <- dds[ rowSums(counts(dds)) > 100, ]
  return(dds)
}
