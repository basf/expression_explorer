output$set1 <- DT::renderDataTable({
    validate(need(!is.null(dataset$meta),message=FALSE))
    dd <- dataset$meta
    dd$Timepoint <- as.character(dd$Timepoint)
    dd = dd[,input$metaShow]
    dd
},server=FALSE,rownames=FALSE,width="100%",
selection = list(mode = 'multiple'),filter = list(position = 'top',clear=TRUE))
  
output$set2 <- DT::renderDataTable({
    validate(need(!is.null(dataset$meta),message=FALSE))
    dd <- dataset$meta
    dd$Timepoint <- as.character(dd$Timepoint)
    dd = dd[,input$metaShow]
    dd
},server=FALSE,rownames=FALSE,width="100%",
selection = list(mode = 'multiple'),filter = list(position = 'top',clear=TRUE))
  
set1Proxy <- dataTableProxy('set1')
set2Proxy <- dataTableProxy('set2')
  
output$set1Selected <- renderText({
    validate(need(length(input$set1_rows_selected)>0,message=FALSE))
    set1 <- rownames(dataset$meta)[input$set1_rows_selected]
    set1 <- paste(set1,collapse=',')
    paste('Selection: ',getReadableSet(set1))
})
  
output$set2Selected <- renderText({
    validate(need(length(input$set2_rows_selected)>0,message=FALSE))
    set2 <- rownames(dataset$meta)[input$set2_rows_selected]
    set2 <- paste(set2,collapse=',')
    paste('Selection: ',getReadableSet(set2))
})
  
output$diffExprPlot <- renderPlot({
    validate(need(length(input$allDiffAnalysis_rows_selected)==1,message=FALSE))
    
    df <- dataset$diffs[[input$allDiffAnalysis_rows_selected]]$diff
    
    if (is.null(input$geneDiff_rows_selected) || length(input$geneDiff_rows_selected)==0) {
      p <- ggplot(data.frame(df),aes(x=log2FoldChange,y=-log10(padj)),color=(padj<0.05)) + geom_point() + theme(legend.position='none')
      p <- p + geom_text_repel(data=subset(df,abs(log2FoldChange)>1),aes(label=Synonym))
      return(p)
    }
    
    genes <- rownames(df)[input$geneDiff_rows_selected]
    #Now get all of their counts
    set1 <- unlist(strsplit(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set1,','))
    set2 <- unlist(strsplit(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set2,','))
    
    if (is.null(dataset$fpkm)) {
      nn <- dataset$norm[genes,c(set1,set2)]
    } else {
      nn <- dataset$fpkm[genes,c(set1,set2)]
    }
    expression <- data.frame(nn,check.names=F)
    expression$Gene <- rownames(expression)
    expression <- melt(expression)
    
    sub <- expression
    sub$syn = dataset$ann[as.character(sub$Gene),]$synonym
    sub$syn[is.na(sub$syn)] = as.character(sub$Gene)[is.na(sub$syn)]
    sub$syn[sub$syn==''] = as.character(sub$Gene)[sub$syn=='']
    #sub <- subset(dataset$expression,Gene %in% genes & variable %in% c(set1,set2))
    
    sub$Set <- 'none'
    sub$Set[sub$variable %in% set1] <- getReadableSet(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set1)
    sub$Set[sub$variable %in% set2] <- getReadableSet(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set2)
    
    if (!input$showHeat_diff) {
      p <- ggplot(sub,aes(x=Set,y=value+1,fill=syn)) + geom_boxplot() + scale_y_log10() + theme(legend.position='bottom')
      if (is.null(dataset$fpkm)) {
        p <- p + ylab('Norm. count') 
      } else {
        p <- p + ylab('FPKM') 
      }
    } else {
      sub = aggregate(value ~ syn + Set, sub, mean)
      p <- ggplot(sub,aes(x=Set,y=syn,fill=value)) + geom_tile() + xlab('Sample') + ylab('Gene')
      p <- p + scale_fill_gradientn( trans = "log", colors = c('white','blue','red'))
    }
    
    p <- p + theme_bw()
    
    return(p)
    
})

output$dwnDiffVenn <- downloadHandler(
    filename = function() { paste('Diff_expression_compasison', '.csv', sep='') },
    content = function(file) {
	#There are multiple selections. Need to join these things.
	df <- allTable()
	cx <- colnames(df)
	for (comp in input$allDiffAnalysis_rows_selected) {
	  sub <- dataset$diffs[[comp]]$diff
	  nn <- paste(getReadableSet(dataset$diffs[[comp]]$set1),' vs ',getReadableSet(dataset$diffs[[comp]]$set2),sep='')
	  
	  df <- cbind(df,sub[rownames(df),]$log2FoldChange)
	  cx <- c(cx,paste(nn,'log2FoldChange'))
	  df <- cbind(df,sub[rownames(df),]$padj)
	  cx <- c(cx,paste(nn,'pval'))
	}
	colnames(df) <- cx
	#Replace commas in "Annotation" since we're writing a csv file:
	df$Annotation <- gsub(',',';',as.character(df$Annotation))
	write.table(df,file,quote=F,sep=',',row.names=F)
    }
)
 
output$dwnDiffBox <- downloadHandler(
    filename = function() { paste('Diff_expression_boxplot', '.pdf', sep='') },
    content = function(file) {
	  validate(need(input$geneDiff_rows_selected,message=FALSE))
	  validate(need(length(input$allDiffAnalysis_rows_selected)==1,message=FALSE))

      df <- dataset$diffs[[input$allDiffAnalysis_rows_selected]]$diff
      genes <- rownames(df)[input$geneDiff_rows_selected]
      #Now get all of their counts
      set1 <- unlist(strsplit(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set1,','))
      set2 <- unlist(strsplit(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set2,','))
      
      nn <- dataset$norm[genes,c(set1,set2)]
      expression <- data.frame(nn,check.names=F)
      expression$Gene <- rownames(expression)
      expression <- melt(expression)
      
      sub <- expression
      #sub <- subset(dataset$expression,Gene %in% genes & variable %in% c(set1,set2))
      
      sub$Set <- 'none'
      sub$Set[sub$variable %in% set1] <- getReadableSet(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set1)
      sub$Set[sub$variable %in% set2] <- getReadableSet(dataset$diffs[[input$allDiffAnalysis_rows_selected]]$set2)
      
      p <- ggplot(sub,aes(x=Set,y=value+1,fill=Gene)) + geom_boxplot() + scale_y_log10() + theme(legend.position='bottom') + theme_bw()
      ggsave(file, plot = p, device = "pdf")
    }
)
  
output$dwnDiffVennSub <- downloadHandler(
    filename = function() { paste('Diff_expression_sub', '.csv', sep='') },
    content = function(file) {
	df <- getAllOverlappingVennGenes()
	
	write.table(df,file,quote=F,sep='\t',row.names=F)
    }
)  

#Update the selection box!
observe({
	  rows_selected = input$allDiffAnalysis_rows_selected
	  
	  for( i in rows_selected) {
		if (i > length(dataset$diffs)) {
			next
		}
        if (!is.null(dataset$diffs[[i]]$diff)) {
				next
			}
			withProgress(message= 'Loading dataset...', value=0.5, {
			dataset$diffs[[i]]$diff <<- as.data.frame(readr::read_csv(dataset$diffs[[i]]$path))

			rownames(dataset$diffs[[i]]$diff) <- dataset$diffs[[i]]$diff$Gene
		})
	  }
	  
	  if (length(rows_selected) > 1) {
		  nms = NULL
		  for (k in input$allDiffAnalysis_rows_selected) {
			  dd = dataset$diffs[[k]]
			  nms = c(nms,paste(getReadableSet(dd$set1),' vs ',getReadableSet(dd$set2),sep=''))
		  }
		  names(rows_selected) = nms
		  updateSelectInput(session,'diff_include',choices=rows_selected,selected=rows_selected)
		  updateSelectInput(session,'diff_exclude',choices=rows_selected,selected=NULL)
	  }
})
 
getAllOverlappingVennGenes <- reactive({
    validate(need(length(input$allDiffAnalysis_rows_selected)>=2,message=FALSE))
    #Cause order is kept in the table, so these should match the list!
    sets_to_include = as.numeric(input$diff_include)
    sets_to_exclude = as.numeric(input$diff_exclude)
    
    keep <- NULL   
    for (k in sets_to_include) {
      dd <- dataset$diffs[[k]]
      relevantSub = subset(dd$diff,	as.numeric(baseMean) > input$baseLimit & 
			abs(as.numeric(log2FoldChange)) > input$log2Limit)
	  if (! input$diff_noSig) {
		  relevantSub = subset(relevantSub,as.numeric(padj)<input$qvalLimit)
	  }
      keep <- c(keep,rownames(relevantSub))
    } 
    
    drop <- NULL
    for (k in sets_to_exclude) {
      dd <- dataset$diffs[[k]]
      dropSub = subset(dd$diff, as.numeric(baseMean) > input$baseLimit 
			& abs(as.numeric(log2FoldChange)) > input$log2Limit)
	  if (! input$diff_noSig) {
		  dropSub = subset(dropSub,as.numeric(padj)<input$qvalLimit)
	  }
      drop <- c(drop,rownames(dropSub))
    }
    
    keep <- names(which(table(keep)==length(sets_to_include)))
    toRem = which(keep %in% drop)
    if (length(toRem) > 0) {
		keep <- keep[-which(keep %in% drop)]
	}
    
    df <- data.frame(row.names=keep)
    comp <- NULL
    for (k in sets_to_include) {
      dd <- dataset$diffs[[k]]
      df <- cbind(df,dd$diff[keep,]$log2FoldChange)
      comp <- c(comp,paste(getReadableSet(dd$set1),' vs ',getReadableSet(dd$set2),sep=''))
    }
    colnames(df) <- comp
    df$Annotation <- gsub(',',';',as.character(dataset$ann[rownames(df),]$annotation))
            
    if (is.null(dataset$ecs) && file.exists(paste(baseFolder,'/',dataset$loaded,'/',"EC_numbers.tab",sep=''))) {
		ecs <- read.table(paste(baseFolder,'/',dataset$loaded,'/',"EC_numbers.tab",sep=''),header=T,row.names=NULL,sep='\t',quote='',comment='',stringsAsFactors=F)
		ecs <- aggregate(EC ~ gene,ecs,function(x) {
					a = unique(ec2AllPath[x,]$V1)
					a = a[!is.na(a)]
					paste(a,collapse=';')
				})
		rownames(ecs) = ecs$gene
		df$Pathway = ecs[rownames(df),]$EC
	} else if (!is.null(dataset$ecs)) {
		df$Pathway = dataset$ecs[rownames(df),]$EC
	}
    
    df
})
 
output$overlappingGene <- DT::renderDataTable({
    getAllOverlappingVennGenes() #* exprFunc 
})

diffData <- observeEvent(input$getDiff, {
    
    validate(need(!is.na(dataset$norm),message=FALSE))
    set1 <- rownames(dataset$meta)[input$set1_rows_selected]
    set2 <- rownames(dataset$meta)[input$set2_rows_selected]

    validate(need(length(set1)>0 & length(set2)>0,"Please select samples to be compared"))

    set1_name <- paste(set1,collapse=',')
    set2_name <- paste(set2,collapse=',')
    DEname <- paste0(set1_name, "_vs_", set2_name)
    
    if(!is.null(dataset$diffs)){
      
      flag <- data.frame(set1 = NA, set2 = NA, s1name = NA, s2name = NA, both = NA)
      
      for(i in 1:length(dataset$diffs)){
        flag[i,1] <- dataset$diffs[[i]]$set1
        flag[i,2] <- dataset$diffs[[i]]$set2
        flag[i,3] <- ifelse(any(flag[i,1:2] %in% set1_name), TRUE, FALSE)
        flag[i,4] <- ifelse(any(flag[i,1:2] %in% set2_name), TRUE, FALSE)
        flag[i,5] <- ifelse(flag[i,3] == TRUE && flag[i,4] == TRUE, TRUE, FALSE)
      }
      flag <- ifelse(any(flag[,5] == TRUE), TRUE, FALSE)
    
    }else{
      flag <- FALSE
    }
    
    if(flag == FALSE){

      dd <- makeDESeqObject(dataset$norm,set1,set2)
      
      validate(need(!is.null(dd),"Both sets should contain at least one valid sample identifier!"))
      
      withProgress(message= 'Computing... this may take a while', value=0.5, {
		  
		#Do we only have 1 sample on each side? Then, forget DESeq, there's no point
	    if (length(set1)==1 & length(set2)==1) {
			dds <- dataset$norm[ rowSums(dataset$norm) > 100, ]
			out <- data.frame(
				row.names = rownames(dds),
				baseMean = rowMeans(dds[,c(set1,set2)]),
				log2FoldChange = log2(dds[,set2]/dds[,set1]),
		 		padj = 1)
		 	out = subset(out,!is.na(out$log2FoldChange))
		 	out$log2FoldChange[is.infinite(out$log2FoldChange) & out$log2FoldChange > 0] = max(out$log2FoldChange[is.finite(out$log2FoldChange)])
		 	out$log2FoldChange[is.infinite(out$log2FoldChange) & out$log2FoldChange < 0] = min(out$log2FoldChange[is.finite(out$log2FoldChange)])
		 	#Convert the Inf and -Inf
		} else {
			res <- DESeq(dd)
			out <- results(res)
		}
        
      })

      sig <- out
      sig <- data.frame(sig)

      showNotification(paste("Found ",nrow(subset(sig,as.numeric(padj)< 0.05))," genes",sep=''))
      sig$Gene <- rownames(sig)
      sig <- sig[,c("Gene","baseMean","log2FoldChange","padj")] 
      sig[,"baseMean"] <- round(sig[,"baseMean"],2)
      sig[,"log2FoldChange"] <- round(sig[,"log2FoldChange"],2)
      #sig[,4] <- formatC(sig[,4],format='e',5)
      sig[,"padj"] <- format.pval(sig[,"padj"], digits=5, eps=1e-16)
      sig[,"padj"] <- as.numeric(gsub("< 1e-16","0",sig[,"padj"]))
      sig$Annotation <- as.character(dataset$ann[rownames(sig),]$annotation)
      sig$Synonym <- as.character(dataset$ann[rownames(sig),]$synonym)
      numDE <- nrow(subset(sig, as.numeric(padj)< 0.05))

      unique_name = uuid::UUIDgenerate()
      old_name = paste0(DEname, "-", numDE ,"_diffExp.csv")
      
      #Keep this result for later
      dataset$diffs[[length(dataset$diffs)+1]] <- list(path=paste0(dataset$folder,"/",unique_name),set1=paste(set1,collapse=','),set2=paste(set2,collapse=','),dgs=numDE,diff=sig)
      dataset$availDiffs[[length(dataset$availDiffs)+1]] <- list(path=paste0(dataset$folder,"/",unique_name),set1=paste(set1,collapse=','),set2=paste(set2,collapse=','),dgs=numDE,diff=sig)

      diffExpr = NULL
      if(file.exists(paste(dataset$folder,'diff_expr.tab',sep='/'))){
        diffExpr = read.table(paste(dataset$folder,'diff_expr.tab',sep='/'),sep='\t')
      }
      diffExpr = rbind(diffExpr,c(unique_name,old_name))
      write.table(diffExpr,paste(dataset$folder,'diff_expr.tab',sep='/'),sep='\t',quote=F)
      readr::write_csv(sig, paste0(dataset$folder,"/",unique_name))
      
    }
    
    #Clear selection
    set1Proxy %>% selectRows(NULL)
    set2Proxy %>% selectRows(NULL)
    
    #Toggle deep dive panel is not toggled yet
    updateCollapse(session,'diffPanel',open=c('doDiff','extra'))
    
})

getReadableSet <- function(set) {
    #Try and collapse names so that it all looks nicer in the table
    set <- unlist(strsplit(set,',')) 
    ss <- dataset$meta[set,]
    samples <- unique(ss$Sample)
    readable = NULL
    for (s in samples){
      ss_sub <- subset(ss,Sample == s)
      time = sort(unique(ss_sub$Timepoint))
      str_time = NULL
      for (t in time) {
	tt_sub <- subset(ss_sub,Timepoint==t)$Replicate
	str_time = c(str_time,paste('T',t,'_','R',paste(tt_sub,collapse=';'),sep=''))
      }
      readable <- c(readable,paste(s,'(',paste(str_time,collapse=';'),')',sep=''))
    }
    return(paste(readable,collapse=' '))
}
  
revertReadableSet <- function(readable){
    #Get sample names
    set <- unlist(strsplit(set1,"\\(T"))
    set <- unlist(strsplit(set,"\\) "))
    samples <- set[grep("_R", set, invert = TRUE)]
    meta <- dataset$meta
    
    used_meta <- list()
    
    for(s in 1:length(samples)){
      sub <- gsub("\\)", "", gsub("^","T",set[which(samples[s] == set) + 1]))
      tt_sub <- unlist(strsplit(sub, "T"))
      tt_sub <- tt_sub[grep("R", tt_sub)]
      length_tt <- length(tt_sub)
      tt <- sapply(strsplit(tt_sub[grep("R", tt_sub)], "_R"), "[[",1)
      
      ss_sub <- gsub("^",";", sapply(strsplit(tt_sub, "_R"), "[[",2))
      
      time_rep <- list()
      
      for(r in 1:length(tt)){
        ss <- unlist(strsplit(ss_sub[r], ";"))
        ss <- ss[2:length(ss)]
        time_rep[[r]] <- data.frame(sample = NA, timepoint = rep(tt[r], length(ss)), replicate = ss, stringsAsFactors = F)
      }
      t_r <- do.call(rbind,time_rep)
      t_r$sample <- samples[s]
      used_meta[[s]] <- t_r
    }
    
    used_meta <- do.call(rbind,used_meta)
    
    keep <- meta[which((meta$Sample %in% used_meta$sample) & (meta$Timepoint == used_meta$timepoint) & (meta$Replicate == used_meta$replicate)),]
    
    return(paste(rownames(keep),collapse = ","))
}
 
output$allDiffAnalysis <- DT::renderDataTable({
    #req(!is.null(input$availDiffs))
     validate(need(!is.null(dataset$availDiffs),message=FALSE))
    df <- NULL
    for (i in 1:length(dataset$availDiffs)) {
      if(is.na(dataset$availDiffs[[i]]$set1)){
        next
      }
      df <- rbind(df,c(getReadableSet(dataset$availDiffs[[i]]$set1),getReadableSet(dataset$availDiffs[[i]]$set2),dataset$availDiffs[[i]]$dgs)) 
    }
    
    df <- data.frame(df,stringsAsFactors=F)
    colnames(df) <- c('Set1','Set2','Number of diff. genes')
    df
},rownames=FALSE,server=TRUE)

diffGeneTableData <- reactive({
	req(length(input$allDiffAnalysis_rows_selected) == 1, !is.null(input$allDiffAnalysis_rows_selected))
    showDG <- dataset$diffs[[input$allDiffAnalysis_rows_selected]]$diff

    isolate({
    	if (is.null(dataset$ecs) && file.exists(paste(baseFolder,'/',dataset$loaded,'/',"EC_numbers.tab",sep=''))) {
			withProgress(message= 'Loading and formatting EC numbers...', value=0.5, {
				ecs <- read.table(paste(baseFolder,'/',dataset$loaded,'/',"EC_numbers.tab",sep=''),header=T,row.names=NULL,sep='\t',quote='',comment='')
				ecs <- aggregate(EC ~ gene,ecs,function(x) {
							sbs = subset(ec2Path,V2 %in% x)
							paste(path2Name[unique(sbs$V1),]$V2,collapse=';')
						})
				rownames(ecs) = ecs$gene
				dataset$ecs = ecs
				showDG$Pathway = ecs[rownames(showDG),]$EC
			})
    	} else if (!is.null(dataset$ecs)) {
			showDG$Pathway = dataset$ecs[rownames(showDG),]$EC
		}
    })
    showDG
})
 
output$geneDiff <- DT::renderDataTable({
    
    req(length(input$allDiffAnalysis_rows_selected) == 1, !is.null(input$allDiffAnalysis_rows_selected))
    diffGeneTableData()
    
},rownames=FALSE,server=TRUE,filter = list(position = 'top',clear=TRUE))
  
geneDiffProxy <- dataTableProxy('geneDiff')
allDiffAnalysisProxy <- dataTableProxy('allDiffAnalysis')

output$diffHierarchicalClustering <- renderPlot({
  validate(need(length(input$allDiffAnalysis_rows_selected)==1,message=FALSE))
  validate(need(length(input$geneDiff_rows_selected>1),message=FALSE))
  
  df <- dataset$diffs[[input$allDiffAnalysis_rows_selected]]$diff
  genes <- rownames(df)[input$geneDiff_rows_selected]
  
  if (length(genes) < 2) {
    return()
  }
  
  #Get correlation of these genes across everything
  all <- dataset$fpkm
  all <- all[genes,]
  all <- sqrt(all)
  rownames(all) = dataset$ann[rownames(all),]$synonym
  c <- cor(t(all),method='pearson')
  h = hclust(as.dist(1-c)) #This might be a stupid distance, but ... meh
  plot(h) #This can be 1000 times prettier, but i don't have time for it now.
})

output$dwnPlot <- downloadHandler(
    filename = function() { paste('Expression_plot', '.pdf', sep='') },
    content = function(file) {
      df <- dataToPlot()
      p <- ggplot(df,aes(x=timepoint,y=value,group=Gene,color=syn,alpha=highlight)) + geom_point() + geom_line() + scale_alpha_manual(values=c(1,0.3),guide = FALSE)
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10), legend.position='bottom') + xlab('Timepoint') + scale_y_log10()
      p <- p + geom_errorbar(aes(ymin=value-(se/2),ymax=value+(se/2)),width=.1)
      p <- p + facet_grid(~ sample) + theme_bw()
      if (is.null(dataset$fpkm)) {
		p <- p + ylab('Norm. count') 
	  } else {
		p <- p + ylab('FPKM') 
	  }
      ggsave(file, plot = p, device = "pdf")
    }
)
 
output$dwnRawGene <- downloadHandler(
    filename = function() { paste('Expression_data', '.csv', sep='') },
    content = function(file) {
      df <- dataToPlot()
      df = df[,1:4]
      if (!is.null(dataset$fpkm)) {
	      colnames(df)[4] = 'Fpkm'
	  } else {
		  colnames(df)[4] = 'Counts'
	  }
      write.table(df,file,quote=F,sep=',',row.names=F)
    }
)
 
output$dwnDiff <- downloadHandler(
    filename = function() { paste('Diff_expression', '.csv', sep='') },
    content = function(file) {
		df <- diffGeneTableData()
		df$Annotation <- gsub(',',';',as.character(df$Annotation))
		write.table(df,file,quote=F,sep=',',row.names=F)
      }
)
 
output$dwnDiffFilt <- downloadHandler(
    filename = function() { paste('Diff_expression_filtered', '.csv', sep='') },
    content = function(file) {
		df <- diffGeneTableData()
		df <- df[input$geneDiff_rows_all,]
		df$Annotation <- gsub(',',';',as.character(df$Annotation))
		write.table(df,file,quote=F,sep=',',row.names=F)
    }
)
observeEvent(input$toggleFilter1, {
    set1Proxy %>% selectRows(c(input$set1_rows_selected,input$set1_rows_all))
})
  
observeEvent(input$clearFilter1, {
    set1Proxy %>% selectRows(NULL)
})
  
observeEvent(input$toggleFilter2, {
    set2Proxy %>% selectRows(c(input$set2_rows_selected,input$set2_rows_all))
})
  
observeEvent(input$clearFilter2, {
    set2Proxy %>% selectRows(NULL)
})
  
observeEvent(input$clearDiff, {
    geneDiffProxy %>% selectRows(NULL)
})
  
observeEvent(input$deleteDiffAnalysis, {

    allDiffAnalysisProxy %>% selectRows(NULL)
    row_to_del <- as.numeric(gsub("Row","",input$allDiffAnalysis_rows_selected))
    rm_res <- file.remove(dataset$diffs[[row_to_del]]$path)

    if (rm_res) {
      diffExpr = read.table(paste(dataset$folder,'diff_expr.tab',sep='/'),sep='\t')
      pp = strsplit(dataset$diffs[[row_to_del]]$path,'/')[[1]]
      pp = pp[length(pp)]
      toRem = grep(pp,diffExpr[,1])
      if (length(toRem) == 1) {
        diffExpr = diffExpr[-toRem,,drop=F]
        if (nrow(diffExpr) == 0) {
          file.remove(paste(dataset$folder,'diff_expr.tab',sep='/'))
        } else {
          write.table(diffExpr,paste(dataset$folder,'diff_expr.tab',sep='/'),sep='\t',quote=F)
        }
      }
      
      dataset$availDiffs <- dataset$availDiffs[-row_to_del]
      dataset$diffs <- dataset$diffs[-row_to_del]
      if(length(dataset$availDiffs) == 0){
        dataset$availDiffs <- NULL
        dataset$diffs <- NULL
      }
      
    } else {
      # show some rror? e.g. shinytoastr or shinyBS::altert
      showModal(modalDialog(
        title = "This file doesn't exist",
        footer = tagList(
          modalButton("Close")
        )
      ))
    }
})

output$vennPlot <- renderPlot({
    validate(need(length(input$allDiffAnalysis_rows_selected)>=2,message=FALSE))
    validate(need(length(input$allDiffAnalysis_rows_selected)<6,message='No way! Are you crazy...'))
    #Cause order is kept in the table, so these should match the list!
    all <- list()
    names <- NULL
    for (k in input$allDiffAnalysis_rows_selected) {
      dd <- dataset$diffs[[k]]
      relSub = subset(dd$diff, as.numeric(baseMean) > input$baseLimit & abs(as.numeric(log2FoldChange)) > input$log2Limit)
      if (!input$diff_noSig) {
		  relSub = subset(relSub,as.numeric(padj)<input$qvalLimit)
	  }
      all[[length(all)+1]] <- rownames(relSub)
      names <- c(names,paste(getReadableSet(dd$set1),getReadableSet(dd$set2),sep=' vs '))
    } #* renderPlot
    names(all) <- names
    #grid.draw(myvenn(all,filename=NULL,imagetype='png',fill=rainbow(length(input$allDiffAnalysis_rows_selected))))
    v <- myvenn(all,filename=NULL,imagetype='png',fill=rainbow(length(names)),category.names=rep('',length(names)))
    cols <- rainbow(length(names))
    lg <- legendGrob(labels=names(all), pch=rep(19,length(names(all))),
                 gp=gpar(col=cols, fill="gray"),
                 byrow=TRUE)
    grid.arrange(gTree(children=gList(v)),lg,nrow=2,heights=c(5,1))
  })
 
  output$dwnDiffVennPlot <- downloadHandler(
    filename = function() { paste('Diff_expression_venn', '.pdf', sep='') },
    content = function(file) {
      #Cause order is kept in the table, so these should match the list!
      all <- list()
      names <- NULL

      for (k in input$allDiffAnalysis_rows_selected) {
	dd <- dataset$diffs[[k]]
	relSub = subset(dd$diff, as.numeric(baseMean) > input$baseLimit & abs(as.numeric(log2FoldChange)) > input$log2Limit)
      if (!input$diff_noSig) {
		  relSub = subset(relSub,as.numeric(padj)<input$qvalLimit)
	  }
      all[[length(all)+1]] <- rownames(relSub)
	names <- c(names,paste(getReadableSet(dd$set1),getReadableSet(dd$set2),sep=' vs '))
      }
      names(all) <- names
      v <- myvenn(all,filename=NULL,imagetype='png',fill=rainbow(length(names)),category.names=rep('',length(names)))
      cols <- rainbow(length(names))
      lg <- legendGrob(labels=names(all), pch=rep(19,length(names(all))),
		  gp=gpar(col=cols, fill="gray"),
		  byrow=TRUE)
      p <- grid.arrange(gTree(children=gList(v)),lg,nrow=2,heights=c(5,1))
      ggsave(file, plot = p, device = "pdf")
    }
)

observeEvent(input$enrichA, {
    validate(need(length(input$allDiffAnalysis_rows_selected)==1,message=FALSE))
    #Check for categorical annotations for this analysis
    f <- list.files(paste(baseFolder,dataset$loaded,sep='/'),'*_enrich.tab')
    if (length(f) == 0) {
	showModal(modalDialog(
		  title = "Enrichment analysis",
		  "There is no annotation available for this dataset",
		  footer = tagList(
		    modalButton("Close")
		  )
	))
    } else {
	n <- gsub("_enrich.tab","",f)
        showModal(modalDialog(
		  title = "Enrichment analysis",
		  selectInput("enrich_pick","Please select the annotation to perform enrichment analysis with: ",n),
		  footer = tagList(
		    actionButton("enrich_run","Compute Enrichment"),
		    modalButton("Close")
		  )
	))
    }
})

observeEvent(input$enrichB, {
    validate(need(length(input$allDiffAnalysis_rows_selected)>1,message=FALSE))
    #Check for categorical annotations for this analysis
    f <- list.files(paste(baseFolder,dataset$loaded,sep='/'),'*_enrich.tab')
    if (length(f) == 0) {
	showModal(modalDialog(
		  title = "Enrichment analysis",
		  "There is no annotation available for this dataset",
		  footer = tagList(
		    modalButton("Close")
		  )
	))
    } else {
	n <- gsub("_enrich.tab","",f)
        showModal(modalDialog(
		  title = "Enrichment analysis",
		  selectInput("enrich_pick","Please select the annotation to perform enrichment analysis with: ",n),
		  footer = tagList(
		    actionButton("enrich_runB","Compute Enrichment"),
		    modalButton("Close")
		  )
	))
    }
})
 
enrich <- reactiveValues(computed=NULL,annotations=NULL)

observeEvent(input$enrich_run, {
    toUse <- input$enrich_pick
    #read the input for this
    withProgress(message= 'Reading enrichment annotation...', value=0.2, {
	enrichData <- read.table(paste(baseFolder,'/',dataset$loaded,'/',toUse,"_enrich.tab",sep=''),header=T,row.names=1,sep='\t',quote='',comment='')
    })
    
    ####Parse this weird format. For now, this is a silly custom thing.
    withProgress(message= 'Parsing annotation...', value=0.4, {
      ann <- unlist(strsplit(as.character(enrichData[["annotation"]]),','))
      desc <- gsub('\"','',unlist(strsplit(as.character(enrichData[["description"]]),',')))
      dd <- data.frame(ann=ann,desc=desc,stringsAsFactors=F)
      dd <- subset(dd,!is.na(ann))
      #Here we have the background!
      cts <- aggregate(desc~ann,dd,length)
      rownames(cts) <- cts$ann
      dd <- dd[!duplicated(dd$ann),]
      rownames(dd) <- dd$ann
    })
    ####Get the selected list
    full <- dataset$diffs[[input$allDiffAnalysis_rows_selected]]$diff
    #Get the filtered stuff.
    full <- full[input$geneDiff_rows_all,]
    
    full$ann <- as.character(enrichData[full$Gene,]$annotation)
    enrich$annotations = full
    myann <- unlist(strsplit(full$ann,','))
    myann <- myann[!is.na(myann)]
    #These are all my categories
    mytable <- table(myann)
    
    withProgress(message= 'Computing enrichment...', value=0.8, {
      res <- NULL
      for (t in names(mytable)) {
	cont_table = matrix( c(mytable[t], sum(mytable), cts[t,]$desc, sum(cts$desc)), 2, 2)
	p <- fisher.test(cont_table,alternative='greater')
	res <- rbind(res,c(t,mytable[t],cts[t,]$desc,p$p.value,p$estimate))
      }
    })
    
    res <- data.frame(res,stringsAsFactors=F)
    res$adj <- p.adjust(as.numeric(res[,4]))
    colnames(res) <- c('Term','Count','Background_count','pval','odds','padj')
    res$Description <- dd[res$Term,]$desc
    enrich$computed <- res[,c(1,7,6,5,2,3,4)]
    removeModal()
    updateNavbarPage(session,'nPage',selected='Enrichment Analysis')
})
 
observeEvent(input$enrich_runB, {
    toUse <- input$enrich_pick
    #read the input for this
    withProgress(message= 'Reading enrichment annotation...', value=0.2, {
	enrichData <- read.table(paste(baseFolder,'/',dataset$loaded,'/',toUse,"_enrich.tab",sep=''),header=T,row.names=1,sep='\t',quote='',comment='')
    })
    
    ####Parse this weird format. For now, this is a silly custom thing.
    withProgress(message= 'Parsing annotation...', value=0.4, {
      ann <- unlist(strsplit(as.character(enrichData[["annotation"]]),' '))
      desc <- gsub('\"','',unlist(strsplit(as.character(enrichData[["description"]]),';')))
      dd <- data.frame(ann=ann,desc=desc,stringsAsFactors=F)
      dd <- subset(dd,!is.na(ann))
      #Here we have the background!
      cts <- aggregate(desc~ann,dd,length)
      rownames(cts) <- cts$ann
      dd <- dd[!duplicated(dd$ann),]
      rownames(dd) <- dd$ann
    })
    ####Get the selected list
    full <- getAllOverlappingVennGenes()
    
    full$ann <- as.character(enrichData[rownames(full),]$annotation)
    enrich$annotations <- full
    myann <- unlist(strsplit(full$ann,' '))
    myann <- myann[!is.na(myann)]
    #These are all my categories
    mytable <- table(myann)
    
    withProgress(message= 'Computing enrichment...', value=0.8, {
      res <- NULL
      for (t in names(mytable)) {
	cont_table = matrix( c(mytable[t], sum(mytable), cts[t,]$desc, sum(cts$desc)), 2, 2)
	p <- fisher.test(cont_table,alternative='greater')
	res <- rbind(res,c(t,mytable[t],cts[t,]$desc,p$p.value,p$estimate))
      }
    })
    
    res <- data.frame(res,stringsAsFactors=F)
    res$adj <- p.adjust(as.numeric(res[,4]))
    colnames(res) <- c('Term','Count','Background_count','pval','odds','padj')
    res$Description <- dd[res$Term,]$desc
    enrich$computed <- res[,c(1,7,6,5,2,3,4)]
    removeModal()
    updateNavbarPage(session,'nPage',selected='Enrichment Analysis')
})

output$overlapDiffPlot <- renderPlot({
    validate(need(length(input$allDiffAnalysis_rows_selected)>=2,message=FALSE))
    #Cause order is kept in the table, so these should match the list!
    keep <- NULL

    for (k in input$allDiffAnalysis_rows_selected) {
      dd <- dataset$diffs[[k]]
      relSub = subset(dd$diff, as.numeric(baseMean) > input$baseLimit & abs(as.numeric(log2FoldChange)) > input$log2Limit)
      if (!input$diff_noSig) {
		relSub = subset(relSub,as.numeric(padj)<input$qvalLimit)
	  }
      keep <- c(keep,rownames(relSub))
    }
     #*  renderPlot 
    keep <- names(which(table(keep)==length(input$allDiffAnalysis_rows_selected)))
    df <- data.frame(row.names=keep)
    comp <- NULL
    for (k in input$allDiffAnalysis_rows_selected) {
      dd <- dataset$diffs[[k]]
      df <- cbind(df,dd$diff[as.character(keep),]$log2FoldChange)
      comp <- c(comp,paste(getReadableSet(dd$set1),' vs ',getReadableSet(dd$set2),sep=''))
    }
    colnames(df) <- comp
    plot(df)
  })
 
output$allOverlapPlot <- renderPlot({
    validate(need(length(input$allDiffAnalysis_rows_selected)>=2,message=FALSE))
    validate(need(length(input$overlappingGene_rows_selected)>0,message=FALSE))
    #Cause order is kept in the table, so these should match the list!
    allSamples <- NULL

    for (k in input$allDiffAnalysis_rows_selected) {
      dd <- dataset$diffs[[k]]
      smpl1 <- unlist(strsplit(dd$set1,','))
      smpl2 <- unlist(strsplit(dd$set2,','))
      d1 <- data.frame(sample=smpl1,set=getReadableSet(smpl1))
      d2 <- data.frame(sample=smpl2,set=getReadableSet(smpl2))
      allSamples <- rbind(allSamples,d1)
      allSamples <- rbind(allSamples,d2)
    }
    
    #rownames(allSamples) <- as.character(allSamples$sample)
    over <- getAllOverlappingVennGenes()[input$overlappingGene_rows_selected,,drop=F]
    nn <- dataset$norm[rownames(over),unique(as.character(allSamples$sample))]
    nn$Gene <- rownames(nn)
    nn <- melt(nn)
    
    nn$set <- allSamples[match(as.character(nn$variable),as.character(allSamples$sample)),]$set
    
    ggplot(nn,aes(x=set,y=value,fill=Gene)) + geom_boxplot() + scale_y_log10()
})

observeEvent(input$pathwayS, {
    validate(need(length(input$allDiffAnalysis_rows_selected)>0,message=FALSE))
    updateNavbarPage(session,'nPage',selected='Pathway View')
})
