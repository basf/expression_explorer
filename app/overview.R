############ Samples table

output$samples <- DT::renderDataTable({	  
    validate(need(!is.null(dataset$meta),message=FALSE))
    #dd <- dataset$meta[,c('Sample','Timepoint','Replicate')]
    dd <- dataset$meta
    dd$Timepoint <- as.character(dd$Timepoint)
    dd = dd[,input$metaShow] #Subset this!
    datatable(dd,rownames=input$showRawNames,filter = list(position = 'top',clear=TRUE))
})

sampleDataSets <- reactive({
  dd <- dataset$sampleSets
  dd <- aggregate(Sample ~ Set, dd, length)
  colnames(dd) <- c('Set','Sample count')
  return(dd)
})

output$sample_sets <- DT::renderDataTable({	  
  validate(need(!is.null(dataset$sampleSets),message='No sample sets defined. Go to the Metadata tab to define sample sets'))
  dd <- sampleDataSets()
  datatable(dd,rownames=FALSE)
})

samplesProxy <- dataTableProxy('samples')

observeEvent(input$clearSample, {
    samplesProxy %>% selectRows(NULL)
})
  
observeEvent(input$selectAllSample, {
	if (nrow(dataset$meta) > 150) {
		showModal(modalDialog(
		  title = "Too much data",
		  "You are trying to load too many samples at the same time. This will make everything super slow. Please subselect the relevant samples by using the 'Select Filtered' option",
		  footer = tagList(
			modalButton("Close")
		  )
	))
	} else {
	  samplesProxy %>% selectRows(1:nrow(dataset$meta))
	}
})

observeEvent(input$sample_sets_rows_selected, ignoreNULL = FALSE, {
  #Get the samples in all of these selected sets.
  if (length(input$sample_sets_rows_selected) == 0) {
    samplesProxy %>% selectRows(NULL)
    return()
  }
  sets = sampleDataSets()[input$sample_sets_rows_selected,]$Set
  samples = subset(dataset$sampleSets, Set %in% sets)$Sample
  selection = unique(match(samples,rownames(dataset$meta)))
  samplesProxy %>% selectRows(selection)
})

observeEvent(input$selectFilteredSample, {
	samplesProxy %>% selectRows(c(input$samples_rows_selected,input$samples_rows_all))
})

############ PCA plot

pca <- reactive({
    validate(need(input$samples_rows_selected,message='Select the samples which you want to see'))
    
    final <- log10(dataset$norm[,rownames(dataset$meta)[input$samples_rows_selected]]+1)
    validate(need(ncol(final)>0,message=FALSE))
    
    pp <- prcomp(t(final))
    
    eig <- pp[["sdev"]]
    eig <- eig/sum(eig)*100
    
    sizes <- pp[["center"]]
    rot <- data.frame(pp[["rotation"]][,1:2])
    rot <- rot*200
    
    toUse <- sort(sizes,decreasing=T)[1:input$loads]
    rot <- rot[names(toUse),]
    
    pos <- data.frame(pp[["x"]][,1:2])
    pos <- cbind(pos,dataset$meta[rownames(pos),])
        
    list(pos,rot,eig)
})
  
output$plot3 <- renderPlot({
  
    res <- pca()
    pos <- res[[1]]
    rot <- res[[2]]
    eig <- res[[3]]
    
    pos_name = paste(pos$Sample,pos$Timepoint,pos$Replicate,sep='_')
    pos[[input$pca_color]] = as.factor(pos[[input$pca_color]])
    
    p <- ggplot() + geom_point(data=pos,aes_string(x="PC1",y="PC2",color=input$pca_color),size=3) + geom_text_repel(data=pos,aes(x=PC1,y=PC2),label=pos_name) 
    p <- p + theme_bw() + geom_segment(data=rot,x=0,y=0,xend=rot[,1],yend=rot[,2],arrow=arrow(length = unit(0.15, "cm")),alpha=0.5,color='red')
    p <- p + geom_text(data=rot,aes(x=PC1,y=PC2),label=rownames(rot),nudge_y=0.1,nudge_x=-0.2,col='red',alpha=0.5)
    p <- p + xlab(sprintf("%3.2f%%",eig[1])) + ylab(sprintf("%3.2f%%",eig[2]))
    
    if (input$show_pca_group) {
      p <- p + stat_ellipse(data=pos,aes_string(color=input$pca_color,x="PC1",y="PC2"))
    }
    
    p
    
})

output$dwnAllGene <- downloadHandler(
  filename = function() { 'Count_data.tab' },
  content = function(file) {
    write.table(dataset$fpkm,file,quote=F,sep='\t')
  }
)

output$dwnPCA <- downloadHandler(
    filename = function() { paste('PCA_plot', '.pdf', sep='') },
    content = function(file) {
      res <- pca()
      pos <- res[[1]]
      rot <- res[[2]]
      eig <- res[[3]]
      
      pos_name = paste(pos$Sample,pos$Timepoint,pos$Replicate,sep='_')
      
      p <- ggplot() + geom_point(data=pos,aes(x=PC1,y=PC2,color=as.factor(Sample),shape=as.factor(Timepoint)),size=3) + geom_text_repel(data=pos,aes(x=PC1,y=PC2),label=pos_name) 
      p <- p + theme_bw() + geom_segment(data=rot,x=0,y=0,xend=rot[,1],yend=rot[,2],arrow=arrow(length = unit(0.15, "cm")),alpha=0.5,color='red')
      p <- p + geom_text(data=rot,aes(x=PC1,y=PC2),label=rownames(rot),nudge_y=0.1,nudge_x=-0.2,col='red',alpha=0.5)
      p <- p + xlab(sprintf("%3.2f%%",eig[1])) + ylab(sprintf("%3.2f%%",eig[2]))
      p
      ggsave(file, plot = p, device = "pdf")
    }
)

############ Gene table and relevant buttons

baseTable <- reactive({
	dd = dataset$norm
	alltab <- data.frame(row.names=rownames(dd),Gene=rownames(dd),stringsAsFactors=F)
	alltab$Ann <- as.character(dataset$ann[as.character(alltab$Gene),]$annotation)
	alltab$Syn <- as.character(dataset$ann[as.character(alltab$Gene),]$synonym)
	colnames(alltab) <- c('Gene','Annotation','Synonym')
	alltab
})

allTable <- reactive({
    
    if (input$updateOnSelect && length(input$samples_rows_selected) > 0) {
		if (!is.null(dataset$fpkm)) {
			dd = dataset$fpkm[,rownames(dataset$meta)[input$samples_rows_selected],drop=F]	
		} else {
			dd = dataset$norm[,rownames(dataset$meta)[input$samples_rows_selected],drop=F]
		}
		counts = rowSums(dd)/ncol(dd)
		alltab = baseTable()
		alltab$Mean <- counts[rownames(alltab)]
		colnames(alltab) <- c('Gene','Annotation','Synonym','Mean_value')
		return (alltab)
    
	}
	return(baseTable())
})

output$genesInRange <- DT::renderDataTable({
    t <- allTable()
    validate(need(!is.null(t),message='You should load a dataset first!'))
    tt = datatable(t, rownames=FALSE,
    editable=TRUE,
    extensions = 'Buttons',
    filter = list(position = 'top',clear=FALSE),
    options = list(
      search = list(regex = TRUE, caseInsensitive = TRUE),
      pageLength = 10,
      dom = 'lBrtip',
      buttons = c('excel', 'pdf', 'print')))
    if ("Mean_value" %in% colnames(t)) {
      tt = tt %>% formatRound('Mean_value', 2)
   }
   tt
  },server=TRUE )

proxy = dataTableProxy('genesInRange')
  
observeEvent(input$genesInRange_cell_edit, {
	  info = input$genesInRange_cell_edit
	  
	  if (info$col == 0) {
		  return(NULL)
	  }	  
	  
	  i = info$row
	  v = info$value
	  from = allTable()
	  name = rownames(from)[i]
	  dataset$ann[name,info$col] = v
	  isolate({
        write.table(dataset$ann,paste(baseFolder,dataset$loaded,'ann.tab',sep='/'),sep='\t',quote=F)
	  })
})
  
observeEvent(input$clear1, {
    proxy %>% selectRows(NULL)
})

observeEvent(input$pre_def_gene_list, {
  sel = input$pre_def_gene_list
  if (sel != "") {
    genes = subset(dataset$geneSets,Set==sel)$Gene
    updateTextAreaInput(session,'gene',value=paste(genes,collapse=' '))
  }
})

observeEvent(input$saveList, {
  setName = input$set_name
  if (setName == "") {
    showModal(modalDialog(
      title = "Info missing",
      "Please provide a name for this gene set.",
      footer = tagList(
        modalButton("Close")
      )
    ))
  } else {
    
    #Get all the genes
    gg <- gsub('\n',' ',input$gene)
    gg <- gsub('\t',' ',gg)
    gg <- unlist(strsplit(gg,' '))
    if (length(gg)==0) {#There's nothing to save here
      return()
    }
    new = dataset$geneSets
    if (!is.null(new)) {
      new = subset(dataset$geneSets,Set != setName)#We're overwriting whatever was here before
      new = rbind(new,data.frame(Gene=gg,Set=setName))
    } else {
      new = data.frame(Gene=gg,Set=setName)
    }
    #Now write this
    write.table(new,paste(dataset$folder,'genesets.tab',sep='/'),quote=F,sep='\t',row.names=F)
    #Update our object
    dataset$geneSets = new
    #Update the dropdown
    updateSelectInput(session,'pre_def_gene_list',choices=c("",unique(dataset$geneSets$Set)),selected="")
    updateTextInput(session,'set_name',value="")
  }
})

############# Gene plots

selectedData <- reactive({
    validate(need(length(input$samples_rows_selected)>0,message='Please select the sample set first'))
    all <- allTable()
    validate(need(!is.null(all),message=FALSE))
    gg <- gsub('\n',' ',input$gene)
    gg <- gsub('\t',' ',gg)
    gg <- unlist(strsplit(gg,' '))   
    gg <- c(gg,all[input$genesInRange_rows_selected,]$Gene)
    
    validate(need(length(gg)>0,message=FALSE))
    
    if (!is.null(dataset$fpkm)) {
		nn <- dataset$fpkm[gg,rownames(dataset$meta)[input$samples_rows_selected],drop=F]
	} else {
		nn <- dataset$norm[gg,rownames(dataset$meta)[input$samples_rows_selected],drop=F]
	}
    span <- NULL
    if (!is.null(dataset$span)) {
      span = dataset$span[gg,rownames(dataset$meta)[input$samples_rows_selected],drop=F]
      span$Gene <- rownames(span)
      span = melt(span)
    }
    
    mut <- NULL
    if (!is.null(dataset$mutations)) {
      mut = dataset$mutations[gg,rownames(dataset$meta)[input$samples_rows_selected],drop=F]
      mut$Gene = gg
      mut = melt(mut,id.vars='Gene')
      mut$value[is.na(mut$value)] = 'ref' #Don't specifically "name" the reference
    }
    
    expression <- data.frame(nn,check.names=F)
    expression$Gene <- rownames(expression)
    expression <- melt(expression)
    
    sub <- expression
    if (nrow(sub) < 1) {
      return(NULL)
    }
    sub$sample <- dataset$meta[as.character(sub$variable),]$Sample
    sub$sample <- factor(sub$sample,levels=unique(dataset$meta$Sample),ordered=T)
    if (input$splitRepl) {
      sub$sample <- paste(dataset$meta[as.character(sub$variable),]$Sample,dataset$meta[as.character(sub$variable),]$Replicate,sep='_')
      sub$sample <- factor(sub$sample,levels=unique(sub$sample),ordered=T)
    }
    
    sub$timepoint <- as.factor(dataset$meta[as.character(sub$variable),]$Timepoint)
    sub$replicate <- as.factor(dataset$meta[as.character(sub$variable),]$Replicate)
      
    se <- aggregate(value ~ sample + timepoint + Gene,sub,sd)
    mmin <- aggregate(value ~ sample + timepoint + Gene,sub,min)
    mmax <- aggregate(value ~ sample + timepoint + Gene,sub,max)
    m <- aggregate(value ~ sample + timepoint + Gene,sub,mean)
    if (!is.null(span)) {
      sub$span <- span$value
      m_span <- aggregate(span ~ sample + timepoint + Gene,sub,min)#Worst case span
      m$frag <- m_span$span < 0.2 #This is where the threshold is
    }
    
    m$se <- se$value
    m$mmin <- mmin$value
    m$mmax <- mmax$value

    m$syn = dataset$ann[as.character(m$Gene),]$synonym
    m$syn[is.na(m$syn)] = as.character(m$Gene)[is.na(m$syn)]
    m$syn[m$syn==''] = as.character(m$Gene)[m$syn=='']
    
    #Here we mess with the synonym if we have variant information
    if (!is.null(mut)) {
      sub$variant <- mut$value
      m_mut <- aggregate(variant ~ sample + timepoint + Gene,sub,function(x) {paste(sort(unique(x)),collapse='_')})
      is_var <- which(m_mut$variant != 'ref')
      if (length(is_var) > 0) {
        m$syn[is_var] = paste(m$syn[is_var],m_mut$variant[is_var],sep='_')
      }
    }
    
    m$type = 'Rna'
    m
    
    
    if (!is.null(dataset$prot) && !input$hideProt) {### Add protein data to this!
      toSelect = rownames(dataset$meta)[input$samples_rows_selected]
      toSelect = toSelect[which(toSelect %in% colnames(dataset$prot))]
      if (length(toSelect) == 0) {
        return(m)
      }
      prot <- dataset$prot[gg[which(gg %in% rownames(dataset$prot))],toSelect,drop=F]
      if (nrow(prot) == 0) {
        return(m)
      }
      expression <- data.frame(prot,check.names=F)
      expression$Gene <- rownames(expression)
      expression <- melt(expression)
      
      sub <- expression
      if (nrow(sub) < 1) {
        return(m)
      }
      sub$sample <- dataset$meta[as.character(sub$variable),]$Sample
      sub$sample <- factor(sub$sample,levels=unique(dataset$meta$Sample),ordered=T)
      
      if (input$splitRepl) {
        sub$sample <- paste(dataset$meta[as.character(sub$variable),]$Sample,dataset$meta[as.character(sub$variable),]$Replicate,sep='_')
        sub$sample <- factor(sub$sample,levels=unique(sub$sample),ordered=T)
      }
      
      sub$timepoint <- as.factor(dataset$meta[as.character(sub$variable),]$Timepoint)
      sub$replicate <- as.factor(dataset$meta[as.character(sub$variable),]$Replicate)
      
      se <- aggregate(value ~ sample + timepoint + Gene,sub,sd)
      mmin <- aggregate(value ~ sample + timepoint + Gene,sub,min)
      mmax <- aggregate(value ~ sample + timepoint + Gene,sub,max)
      mp <- aggregate(value ~ sample + timepoint + Gene,sub,mean)
      if (!is.null(span)) {#Just faking this in there
        mp$frag = FALSE
      }
      
      mp$se <- se$value
      mp$mmin <- mmin$value
      mp$mmax <- mmax$value
      
      mp$syn = dataset$ann[as.character(mp$Gene),]$synonym
      mp$syn[is.na(mp$syn)] = as.character(mp$Gene)[is.na(mp$syn)]
      mp$syn[mp$syn==''] = as.character(mp$Gene)[m$syn=='']
      
      mp$type = 'Protein'
      m = rbind(m,mp)
    }
    
    return(m)
    
})

highlight <- reactiveValues(metaUpdate=FALSE,
                    networkUsed=NA)
    
dataToPlot <- reactive({
    df <- selectedData()
    validate(need(!is.null(df),message=FALSE))
    
    if (nrow(df) != length(highlight$clicked)) {#Reset! This is probalby a new selection
      highlight$clicked <- TRUE
    }
    df$highlight <- highlight$clicked
    
    df
})

output$plot1 <- renderPlot({
    df <- dataToPlot()
    validate(need(nrow(df)>0,"No genes selected..."))
    
    ####There is a problem with the y-axis ticks of the values are within one log order (which is silly)
    r <- range(log10(c(df$value+1,df$mmin+1,df$mmax+1)))
    #Get closest integer to each
    r[1] <- floor(r[1])
    r[2] <- ceiling(r[2])
    lim <- 10^c(r[1],r[2])
    
    if(!input$showHeat) {#Ok, we're doing normal plotting!
        if (!input$showFlip) {
          ####
          if (input$showTime) {
            p <- ggplot(df,aes(x=timepoint,y=value+1,group=Gene,color=syn,alpha=highlight)) + geom_line() + scale_alpha_manual(values=c(1,0.3),guide = 'none')
          } else {
            p <- ggplot(df,aes(x=as.numeric(as.character(timepoint)),y=value+1,group=Gene,color=syn,alpha=highlight)) + geom_line() + scale_alpha_manual(values=c(1,0.3),guide = 'none')
          }
          
          if ("frag" %in% colnames(df)) {
            p <- p + geom_point(aes(shape=frag), size=3) + scale_shape_manual(values = c(16,1))
          } else {
            p <- p + geom_point(size=3)
          }
          
          p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10), legend.position='bottom') + xlab('Timepoint')
          p <- p + geom_errorbar(aes(ymin=mmin,ymax=mmax),width=.1)
          
          if (input$showLog) {
      		p <- p + scale_y_log10()
      	}
          if (is.null(dataset$fpkm)) {
      		p <- p + ylab('Norm. count') 
      	} else {
      		p <- p + ylab('FPKM') 
      	}
      	
          p = p + facet_grid(type ~ sample,scales='free_y') + theme_bw()
        }else {
          if (input$showTime) {
            p <- ggplot(df,aes(x=as.factor(timepoint),y=value+1,group=sample,color=sample,alpha=highlight)) + 
                  geom_line() + 
                  scale_alpha_manual(values=c(1,0.3),guide = FALSE)
          } else {
            p <- ggplot(df,aes(x=as.numeric(as.character(timepoint)),y=value+1,group=sample,color=sample,alpha=highlight)) + geom_line() + 
              scale_alpha_manual(values=c(1,0.3),guide = FALSE)
          }
          
          if ("frag" %in% colnames(df)) {
            p <- p + geom_point(aes(shape=frag), size=3) + scale_shape_manual(values = c(16,1))
          } else {
            p <- p + geom_point(size=3)
          }
          
          p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10), legend.position='bottom') + xlab('Timepoint')
          p = p + geom_errorbar(aes(ymin=mmin,ymax=mmax),width=.1)
          if (input$showLog) {
            p <- p + scale_y_log10()
          }
          if (is.null(dataset$fpkm)) {
            p <- p + ylab('Norm. count') 
          } else {
            p <- p + ylab('FPKM') 
          }
          
          p = p + facet_grid(type ~ syn,scales='free_y') + theme_bw()
        }
    } else {
      df$s = paste(df$sample,df$timepoint,sep='_')
      df$s = factor(df$s,levels=unique(df$s),ordered = T)
      ### Get the order of the labels, so that things cluster nicely together.
      ddist = dcast(df, syn ~ s)
      rownames(ddist) = ddist$syn
      ddist = ddist[,2:ncol(ddist)]
      ddist = hclust(dist(log10(ddist+1)))
      ord = ddist$labels[ddist$order]
      ## Force this order
      df$syn = factor(df$syn,levels=ord,ordered=T)
      p <- ggplot(df,aes(x=s,y=syn,fill=value)) + geom_tile() + xlab('Sample') + ylab('Gene')
      if (input$showLog) {#Log this stuff too
        p <- p + scale_fill_gradientn( trans = "log", colors = c('white','red'))
      } else {
        p <- p + scale_fill_gradientn( colors = c('white','red'))
      }
      p = p + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=10), 
                                 axis.text.y = element_text(size=10))
      p = p + facet_grid(. ~ sample,scales='free_x')
    }
    p + theme(strip.text.x = element_text(size = 15,angle=90),axis.text = element_text(size=15))
})

observeEvent(eventExpr = input$plot1_click, handlerExpr = {
      d <- dataToPlot()
      rr <- getPoint(input$plot1_click,d)
      if (nrow(rr)==1) {      
			if (length(highlight$clicked)==1) {
			  highlight$clicked <- rep(TRUE,nrow(d))
			}
			highlight$clicked[which(d$Gene==rr$Gene)] <- !highlight$clicked[which(d$Gene==rr$Gene)]
      }
      
})

output$hover_info <- renderUI({
    hover <- input$plot1_hover
    if (is.null(hover$x)) return()
    d <- dataToPlot()
    if (nrow(d) < 1) return()
    rr <- getPoint(hover,d)
    
    if (nrow(rr) < 1) return()
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - log10(hover$y)) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px;")
    
    # actual tooltip created as wellPanel
    rr$Annotation <- allTable()[which(allTable()$Gene==rr$Gene),]$Annotation
    wellPanel(
      style = style,
      p(HTML(paste0(rr$syn,"<br/>",rr$Annotation,"<br/>")))
    )
})

getJbrowseLink <- function (genome, chromosome, range_start, range_end, hight_start, 
                            hight_end, tracks = c("DNA", "gene_model"), baseurl, 
                            show_navagation = TRUE, show_tracklist = FALSE, show_overview = FALSE) 
{
  baseurl_ <- baseurl
  genome_ <- genome
  chromosome_ <- chromosome
  range_ <- if (missing(range_start) || missing(range_end)) 
    ""
  else paste0(":", parseRange(start = range_start, end = range_end, 
                              resizeFactor = 1.5))
  highlight_ <- if (missing(hight_start) || missing(hight_end)) 
    ""
  else paste0("&highlight=", chromosome, ":", parseRange(start = hight_start, 
                                                         end = hight_end, resizeFactor = 1))
  tracks_ <- paste(unique(tracks), collapse = ",")
  navagation_ <- if (show_navagation) 
    ""
  else "&nav=0"
  tracklist_ <- if (show_tracklist) 
    ""
  else "&tracklist=0"
  overview_ <- if (show_overview) 
    ""
  else "&overview=0"
  sprintf("%s/?data=references/%s&loc=%s%s&tracks=%s%s%s%s%s", baseurl_, 
          genome_, chromosome_, range_, tracks_, highlight_, navagation_, 
          tracklist_, overview_)
}

output$jBrowse_cov <- renderJbrowse({
  smples = dataset$meta[input$samples_rows_selected,]
  toUse = rownames(smples)
  all <- allTable()
  selectedGene = all[input$genesInRange_rows_selected,]$Gene[1]
  location = dataset$geneCoordinates[selectedGene,]
  chromosome = strsplit(as.character(location[1]),";")[[1]][1]
  start = as.numeric(strsplit(as.character(location[2]),";")[[1]][1])
  end = as.numeric(strsplit(as.character(location[3]),";")[[1]][1])
  url <- getJbrowseLink(genome = dataset$jbrowse, chromosome = chromosome, range_start=start, range_end=end, baseurl=dataset$jbrowse_server,
                        tracks = c("Annotation",toUse), show_navagation = TRUE, show_tracklist = FALSE, show_overview = FALSE)
  print(url)
  iframeJbrowse(url)
})
