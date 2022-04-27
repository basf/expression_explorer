computedNetworks <- reactive({
    
    place = paste(baseFolder,dataset$loaded,'availableNetworks.tab',sep='/')
    #Check for all computed correlation networks
    if (file.exists(place)) {
        df <- read.table(place,header=T,row.names=1,check.names=F,sep='\t',stringsAsFactors=F)
    } else {#ok, this probably only has the standard stuff that we historically precomputed so let's make this file
        df <- data.frame(row.names=1,fileDesignator='',samplesUsed='All',minNegative=-0.95,minPositive=0.95,name='All')
        write.table(df,place,quote=F,sep='\t')
    }

    df$inUse = FALSE

    if (is.na(highlight$networkUsed)) {
        df$inUse[which(is.na(df$fileDesignator))] = TRUE
    } else {
        df$inUse[which(df$fileDesignator==highlight$networkUsed)] = TRUE
    }
    
    return(df)
})

output$network_list <- DT::renderDataTable({
    t <- allTable()
    validate(need(!is.null(t),"You should load a dataset first!"))
    t
  },server=TRUE,
    rownames=FALSE,
    filter = list(position = 'top',clear=FALSE)
   )

proxy1 = dataTableProxy('network_list')
 
observeEvent(input$clear2, {
    proxy1 %>% selectRows(NULL)
})

genesToShow <- reactive({
    validate(need(length(input$network_list_rows_selected)>0,"Select some genes below, see what happens..."))
    genes <- allTable()[input$network_list_rows_selected,]$Gene
    genes
})
 
output$net_genes <- renderText({
    gg <- genesToShow()
    eData <- subset(dataset$edgeData,source %in% gg | target %in% gg)
    nData <- subset(dataset$nodeData,name %in% gg | name %in% eData$source | name %in% eData[["target"]])
    
    as.character(nData$name)
})
 
output$network <- renderCyjShiny({
  
    genes <- genesToShow()
    eData <- subset(dataset$edgeData,source %in% genes | target %in% genes)
    nData <- subset(dataset$nodeData,name %in% genes | name %in% eData$source | name %in% eData[["target"]])
    nData$tooltip = paste(nData$name,allTable()[as.character(nData$name),]$Annotation)
    
    #There is no alpha support
    #nData$color = substr(map2color(log10(allTable()[nData$name,][['Average Normalized Count']]),rainbow(20)),1,7)
    
    if (nrow(eData) > 0) {
      eData$color = 'blue'
      eData$color[eData$score<0] = 'red'
    } else {
      eData <- data.frame(matrix(ncol=3,nrow=0))
      colnames(eData) = c('source','target','score')
    }
    
    eData$source = as.character(eData$source)
    eData$target = as.character(eData$target)
    eData$interaction = eData$score
    
    nData$id = as.character(nData$id)
    nData$name = as.character(nData$name)
    rownames(nData) = nData$name
    nData$tooltip = as.character(nData$tooltip)
    
    graph.json <- dataFramesToJSON(eData,nData)
    net = cyjShiny(graph=graph.json, layoutName="cola", style_file="graph_style.js")
               
    return(net)
    
})

observeEvent(input$deleteNet, {

    row_to_del <- as.numeric(gsub("Row","",input$availableNetworks_rows_selected))
    net <- computedNetworks()
    run <- net$fileDesignator[row_to_del]
    namesRun <- net$name[row_to_del]

    if(net$samplesUsed[row_to_del] != "All"){
      rm_res <- file.remove(paste(baseFolder,dataset$loaded,paste('edgeData', run,'_',namesRun,'.tab',sep=''),sep='/'))
      
      if (rm_res) {
        net <- net[-row_to_del,]
        write.table(net,paste(baseFolder,dataset$loaded,'availableNetworks.tab',sep='/'),quote=F,sep='\t')  
        # this triggers computedNetworks to update the table
        highlight$networkUsed <- run
        
      } else {
        # show some rror? e.g. shinytoastr or shinyBS::altert
        showModal(modalDialog(
          title = "This file doesn't exist",
          footer = tagList(
            modalButton("Close")
          )
        ))}
      }
})

output$availableNetworks <- DT::renderDataTable({
    net = computedNetworks()
    df = net[,c("samplesUsed", "minNegative", "minPositive", "inUse")]
    df
    
},rownames=FALSE,server=TRUE,selection = list(mode = 'single'))
  
output$samplesNetwork <- DT::renderDataTable({
    validate(need(!is.null(dataset$meta),message=FALSE))
    dd <- dataset$meta[,c('Sample','Timepoint','Replicate')]
    dd$Timepoint <- as.character(dd$Timepoint)
    dd
  },server=FALSE,rownames=FALSE,filter = list(position = 'top',clear=TRUE))
 
observeEvent(input$useNet, {
    if (length(input$availableNetworks_rows_selected) != 1) {
        showModal(modalDialog(
            title = "Network selection",
            "Please select a network you would like to visualize",
            footer = tagList(
                modalButton("Close")
            )
        ))
    } else {
        net <- computedNetworks()
        run <- net$fileDesignator[input$availableNetworks_rows_selected]
        
        if(!is.null(net$name[input$availableNetworks_rows_selected])){
          namesRun <- net$name[input$availableNetworks_rows_selected]
        }else{
          namesRun <- revertReadableSet(net$samplesUsed[input$availableNetworks_rows_selected])
          net$name[input$availableNetworks_rows_selected] <- namesRun
        }
        if (is.na(run)) {
            run <- '';
        }
        
        withProgress(message= 'Loading network...', value=0.5, {
          dataset$edgeData <<- read.table(paste(baseFolder,dataset$loaded,paste('edgeData', run,'_',namesRun,'.tab',sep=''),sep='/') ,header=T,sep='\t',quote='')
        })
        
        highlight$networkUsed <- run
    }
})
 
observeEvent(input$computeNet, {
    
    if (length(input$samplesNetwork_rows_selected) < 2) {
        showModal(modalDialog(
            title <- "Correlation error",
            "Please select at least two samples in order to perform correlation analysis",
            footer = tagList(
                modalButton("Close")
            )
        ))
    } else {
        minCor <- input$minCorr
        maxCor <- input$maxCorr
        
        withProgress(message= 'Computing network...', value=0.5, {
            final <- sqrt(dataset$norm)
            use <- rownames(dataset$meta)[input$samplesNetwork_rows_selected]
            final <- final[,use]
            #Remove low abundant stuff
            final <- final[which(rowSums(final)>10),]
            #Get correlation matrix
            c <- cor(t(final),method='pearson')
            incProgress(2/10)
            
            #Make list of all paths in this
            is = (c >= maxCor | c <= minCor)
            is[lower.tri(is,diag=T)] <- NA
            
            idx = which(is,arr.ind=T)
            
            edgeData <- data.frame(source=rownames(idx), target=rownames(is)[idx[,2]], score=c[which(is)], stringsAsFactors=FALSE)
            incProgress(2/10)
            
            net <- computedNetworks()
            if(length(net$fileDesignator) == 2){
              run <- 1
            }else{
              run <- (max(na.omit(net$fileDesignator)) + 1)
            }

            namesRun <- paste(use,collapse=',')# net$samplesUsed

            write.table(edgeData,paste(baseFolder,dataset$loaded,paste('edgeData', run,'_',namesRun,'.tab',sep=''),sep='/') ,quote=F,sep='\t')
            
            incProgress(1/10)
            
            net <- rbind(net,c(run,getReadableSet(paste(use,sep=',')),minCor,maxCor,namesRun, FALSE))
            #Update the compute networs file
            write.table(net,paste(baseFolder,dataset$loaded,'availableNetworks.tab',sep='/'),quote=F,sep='\t')
            highlight$networkUsed <- run
        })
        
    }
})
