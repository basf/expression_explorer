shinyServer(function(input,output,session) {
  observe({
    args <- parseQueryString(session$clientData$url_search)
    baseFolder <- ifelse(Sys.getenv("RSTUDIO") == "1", "./data", "/data")
    dir_exists <- fs::dir_exists(fs::path(baseFolder, args$name))
  
      if(isTruthy(dir_exists)) {
      res <- loadDataset(args$name)
      if (is.list(res)) {
        dataset$norm <<- res$norm
        dataset$fpkm <<- res$fpkm
        dataset$ann <<- res$ann
        dataset$edgeData <<- res$edgeData
        dataset$nodeData <<- res$nodeData
        dataset$meta <<- res$meta
        dataset$expression <<- res$expression
        dataset$availDiffs <<- res$diff
        dataset$diffs <<- res$diff
        dataset$loaded <<- args$name
        dataset$folder <<- res$folder
        dataset$jbrowse <<- res$jbrowse
        dataset$prot <<- res$prot
        dataset$span <<- res$span
        dataset$mutations <<- res$mutations

        updateSelectInput(session,'pca_color',choices=colnames(dataset$meta))
        updateNavbarPage(session,'nPage',selected='Dataset overview')
        updateCheckboxGroupInput(session,"metaShow",choices=colnames(dataset$meta),selected=c("Sample","Timepoint","Replicate"))
        
      }
    }
  })

  #### Some stuff we'll use in this session
  
  dataset = reactiveValues(now = 0,
    meta = NULL,
    norm = NULL,
    edgeData = NULL,
    nodeData = NULL,
    ann = NULL,
    expression = NULL,
    fpkm = NULL,
    c = NULL,
    diffs = NULL,
    availDiffs = NULL,
    jbrowse=NULL,
    geneSets=NULL,
    sampleSets=NULL,
    ecs = NULL,
    span = NULL,
    mutations=NULL)
  
  ###### Loads stuff for the overview tab
  source(file.path("overview.R"),  local = TRUE)$value
  
  ###### Loads stuff for the differential expression tab
  source(file.path("diffExpr.R"),  local = TRUE)$value
  
  ###### Loads stuff for the network bits
  source(file.path("network_view.R"),  local = TRUE)$value
  
  ###### Loads stuff for metadata
  source(file.path("metadata.R"),  local = TRUE)$value
  
  
  output$dataset_list <- DT::renderDataTable({
    ### Just get all the samples info we have
    allData
  },rownames=FALSE,selection='single')
    
  observeEvent(input$load, {
    #Get selected folder
    sel <- input$dataset_list_rows_selected
    if (length(sel) != 1) {
      showModal(modalDialog(
		  title = "Select dataset you want to load",
		  footer = tagList(
		    modalButton("Close")
		  )
	))
    } else {
      withProgress(message= 'Loading dataset...', value=0.5, {
        res <- loadDataset(allData[sel,1])
      })
      if (is.list(res)) {
      	dataset$norm <<- res$norm
      	dataset$ann <<- res$ann
      	dataset$edgeData <<- res$edgeData
      	dataset$nodeData <<- res$nodeData
      	dataset$meta <<- res$meta
      	dataset$expression <<- res$expression
      	dataset$fpkm <<- res$fpkm
      	dataset$availDiffs <<- res$diff
      	dataset$diffs <<- res$diff
        dataset$loaded <<- allData[sel,1]
      	dataset$folder <<- res$folder
      	dataset$jbrowse <<- res$jbrowse
      	dataset$geneSets <<- res$geneSets
      	dataset$sampleSets <<- res$sampleSets
      	dataset$prot <<- res$prot
      	dataset$geneCoordinates <<- res$geneCoordinates
      	dataset$span <<- res$span
      	dataset$mutations <<- res$mutations
      	
      	updateSelectInput(session,'pca_color',choices=colnames(dataset$meta))
      	updateSelectInput(session,'pre_def_gene_list',choices=c("",unique(dataset$geneSets$Set)),selected="")
      	updateNavbarPage(session,'nPage',selected='Dataset overview')
      	updateCheckboxGroupInput(session,"metaShow",choices=colnames(dataset$meta),selected=c("Sample","Timepoint","Replicate"))
      	toggleDropdownButton("dropSampleSelect",session)
      	
      } else {
	showModal(modalDialog(
		  title = "Sorry, something went wrong when loading data...",
		  res,
		  footer = tagList(
		    modalButton("Close") 
		  )
	))
      }
    }
    
  })
  
  ##### This is the enrichment tab
  
  output$Enrich_res_table <- DT::renderDataTable({
    t <- enrich$computed
    t$Description = goDescriptions[as.character(t$Term),] #Overwrite the description to make more sense
    validate(need(!is.null(t),"No enrichment to show. Go compute one on a differential expression analysis"))
    t
  },server=TRUE,
    rownames=FALSE,
    filter = list(position = 'top',clear=FALSE),
   )
   
  output$Enrich_res_go <- DT::renderDataTable({
    t <- enrich$annotations
    validate(need(!is.null(t),NULL))
    t
  },server=TRUE,
    rownames=TRUE
   )
 
})


#shinyApp(ui,server)
