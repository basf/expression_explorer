output$dwn_ident <- downloadHandler(
    filename = function() { paste('Meta_identifiers', '.csv', sep='') },
    content = function(file) {
        df <- dataset$meta[,c('Sample'),drop=F]
        df$Unique_ID = rownames(df)
        write.table(df,file,quote=F,sep='\t',row.names=F)
    }
)

output$Metadata_table <- DT::renderDataTable({
        dd <- dataset$meta
        if (!is.null(dataset$sampleSets)) {
          m <- dataset$sampleSets
          m <- aggregate(Set ~ Sample, m, function(x) paste(x,collapse=' | '))
          rownames(m) <- m$Sample
          dd$Set <- m[rownames(dd),]$Set
        }
        dd
    },server=TRUE,rownames=TRUE,filter = list(position = 'top',clear=FALSE))
 
observe({
	inFile <- input$meta_file

	if (is.null(inFile))
		return(NULL)
		
	df = read.table(inFile$datapath,header=T,sep='\t',check.names=F)
	if ("Unique_ID" %in% colnames(df)) {
		rownames(df) = df$Unique_ID
	} else {
		showModal(modalDialog(
			title = "Error",
			"Sorry, file is not in the right format.",
			footer = tagList(
				modalButton("Close")
			)
		))
		return(NULL)
	}
	
	#It is very important that we isolate this stuff
	isolate({
		#Keep only stuff that makes sense
		df = df[which(rownames(df) %in% rownames(dataset$meta)),]
		#Remove the unique ID
		df = df[,-grep('Unique_ID',colnames(df)),drop=F]
		#Add this
		keep = colnames(df)[which(!colnames(df) %in% colnames(dataset$meta))]
		df = df[,keep,drop=F]
		dataset$meta = cbind(dataset$meta,df[rownames(dataset$meta),,drop=F])
	})
       
})
 
observeEvent(input$meta_save,{
    #Over write the metadata file!
    isolate({
        write.table(dataset$meta,paste(baseFolder,dataset$loaded,'meta.tab',sep='/'),sep='\t',quote=F)
    })
    showModal(modalDialog(
        title = "Message",
        "Metadata successfully updated",
        footer = tagList(
            modalButton("Close")
        )
    ))
})

observeEvent(input$saveSampleSet, {
  setName = input$sampleSet_name
  if (setName == "") {
    showModal(modalDialog(
      title = "Info missing",
      "Please provide a name for this sample set.",
      footer = tagList(
        modalButton("Close")
      )
    ))
  } else {
    
    #Get all the samples
    samples = rownames(dataset$meta)[input$Metadata_table_rows_selected]
    if (length(samples)==0) {#There's nothing to save here
      return()
    }
    new = dataset$sampleSets
    if (!is.null(new)) {
      new = subset(dataset$sampleSets,Set != setName)#We're overwriting whatever was here before
      new = rbind(new,data.frame(Sample=samples,Set=setName))
    } else {
      new = data.frame(Sample=samples,Set=setName)
    }
    #Now write this
    write.table(new,paste(dataset$folder,'samplesets.tab',sep='/'),quote=F,sep='\t',row.names=F)
    #Update our object
    dataset$sampleSets = new
    updateTextInput(session,'set_name',value="")
  }
})
