
jscode <- "Shiny.addCustomMessageHandler('pathwayRedirect', function(url) {window.open(url);});"

ui <- navbarPage('Expression Explorer',id='nPage',
  tabPanel("Available datasets",
    DT::dataTableOutput("dataset_list"),
    actionButton('load','Load dataset')
  ),
  tabPanel("Dataset overview",
           fluidRow(
    column(4,dropdownButton(label = "Select samples...",inputId="dropSampleSelect",
                   icon = icon("table"),
                   circle=FALSE,
	  DT::dataTableOutput("samples"),
	  hr(),
	  fluidRow(
	    column(actionButton('clearSample', 'Clear Sel.'),width=4),
	    column(actionButton('selectFilteredSample', 'Sel. Filtered'),width=4),
	    column(actionButton('selectAllSample', 'Sel. All'),width=4)
	  ),
	  fluidRow(
		#column(bookmarkButton(),width=4)
		column(checkboxInput('showRawNames', 'Show sample names', value = FALSE, width = NULL),width=12)
	  )	  
	)),
	column(4,dropdownButton(label = "Select sets...",
	                        icon = icon("table"),
	                        circle=FALSE,
	                        DT::dataTableOutput("sample_sets")	  
	)),
	column(4,dropdownButton(label = 'Metadata columns...',
	               icon = icon("table"),
	               circle = FALSE,
	               checkboxGroupInput("metaShow", "Metadata to show:",
	                                  c("Sample","Timepoint"),selected=c("Sample","Timepoint"))
	  
	))),
	  br(),
	  bsCollapse(open=c('doOver'),multiple=FALSE,id='overPanel',
	    bsCollapsePanel("Samples Overview",value='doOver',
	      plotOutput('plot3',height='600px'),
	      hr(),
	      fluidRow(
			column(
				sliderInput('loads',"Number of top loading components:",min=10,max=50,value=10),
			width=4),
			column(selectInput('pca_color','Color by:',choices=c('sample','timepoint')),width=3),
			column(downloadButton('dwnPCA','Download PCA'),width=2),
			column(2,checkboxInput('show_pca_group','Highlight groups',value=F))
	      )
	    ),
	    bsCollapsePanel("Gene expression",value='extra',
	      div(
		style="position:relative",
		plotOutput('plot1',height='600px',hover=hoverOpts(id="plot1_hover",delay=100,delayType='debounce'),click='plot1_click'),
		uiOutput("hover_info"),
		  fluidRow(
		      column(checkboxInput('showLog','Scale y-axis log10',value=TRUE),width=2),
		      column(checkboxInput('showHeat','Show as heatmap',value=FALSE),width=2),
		      column(checkboxInput('showFlip','Flip',value=FALSE),width=2,
		             bsTooltip('showFlip','Show sample expression per gene')),
		      column(checkboxInput('showTime','Timepoints as factors',value=TRUE),width=2),
		      column(checkboxInput('hideProt','Hide proteins',value=FALSE),width=2,
		             bsTooltip('hideProt','If selected, will hide the protein expression values.')),
		      column(checkboxInput('splitRepl','Split replicates',value=FALSE),width=2,
		             bsTooltip('splitRepl','Do not aggregate replicates into the same panel. Will treat each replicate as a different sample. Useful for understanding underlying biological variance'))
		  )
	      ),
	      hr(),
	      fluidRow(
            column(3,fluidRow(column(12,downloadButton('dwnPlot','Download Fig.'))),
                    fluidRow(
                      column(6,downloadButton('dwnRawGene','Get Plot Data')),
                      column(6,downloadButton('dwnAllGene','Get All Counts'))
                    ),
                    fluidRow(column(12,selectInput('pre_def_gene_list','Gene sets:',choices=c('')))),
                    fluidRow(column(12,textAreaInput('gene', 'Gene list',width='100%',height='300px'))),
                    fluidRow(column(12,actionButton('clear1', 'Clear Selection'))),
                    fluidRow(column(4,actionButton('saveList','Save gene set')),
                             column(8,textInput('set_name',label=NULL,placeholder="Set name"))
                             )
                   ),
                    
            column(
				fluidRow(
					checkboxInput('updateOnSelect', 'Show mean expression values for each gene', value = FALSE, width = NULL)
				),
				fluidRow(DT::dataTableOutput("genesInRange")),width=9)
	      )
	  ),bsCollapsePanel("Coverage",value='covs',
	       JbrowseOutput("jBrowse_cov", width = "100%", height = "400px")
	   )
	)
  ),
  tabPanel('Differential expression',
      tags$head(tags$script(jscode)),
      bsCollapse(open=c('extra'),multiple=TRUE,id='diffPanel',
	bsCollapsePanel("Perform Differential Expression",value='doDiff',
	tags$html("Select samples you would like to compare."),
	br(),
	tags$html("One set on the left (i.e. Condition1) and one on the right (i.e. Condition2)."),
	br(),
	fluidRow(
	  column(DT::dataTableOutput("set1"),actionButton('toggleFilter1', 'Toggle Filtered'),actionButton('clearFilter1', 'Clear Selection'),
	    br(),br(),
	    textOutput("set1Selected"),
	  width=6),
	  column(DT::dataTableOutput("set2"),actionButton('toggleFilter2', 'Toggle Filtered'),actionButton('clearFilter2', 'Clear Selection'),
	    br(),br(),
	    textOutput("set2Selected"),
	  width=6)
	),
	hr(),
	fluidRow(
	  column(actionButton('getDiff','Calculate',icon("paper-plane"),style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),width=2)
	)),
	bsCollapsePanel("Deep dive...",value='extra',
	tags$html("All performed differential expression analyses will be saved and can be investigated here."),
	br(),
	tags$html("The number of differentially expressed genes is the number of siginificant differences after correction. That is, padj < 0.05"),
	br(),
	fluidRow(
	  column(10,DT::dataTableOutput('allDiffAnalysis')),
	  column(2,actionButton('deleteDiffAnalysis', 'Delete Analysis'))
	),	conditionalPanel(condition="input.allDiffAnalysis_rows_selected.length == 1",
	fluidRow(
	  column(DT::dataTableOutput("geneDiff"),width=8),
	  column(
	    fluidRow(plotOutput('diffExprPlot')),
	    fluidRow(checkboxInput('showHeat_diff','Show as heatmap',value=FALSE)),
	    fluidRow(downloadButton('dwnDiffBox','Figure'),downloadButton('dwnDiff','Table'),downloadButton('dwnDiffFilt','Filtered Table')),
	    fluidRow(actionButton('enrichA','Enrichment Analysis')),
	    width=4)
	    
	),
	fluidRow(
	  column(actionButton('clearDiff','Clear Selection'),width=2)
	),
	fluidRow(
	  column(plotOutput("diffHierarchicalClustering"),width=12)
	)
	),
	
	conditionalPanel(condition="input.allDiffAnalysis_rows_selected.length > 1",
	fluidRow(
	  column(plotOutput('vennPlot'),width=7),
	  column(plotOutput('overlapDiffPlot'),width=5)
	),
	downloadButton('dwnDiffVennPlot','Download Venn'),
	hr(),
	fluidRow(
	  column(sliderInput('qvalLimit','Significance threshold',min=0.01,max=0.2,value=0.05),width=4),
	  column(sliderInput('baseLimit','Min. counts',min=0,max=50000,value=10),width=4),
	  column(sliderInput('log2Limit','Absolute log2 fold difference',min=0,max=10,value=1),width=4)
	  #column(selectInput('venn_side','Group:',choices=c('overlap','left','right')),width=3)
	),
	fluidRow(
	  column(selectInput('diff_include','Include:',choices='',multiple=T),width=4),
	  column(selectInput('diff_exclude','Exclude:',choices='',multiple=T),width=4),
	  column(checkboxInput('diff_noSig', 'Ignore p-vals', value = FALSE, width = NULL),width=4)
	),
	hr(),
	fluidRow(
	  column(DT::dataTableOutput("overlappingGene"),width=12)
	),
	fluidRow(
	  actionButton('enrichB','Enrichment Analysis')
	),
	fluidRow(
	  column(plotOutput('allOverlapPlot'),width=12)
	),
	downloadButton('dwnDiffVennSub','Download Common Set'),
	downloadButton('dwnDiffVenn','Download Full Comparison')
	)
	)
      )
  ),
  tabPanel('Network View',
    bsCollapse(open=c('viewNetwork'),multiple=FALSE,id='nets',
        bsCollapsePanel("Select/Compute Network",value='selectNetwork',
            fluidRow(
                column(8,DT::dataTableOutput("availableNetworks")),
                column(1,actionButton('useNet','Use')),
                column(1,actionButton('deleteNet','Delete'))
            ),
            hr(),
            fluidRow(
                column(6,DT::dataTableOutput("samplesNetwork")),
                column(4,
                    fluidRow(
                        sliderInput('minCorr','Max. Negative Correlation',min=-1,max=-0.5,value=-0.95)
                    ),
                    fluidRow(
                        sliderInput('maxCorr','Min. Positive Correlation',min=0.5,max=1,value=0.95)
                    ),
                    hr(),
                    fluidRow(
                        actionButton('computeNet','Compute Network')
                    )
                )
            )
        ),
        bsCollapsePanel("View Network",value='viewNetwork',
            fluidRow(
            column(cyjShinyOutput("network", height="600px"),width=10),
            column(textOutput("net_genes"), width=2)
            ),
            hr(),
            DT::dataTableOutput("network_list"),
            hr(),
            actionButton('clear2', 'Clear Selection')
            #downloadButton('downNet', 'Download graph')
        )
    )
  ),
  
  tabPanel('Metadata',
    fluidRow(
        column(DT::dataTableOutput('Metadata_table'),width=12)
    ),
    fluidRow(column(4,actionButton('saveSampleSet','Save selection as sample set')),
             column(8,textInput('sampleSet_name',label=NULL,placeholder="Sample set name:"))
    ),
    fluidRow(
        column(3,fileInput("meta_file","Upload additional metadata",accept = c('text/csv','.csv'))),
        column(3,downloadButton("dwn_ident","Download identifiers")),
        column(3,actionButton("meta_save","Save metadata changes"))
    )
  ),
  tabPanel('Enrichment Analysis',
	fluidRow(
		column(DT::dataTableOutput('Enrich_res_table'),width=12)
	),
	fluidRow(
		column(DT::dataTableOutput('Enrich_res_go'),width=12)
	)
  )
)
