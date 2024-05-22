
server <- function(input, output, session) {
  
  ## Define max file size
  options(shiny.maxRequestSize=530*1024^2)
  #options(DT.options = list(pageLength = 15))
  options(DT.options = list(scrollX = TRUE))
  source('global.R', local = TRUE)
  
  
  ## Prepare validation Data
  gse87211_validate <- eventReactive(input$buttonValidateGse87211, ignoreNULL = TRUE, ignoreInit = FALSE, {
    if(input$choice=="Select"){
      shiny::validate(
        need(input$choice != "Select", "Please select balanced or full cohort to validate")
      )
    }else{
      cohort  <-  input$choice
    }
    shiny::validate(
      need(is.null(input$file1$datapath) != TRUE, "Please select expression data for validation cohort")
    )
    shiny::validate(
      need(is.null(input$file2$datapath) != TRUE, "Please select meta data for validation cohort")
    )
    expression_data <- read.csv(file = input$file1$datapath, header = TRUE, row.names = 1)
    metadata <- read.csv(file = input$file2$datapath, header = TRUE, row.names = 1)
    nar <- FALSE
    TN_stage <- FALSE
    TN_data <- "factor"
    val_cohort_size <- input$choice
    load(paste0(path_to_discovery_data,"DiscoveryBig.RData"))
    results <- validate_GSE87211(as.matrix(expression_data), metadata, val_cohort_size, combat_type)
    accuracy_table <- as.data.frame(t(results[[3]]$overall))
    return(list(results[[1]], results[[2]], accuracy_table, results[[4]]))
  })

  output$predictionprobs = DT::renderDataTable(DT::datatable({gse87211_validate()[[4]][,c(1:4)]}
                                                             , extensions = 'Buttons'
                                                             , options = list( 
                                                               dom = "Blfrtip"
                                                               , buttons = 
                                                                 list("copy", list(
                                                                   extend = "collection"
                                                                   , buttons = c("csv", "excel", "pdf")
                                                                   , text = "Download"
                                                                 ) ) # end of buttons customization
                                                               
                                                               # customize the length menu
                                                               , lengthMenu = list( c(10, 20, -1) # declare values
                                                                                    , c(10, 20, "All") # declare titles
                                                               ) # end of lengthMenu customization
                                                               , pageLength = 10
                                                               
                                                             ) # end of options
                                                             
  ), server = FALSE
  )
  output$predictionaccuracy = DT::renderDataTable(DT::datatable(gse87211_validate()[[3]], escape=FALSE,
                                                             options = list(pageLength = 5, autoWidth = TRUE,
                                                                            columnDefs = list(list( targets = c(6), width = '200px')),
                                                                            scrollX = TRUE)))

  output$probPlot <- renderPlot({
        return(gse87211_validate()[[1]])
      })
  output$rocPlot <- renderPlot({
    return(gse87211_validate()[[2]])
  })
  
  
  ## Predict label for validation Data
  gse87211_predict <- eventReactive(input$buttonPredictGse87211, ignoreNULL = TRUE, ignoreInit = FALSE, {
    if(input$choice=="Select"){
      shiny::validate(
        need(input$choice != "Select", "Please select balanced or full cohort to validate")
      )
    }else{
      cohort  <-  input$choice
    }
    shiny::validate(
      need(is.null(input$file1$datapath) != TRUE, "Please select expression data for validation cohort")
    )
    shiny::validate(
      need(is.null(input$file2$datapath) != TRUE, "Please select meta data for validation cohort")
    )
    expression_data <- read.csv(file = input$file1$datapath, header = TRUE, row.names = 1)
    metadata <- read.csv(file = input$file2$datapath, header = TRUE, row.names = 1)
    nar <- FALSE
    TN_stage <- FALSE
    TN_data <- "factor"
    val_cohort_size <- input$choice
    load(paste0(path_to_discovery_data,"DiscoveryBig.RData"))
    results <- predict_GSE87211(as.matrix(expression_data), metadata, val_cohort_size, combat_type)
    accuracy_table <- "NA"
    return(list(results[[1]], results[[2]], accuracy_table, results[[4]]))
  })
  
  output$predictionprobs2 = DT::renderDataTable(DT::datatable({gse87211_predict()[[4]][,c(1:4)]}
                                                             , extensions = 'Buttons'
                                                             , options = list( 
                                                               dom = "Blfrtip"
                                                               , buttons = 
                                                                 list("copy", list(
                                                                   extend = "collection"
                                                                   , buttons = c("csv", "excel", "pdf")
                                                                   , text = "Download"
                                                                 ) ) # end of buttons customization
                                                               
                                                               # customize the length menu
                                                               , lengthMenu = list( c(10, 20, -1) # declare values
                                                                                    , c(10, 20, "All") # declare titles
                                                               ) # end of lengthMenu customization
                                                               , pageLength = 10
                                                               
                                                             ) # end of options
                                                             
  ), server = FALSE
  )
}
