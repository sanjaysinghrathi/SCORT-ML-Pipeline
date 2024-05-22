library(shiny)
library(DT)

ui <- fluidPage(
  theme = "bootstrap.css",
  hr(),
  titlePanel(title=div("Colorectal Cancer Stratification App by S:CORT",h6("Developed by Sanjay Rathee & Maintained by Andrew")), windowTitle = ""),  
  hr(),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Upload Validation Cohort"),
        hr(),
        # div(style="display: inline-block;vertical-align:top; width: 47%;",textInput("searchword", " Search Drug Name", value = "", width = NULL, placeholder = NULL)),
        #div(style="display: inline-block;vertical-align:top; width: 47%;",checkboxInput("checkGSE87211",value = TRUE, label = "Select GSE87211")),
        #hr(),
        #h4("Or"),
        fileInput("file1", "Upload Expression CSV File", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        fileInput("file2", "Upload Metadata CSV File", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        hr(),
        div(style="display: inline-block;vertical-align:top; width: 94%;",selectInput('choice', 'Select cohort size', selected = "Select" ,choices = c("Select","Full Cohort", "Balanced Cohort"))),
        div(style="display: inline-block;vertical-align:top; width: 47%;",actionButton("buttonValidateGse87211", "Validate Cohort")),
        div(style="display: inline-block;vertical-align:top; width: 47%;",actionButton("buttonPredictGse87211", "Predict Cohort")),
      ),
      wellPanel(
        h5("Disclaimer: We don not store any data uploaded here."),
        hr(),
        h5("All rights reserved to S:CORT"),
      ),
        
    ),
  
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Instructions",
                           h3(""),
                           h4("Please prepare the validation data as instructed:"),
                           hr(),
                           h5("1. Expression data csv file should have samples as columns and genes (Entrez ID only) as rows."),
                           h5("2. Metadata csv file should have samples as row and attributes(response, t-stage, n-stage) as column."),
                           h5("3. Metadata csv file for predictions should have only two attributes (t-stage and n-stage) as column."),
                           h5("4. Both csv files should have row and col names."),
                           h5("5. Put all t-stage and n-stage values as 1 if no data available."),
                           hr(),
                           #DT::dataTableOutput('sample_exp', width = "99%", height = "20em") %>% withSpinner(color="#0dc5c1"),
                           hr(),
                  ),
                  tabPanel("Validation Plots",
                           h3(""),
                           h5("Plots to show prediction probabilities and roc curve for  uploaded cohort. Please confirm if sample names in columns of expression data are same as sample name in row names of metadata. Metadata should have three columns as response, t-stage, n-stage."),
                           hr(),
                           plotOutput("probPlot", height = "20em") %>% withSpinner(color="#0dc5c1"),
                           hr(),
                           plotOutput("rocPlot", height = "20em") %>% withSpinner(color="#0dc5c1"),
                  ),
                  tabPanel("Validation Tables",
                           # column(8, wellPanel(
                           h5("Plots to show prediction probabilities and roc curve for  uploaded cohort Please confirm if sample names in columns of expression data are same as sample name in row names of metadata. Metadata should have three columns as response, t-stage, n-stage."),
                           hr(),
                           DT::dataTableOutput('predictionprobs', width = "99%", height = "20em") %>% withSpinner(color="#0dc5c1"),
                           hr(),
                           DT::dataTableOutput('predictionaccuracy', width = "99%", height = "10em") %>% withSpinner(color="#0dc5c1"),
                  ),
                  tabPanel("Prediction Tables",
                           # column(8, wellPanel(
                           h5("Plots to show prediction probabilities and roc curve for uploaded cohort. Please confirm if sample names in columns of expression data are same as sample name in row names of metadata. Metadata should have two columns as  t-stage, n-stage."),
                           hr(),
                           DT::dataTableOutput('predictionprobs2', width = "99%", height = "20em") %>% withSpinner(color="#0dc5c1"),
                  )
      )
      
    ),
    
  )
  
  
  
)
