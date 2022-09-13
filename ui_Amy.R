## Load R packages
library(shiny)
library(shinythemes)
library(markdown)


## set working directory
#setwd("~/Desktop/NEW_ChROseq_HIO/miRNA_effect_2021/shiny")


## Load required data
# gene name list include genes with mean of RNA nc > 50 and mean of ChRO TPM > 10 in each of the comparison
gc <- scan('./data/geneNameChoices_HIO.csv', as.character(), sep = '\n')


## miR list include miRs with mean RPMMM > 50 across all stages
mc <- scan('./data/miRNameChoices_HIO.csv', as.character(), sep = '\n')


## Define UI -------------------------------------------------------
UI <- fixedPage(
                theme = shinytheme("simplex"),
                navbarPage(
                  "Post-transcriptional Regulation of Human Intestinal Development (ProHumID)",
                  
                  tabPanel("Introduction", 
                           
                           fluidRow(
                             column(7, offset = 0,
                                    includeMarkdown("md/intro_Amy.md"),
                                    hr()),
                             
                             column(3,
                                    br(), 
                                    br(),
                                    br(),
                                    br(),
                                    img(src='graphical_abstract.png',width="430",height="500"))
                           )
                  ), # tabPanel
                  
                  tabPanel("Data exploration", 
                           
                           sidebarPanel(width = 3, 
                                        style="min-width:300px; max-width:300px", 
                                        
                                        selectInput('moi', 
                                                    HTML('<b>miRNA of interest:<b>'), 
                                                    mc, 
                                                    multiple = F,
                                                    selected = 'miR-182-5p'), 
                                        
                                        selectInput('goi', 
                                                    HTML('<b>Gene of interest:<b>'), 
                                                    gc, 
                                                    multiple = F,
                                                    selected = 'SOX2'),
                                        
                                        h5(HTML('<b>Is the gene a target of the miRNA?<b>')),
                                        h4(textOutput('mirtarget')),
                                        
                                        submitButton("Update View")
                                        
                           ), # sidebarPanel
                           
                           mainPanel(width = 12,
                                     h5(HTML('<b>Plots and statisitcal significance are summarized below:<b>')),
                                     br(),
                                     #HTML('<center><img src="directed_differentiation.png" style="width: 300px"><center>'),
                                     
                                     tabsetPanel(
                                       id = 'dataset',
                                       
                                       tabPanel(
                                         'Definitive endoderm formation', 
                                         
                                         
                                         fixedRow(
                                           div(style = "padding: 10px 0px",
                                           
                                           column(3, 
                                                  div(style = "height:320px; width:250px; background-color: white", 
                                                      plotOutput('expression_miRNA1', width = "250px", height = "250px"),
                                                      h5("Level of miRNA is:"),
                                                      verbatimTextOutput("miR_summary1")
                                                  )
                                           ), 
                                           
                                           column(3, offset = 0, 
                                                  div(style = "height:320px; width:250px; background-color: white", 
                                                      plotOutput('expression_ChRO1', width = "250px", height = "250px"),
                                                  h5("Transcription activity of gene is:"),
                                                  verbatimTextOutput("ChRO_summary1")
                                                  )
                                           ),
                                           
                                           column(3, offset = 0,
                                                  div(style = "height:320px; width:250px; background-color: white", 
                                                      plotOutput('expression_RNA1', width = "250px", height = "250px"),   
                                                  h5("Steady state level of gene is:"),
                                                  verbatimTextOutput("RNA_summary1")
                                                  )
                                           ), 
                                           
                                           column(3, offset = 0, 
                                                  div(style = "height:320px; width:250px; background-color: white", 
                                                      plotOutput('scatter_1', width = "250px", height = "250px"),
                                                  h5("Post-transcriptional status of gene is:"),
                                                  verbatimTextOutput("tf_summary1")
                                                  )
                                           )
                                           
                                        ) # divstyle
                                    ), # fixedRow
                                         
                                         fixedRow(
                                           column(12,
                                                  div(style = "padding: 10px 0px",
                                                      
                                                      h5("Is the miR:target relationship likely present in this stage transition? (based on smRNA-seq and two-factor analysis)"),
                                                      verbatimTextOutput("axis1"), 
                                                      br()
                                                     ) # div style
                                                  ) # column
                                              ) # fixedRow
                                         
                                           ), # tabPanel
                                       
                                       tabPanel(
                                         'Intestinal lineage specification',
                                         
                                         fixedRow(
                                           div(style = "padding: 10px 0px",
                                               
                                               column(3, 
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('expression_miRNA2', width = "250px", height = "250px"),
                                                          h5("Level of miRNA is:"),
                                                          verbatimTextOutput("miR_summary2")
                                                      )
                                               ), 
                                               
                                               column(3, offset = 0, 
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('expression_ChRO2', width = "250px", height = "250px"),
                                                          h5("Transcription activity of gene is:"),
                                                          verbatimTextOutput("ChRO_summary2")
                                                      )
                                               ),
                                               
                                               column(3, offset = 0,
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('expression_RNA2', width = "250px", height = "250px"),   
                                                          h5("Steady state level of gene is:"),
                                                          verbatimTextOutput("RNA_summary2")
                                                      )
                                               ), 
                                               
                                               column(3, offset = 0, 
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('scatter_2', width = "250px", height = "250px"),
                                                          h5("Post-transcriptional status of gene is:"),
                                                          verbatimTextOutput("tf_summary2")
                                                      )
                                               )
                                               
                                           ) # div style
                                         ), # fixeddRow
                                         
                                         fixedRow(
                                           column(12,
                                                  div(style = "padding: 10px 0px",
                                                      
                                                      h5("Is the miR:target relationship likely present in this stage transition? (based on smRNA-seq and two-factor analysis)"),
                                                      verbatimTextOutput("axis2"), 
                                                      br()
                                                  ) # div style
                                           ) # column
                                         ) # fixedRow
                                         
                                       ), # tabPanel
                                       
                                       tabPanel(
                                         'Intestinal regional specification',
                                         
                                         fixedRow(
                                           div(style = "padding: 10px 0px",
                                               
                                               column(3, 
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('expression_miRNA3', width = "250px", height = "250px"),
                                                          h5("Level of miRNA is:"),
                                                          verbatimTextOutput("miR_summary3")
                                                      )
                                               ), 
                                               
                                               column(3, offset = 0, 
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('expression_ChRO3', width = "250px", height = "250px"),
                                                          h5("Transcription activity of gene is:"),
                                                          verbatimTextOutput("ChRO_summary3")
                                                      )
                                               ),
                                               
                                               column(3, offset = 0,
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('expression_RNA3', width = "250px", height = "250px"),   
                                                          h5("Steady state level of gene is:"),
                                                          verbatimTextOutput("RNA_summary3")
                                                      )
                                               ), 
                                               
                                               column(3, offset = 0, 
                                                      div(style = "height:320px; width:250px; background-color: white", 
                                                          plotOutput('scatter_3', width = "250px", height = "250px"),
                                                          h5("Post-transcriptional status of gene is:"),
                                                          verbatimTextOutput("tf_summary3")
                                                      )
                                               )
                                               
                                           ) # div style
                                         ), # fluidRow
                                         
                                         fixedRow(
                                           column(12,
                                                  div(style = "padding: 10px 0px",
                                                      
                                                      h5("Is the miR:target relationship likely present in this stage transition? (based on smRNA-seq and two-factor analysis)"),
                                                      verbatimTextOutput("axis3"), 
                                                      br()
                                                  ) # div style
                                           ) # column
                                         ) # fixedRow
                                         
                                       ) # tabPanel
                                    ) # tabsetPanel         
                           ) # mainPanel
                           
                  ))) # fluidPage

