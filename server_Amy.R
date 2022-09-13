# Load required packages
library(tidyverse)
library(ggpubr)
library(ggrepel)


## Load UI component
source('./ui_Amy.R')


## Load required data
# ChRO-seq expression
tpm1 <- readRDS('./data/ChROseq_TPM_noTRE_DEvshESC.rds')
tpm2 <- readRDS('./data/ChROseq_TPM_noTRE_DuovsDE.rds')
tpm3 <- readRDS('./data/ChROseq_TPM_noTRE_IlevsDuo.rds')

# RNA-seq expression
nc1 <- readRDS('./data/RNAseq_nc_DEvshESC.rds')
nc2 <- readRDS('./data/RNAseq_nc_DuovsDE.rds')
nc3 <- readRDS('./data/RNAseq_nc_IlevsDuo.rds')

# two-factor 
tf1 <- readRDS('./data/TF_DESeq2_noTRE_DEvshESC.rds')
tf2 <- readRDS('./data/TF_DESeq2_noTRE_DuovsDE.rds')
tf3 <- readRDS('./data/TF_DESeq2_noTRE_IlevsDuo.rds')

# small RNA-seq expression
rpmmm <- readRDS('./data/smRNAseq_all_stages_nc.rds')

# DESeq2 summary
miR_summary <- readRDS('./data/smRNAseq_DESeq2_summary.rds')
RNA_summary <- readRDS('./data/RNAseq_DESeq2_summary.rds')
ChRO_summary <- readRDS('./data/ChROseq_DESeq2_summary.rds')


## Load functions
# mir target function 
mir2mirFam <- readRDS('./data/miR2miR_Family.rds')

# Source functions
source("./functions/subsetDF_2.R")
source("./functions/mirTarget.R")
source("./functions/plotting.R")


# establish shiny server -------------------------------------------------------
shinyServer <- function(input, output) {
  
  ###############################################
  #
  #    Query miR:target result
  #
  ###############################################
  
  mirtarget <- reactive({
    
    mirTarget.f(input$moi, input$goi, mir2mirFam)
    
  })
  
  output$mirtarget <- renderText({
    mirtarget()
  })
  
  
  ############################################### 
  #
  #    DE formation 
  #
  ############################################### 
  ## ChRO-seq 
  tpm1_out <- reactive({
    subsetGene_DF.f(tpm1, input$goi, "DE", "hESC")
  })
  
  output$expression_ChRO1 <- renderPlot({
    bar.f(tpm1_out(),
          "ChRO-seq", 
          bquote(italic(paste(.(input$goi))) ~ paste(' (TPM)')))
  })
  
  ## RNA-seq 
  nc1_out <- reactive({
    subsetGene_DF.f(nc1, input$goi, "DE", "hESC")
  })
  
  output$expression_RNA1 <- renderPlot({
    bar.f(nc1_out(),
          "RNA-seq", 
          bquote(italic(paste(.(input$goi))) ~ paste(' (nc)')))
  })
  
  
  ## smRNA-seq 
  rpmmm1_out <- reactive({
    subsetmiR_DF.f(rpmmm, input$moi)
  })
  
  output$expression_miRNA1 <- renderPlot({
    bar.mir.f(rpmmm1_out(),
              "smRNA-seq", 
              "hESC", 'DE',
              paste0(input$moi, ' (nc)'))
  })
  
  
  ## scatter plot
  output$scatter_1 <- renderPlot({
    scatter.plot(tf1,
                 input$goi
    )
  })
  
  
  ## DESeq2_summary  
  # smRNA-seq
  miR_summary1 <- reactive({
    
    summary(miR_summary, input$moi, 2)
    
  })
  
  output$miR_summary1 <- renderText({
    miR_summary1()
  })
  
  # RNA-seq
  output$RNA_summary1 <- renderText({
    summary(RNA_summary, input$goi, 2)
  })
  
  # ChRO-seq
  output$ChRO_summary1 <- renderText({
    summary(ChRO_summary, input$goi, 2)
  })
  
  # two factor analysis
  tf_summary1 <- reactive({
    
    tf.summary(tf1, input$goi)
    
  })
  
  output$tf_summary1 <- renderText({
    tf_summary1()
  })
  
  
  ## miR:target relationship
  output$axis1 <- renderText({
    prediction(mirtarget(), miR_summary1(), tf_summary1())
  })
  
  ###############################################
  #
  #    SI lineage formation
  #
  ###############################################
  ## ChRO-seq 
  tpm2_out <- reactive({
    subsetGene_DF.f(tpm2, input$goi, "Duo", "DE")
  })
  
  output$expression_ChRO2 <- renderPlot({
    bar.f(tpm2_out(),
          "ChRO-seq", 
          bquote(italic(paste(.(input$goi))) ~ paste(' (TPM)')))
  })
  
  ## RNA-seq 
  nc2_out <- reactive({
    subsetGene_DF.f(nc2, input$goi, "Duo", "DE")
  })
  
  output$expression_RNA2 <- renderPlot({
    bar.f(nc2_out(),
          "RNA-seq", 
          bquote(italic(paste(.(input$goi))) ~ paste(' (nc)')))
  })
  
  
  ## smRNA-seq 
  rpmmm2_out <- reactive({
    subsetmiR_DF.f(rpmmm, input$moi)
  })
  
  output$expression_miRNA2 <- renderPlot({
    bar.mir.f(rpmmm2_out(),
              "smRNA-seq", 
              "DE", 'Duo',
              paste0(input$moi, ' (nc)'))
  })
  
  
  ## scatter plot
  output$scatter_2 <- renderPlot({
    scatter.plot(tf2,
                 input$goi
    )
  })
  
  
  ## DESeq2_summary  
  # smRNA-seq
  miR_summary2 <- reactive({
    
    summary(miR_summary, input$moi, 3)
    
  })
  
  output$miR_summary2 <- renderText({
    miR_summary2()
  })
  
  # RNA-seq
  output$RNA_summary2 <- renderText({
    summary(RNA_summary, input$goi, 3)
  })
  
  # ChRO-seq
  output$ChRO_summary2 <- renderText({
    summary(ChRO_summary, input$goi, 3)
  })
  
  # two factor analysis
  tf_summary2 <- reactive({
    
    tf.summary(tf2, input$goi)
    
  })
  
  output$tf_summary2 <- renderText({
    tf_summary2()
  })
  
  
  ## miR:target relationship
  output$axis2 <- renderText({
    prediction(mirtarget(), miR_summary2(), tf_summary2())
  })
  
  ###############################################
  #
  #    Regional patterning
  #
  ###############################################
  ## ChRO-seq 
  tpm3_out <- reactive({
    subsetGene_DF.f(tpm3, input$goi, "Ile", "Duo")
  })
  
  output$expression_ChRO3 <- renderPlot({
    bar.f(tpm3_out(),
          "ChRO-seq", 
          bquote(italic(paste(.(input$goi))) ~ paste(' (TPM)')))
  })
  
  ## RNA-seq 
  nc3_out <- reactive({
    subsetGene_DF.f(nc3, input$goi, "Ile", "Duo")
  })
  
  output$expression_RNA3 <- renderPlot({
    bar.f(nc3_out(),
          "RNA-seq", 
          bquote(italic(paste(.(input$goi))) ~ paste(' (nc)')))
  })
  
  
  ## smRNA-seq 
  rpmmm3_out <- reactive({
    subsetmiR_DF.f(rpmmm, input$moi)
  })
  
  output$expression_miRNA3 <- renderPlot({
    bar.mir.f(rpmmm3_out(),
              "smRNA-seq", 
              "Duo", 'Ile',
              paste0(input$moi, ' (nc)'))
  })
  
  
  ## scatter 
  output$scatter_3 <- renderPlot({
    scatter.plot(tf3,
                 input$goi
    )
  })
  
  
  ## DESeq2_summary 
  # smRNA-seq
  miR_summary3 <- reactive({
    
    summary(miR_summary, input$moi, 4)
    
  })
  
  output$miR_summary3 <- renderText({
    miR_summary3()
  })
  
  # RNA-seq
  output$RNA_summary3 <- renderText({
    summary(RNA_summary, input$goi, 4)
  })
  
  # ChRO-seq
  output$ChRO_summary3 <- renderText({
    summary(ChRO_summary, input$goi, 4)
  })
  
  # two factor analysis
  tf_summary3 <- reactive({
    
    tf.summary(tf3, input$goi)
    
  })
  
  output$tf_summary3 <- renderText({
    tf_summary3()
  })
  
  ## miR:target relationship
  output$axis3 <- renderText({
    prediction(mirtarget(), miR_summary3(), tf_summary3())
  })
  
  
}

# Create Shiny object
shinyApp(ui = UI, server = shinyServer)
