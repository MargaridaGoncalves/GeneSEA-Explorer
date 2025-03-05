library(ggplot2)
library(shiny)
library(DT)
library(plotly)
library(gridExtra)
library(dplyr)
library(DESeq2)
library(readxl)
library(readr)
library(bslib)
library(tidyverse)
library(waiter)
library(utils)
library(DOSE)
library(edgeR)
library(sva)
library(affydata)
library(inops)
library(densityClust)
library(wesanderson)
library(paletteer)
library(Glimma)
library(shinyBS)
library(limma)

## Gene Ontology packages

library(BiocGenerics)
library(BiocManager)


library(clusterProfiler)
library(AnnotationDbi)
library(gprofiler2)
library(clustifyr)


library(org.Hs.eg.db) # Human
library(org.Mm.eg.db) # Mouse
library(org.Rn.eg.db) # Rat
library(org.Mmu.eg.db) # Monkey



library(preprocessCore)
library(DEGreport)
library(enrichplot)
library(GOSemSim)
library(pathview)
library(topGO)
library(ggupset)

library(vidger)
library(R.utils)
library(Biobase)


## shiny packages

library(shinyWidgets)

library(shinythemes)

## Shannon PLots

library(GGally)
library(ggbump)
library(genekitr)

# Venn

library(ggVennDiagram)
library(ggvenn)
library(RColorBrewer)


# fonts

library(extrafont)
loadfonts(device = "win", quiet = TRUE) 



#### Ui ####

options(shiny.maxRequestSize=100*1024^2)  # the default maximum file size in shiny is 5 Mb, so this code line allows to change it to 30 Mb


link_shiny <- tags$a(shiny::icon("github"), "Shiny", href = "https://github.com/MargaridaGoncalves", target = "_blank")
#link_posit <- tags$a(shiny::icon("r-project"), "Posit", href = "https://posit.co", target = "_blank")


# Define UI for application that draws a histogram
ui <- navbarPage(title = "GeneSEA Explorer",theme = shinytheme("united"),inverse=TRUE,
                 
                 
                 #tags$style(".row{height: 700px;} .row div:nth-child(1){height: 200%;}"),
                 
                 tabPanel(title = "Home",
                          
                          imageOutput("home_img", height = "1%"),
                          br(),
                          hr(),
                          h4(strong("Project Description")),
                          p(style="text-align: justify; font-size = 25px",
                            "Background: RNA sequencing (RNA-Seq) has become the go-to method for differential gene expression (DGE) analyses in transcriptome research. 
                            Despite the widespread use of R packages like DESeq2 and edgeR for RNA-Seq data analysis, challenges remain in identifying differentially expressed genes (DEG), particularly in selecting suitable normalization methods. 
                            The complexity of these packages often requires programming proficiency, leading researchers to default normalization methods or avoid comparative analyses of DEG results across different techniques.", br(), hr(),
                            
                            "Methods: We introduce a novel web application,", em("GeneSEA Explorer"), ", designed to perform differential gene expression analyses using various normalization methods. 
The application presents outputs through interactive plots and tables, and uses Shannon entropy, a novel approach in the transcriptomics field, to aggregate DEG results, providing statistically supported outcomes.", br(),hr(),

"Results: The ", em("GeneSEA Explorer"), " allows researchers to explore and compare diverse DEGs outcomes across various normalization methods, including less commonly used ones.
The innovative use of Shannon entropy to aggregate all DEG outputs provides an informative selection of DEGs, enhancing researchers’ understanding of their RNA-Seq data results.", br(),hr(),

"Conclusions: The ", em("GeneSEA Explorer"), " is an innovative bioinformatics tool for conducting DGE analyses. Its user-friendly interface enables users to effortlessly explore and analyse diverse DEG outputs. 
The proposed aggregation method, Shannon Entropy Analysis (SEA), aims to minimize challenges for researchers less confident in this domain or those seeking to optimize their time when exploring their data for the first time."
                          ),
br(),
p("Go to",
  a(href = "https://github.com/MargaridaGoncalves",
    "GeneSEA Explorer GitHub Page"),
  "to find more details on the source code."),
p(strong("GeneSEA Explorer is still under continuous development. 
           Please look forward to future updates!")),
hr(),
hr(),
h4(strong("Our Team")),
p(style="text-align: justify; font-size = 25px",
  "Ana Gonçalves", a(href = "https://orcid.org/0009-0001-0800-0019", "https://orcid.org/0009-0001-0800-0019,"),
  "Pedro Macedo", a(href = "https://orcid.org/0000-0002-4371-8069", "https://orcid.org/0000-0002-4371-8069,"),
  "Patrício Costa", a(href = "https://orcid.org/0000-0002-1201-9177", "https://orcid.org/0000-0002-1201-9177,"),
  "Nuno Osório",a(href = "https://orcid.org/0000-0003-0949-5399", "https://orcid.org/0000-0003-0949-5399.")),
hr(),
h4(strong("Acknowledgements")),
hr(),
p(style="text-align: justify; font-size = 25px",
  "This work has been supported by National funds, through the Portuguese Foundation for Science and Technology (FCT) within projects UIDB/50026/2020",
  a(href = "https://doi.org/10.54499/UIDB/50026/2020","(https://doi.org/10.54499/UIDB/50026/2020),"), 
  "UIDP/50026/2020", a(href = "https://doi.org/10.54499/UIDP/50026/2020","(https://doi.org/10.54499/UIDP/50026/2020),"), 
  "LA/P/0050/2020", a(href = "https://doi.org/10.54499/LA/P/0050/2020","(https://doi.org/10.54499/LA/P/0050/2020)"),
  "and 10.54499/CEECINST/00018/2021/CP2806/CT0011", a(href = "https://doi.org/10.54499/CEECINST/00018/2021/CP2806/CT0011","(https://doi.org/10.54499/CEECINST/00018/2021/CP2806/CT0011);"),
  "UIDB/04106/2020", a(href = "https://doi.org/10.54499/UIDB/04106/2020","(https://doi.org/10.54499/UIDB/04106/2020)"),
  "and UIDP/04106/2020", a(href = "https://doi.org/10.54499/UIDP/04106/2020","(https://doi.org/10.54499/UIDP/04106/2020)."),
  hr(),
  "A special thank you to Catarina Lourenço who designed our logo. Go to", a(href = "https://www.behance.net/catariiinalourenco", "Catarina Lourenço's behance account"), "for more informations."),

hr()),




tabPanel(title = "Data Input",
         
         fluidPage(
           
           sidebarLayout(
             sidebarPanel(width = 3,
                          
                          p("Upload txt, xlsx, tsv and csv data files"),
                          br(),
                          fileInput("metadata", label="Input the RNA-Seq count matrix",
                                    accept = c(
                                      'text/csv', ".xlsx",
                                      'text/comma-separated-values',
                                      '.csv', 'text/plain', ".xls"
                                    )),
                          
                          
                          fileInput("coldata", label=" Input the Metadata matrix",
                                    accept = c(
                                      'text/csv', ".xlsx",
                                      'text/comma-separated-values',
                                      '.csv', 'text/plain', ".xls"
                                    )),
                          br(),
                          actionButton("gogo", "Demo dataset"),
                          br(),
                          br(),
                          uiOutput("picker"),
                          br(),
                          uiOutput("Cont"),
                          br(),
                          sliderTextInput(
                            inputId = "Replicates",
                            label = "Number of replicates:", 
                            choices = seq(from = 2,
                                          to = 20,
                                          by = 1),
                            selected = 3,
                            grid = TRUE
                          ),
                          br(),
                          uiOutput("SelectContrast"),
                          textOutput("selected_contr"),
                          br()),
             
             mainPanel(width = 9,
                       
                       DTOutput("Metadatatable1"),
                       DTOutput("Infotable"))
             
             
             
           ))),




tabPanel(title = "Differential Gene Expression Analysis",
         tags$style(".row{height: 700px;} .row div:nth-child(1){height: 200%;}"),
         waiter::autoWaiter(html = spin_loaders(id = 34, color = "#00688B"), color = "white"),
         fluidPage(
           
           
           tabsetPanel(
             id = "t",
             type = "tabs",
             
             
             tabPanel("Volcano Plots", 
                      fluidPage(
                        
                        
                        tags$div(class="clearfix"),
                        fluidRow(
                          
                          tags$div(style = "margin-top:-97em"),
                          column(4, p(strong("No Normalization")),
                                 h5("Number of differentially expressed genes:"),
                                 textOutput("DEnumberNN"),
                                 plotOutput("VolcanoNN")),
                          column(4, 
                                 p(strong("TMM")), 
                                 h5("Number of differentially expressed genes:"),
                                 textOutput("DEnumberTMM"),
                                 plotOutput("VolcanoTMM")),
                          column(4, 
                                 p(strong("TMMwsp")),
                                 h5("Number of differentially expressed genes:"),
                                 textOutput("DEnumberTMMwsp"),
                                 plotOutput("VolcanoTMMwsp")),
                          
                          
                          tags$div(style = "margin-top:-20em"),
                          
                          column(4, 
                                 br(),
                                 p(strong("Deseq2")),
                                 h5("Number of differentially expressed genes:"),
                                 textOutput("DEnumberDeseq2"),
                                 plotOutput("VolcanoDeseq2")),
                          column(4, 
                                 br(),
                                 p(strong("RLE")),
                                 h5("Number of differentially expressed genes:"),
                                 textOutput("DEnumberRLE"),
                                 plotOutput("VolcanoRLE")),
                          column(4, 
                                 br(),
                                 p(strong("PoissonSeq")),
                                 h5("Number of differentially expressed genes:"),
                                 textOutput("DEnumberPoissonSeq"),
                                 plotOutput("VolcanoPoisson"))),
                        
                        
                        tags$div(style = "margin-top:-20em"),
                        
                        
                        column(4, 
                               br(),
                               p(strong("Median")), 
                               h5("Number of differentially expressed genes:"),
                               textOutput("DEnumberMedian"),
                               plotOutput("VolcanoMedian")),
                        column(4, 
                               br(),
                               p(strong("Upper-Quartile")),
                               h5("Number of differentially expressed genes:"),
                               textOutput("DEnumber_UpperQuartile"),
                               plotOutput("VolcanoUpper")),
                        column(4,
                               br(),
                               p(strong("Quantile")),  
                               h5("Number of differentially expressed genes:"),
                               textOutput("DEnumberQuantile"),
                               plotOutput("VolcanoQuantile")),
                        
                        tags$div(style = "margin-top:-20em"),
                        
                        
                        column(4, 
                               br()),
                        column(4, 
                               br()),
                        column(4,
                               br()),
                        
                        tags$div(style = "margin-top:-20em"),
                        
                        
                        column(4, 
                               br(),
                               h6("Select the desired adjusted p-value and log fold change thresholds:"),
                               sliderInput(
                                 inputId = "padj_variable",
                                 label = "adjusted p-value",
                                 min = 0.01,
                                 max = 0.1,
                                 value = 0.05,
                                 step = 0.001
                               ),
                               
                               sliderInput(
                                 inputId = "lfc_variable",
                                 label = "log2 Fold Change",
                                 min = 0,
                                 max = 5,
                                 value = 2,
                                 step = 0.1)),
                        column(4, 
                               br()),
                        column(4,
                               br())
                        
                      )
                      
                      
                      #### Ui results ####                 
             ),
             
             tabPanel("Results", 
                      
                      br(),
                      h6("Select the desired adjusted p-value and log fold change thresholds and the normalization method"),
                      h6("The values selected in this section will be employed in the subsequent Over-representation Analysis (ORA)."),
                      sliderInput(
                        inputId = "FDR",
                        label = "adjusted p-value",
                        min = 0.01,
                        max = 0.1,
                        value = 0.05,
                        step = 0.001
                      ),
                      
                      sliderInput(
                        inputId = "LFC",
                        label = "log2 Fold Change",
                        min = 0,
                        max = 5,
                        value = 2,
                        step = 0.1),
                      uiOutput("Normalization"),
                      br(),
                      br(),
                      p("The following plots show the multi-dimensional scaling (MDS) results."),
                      plotOutput("DataMDS1"),
                      br(),
                      actionButton("glimma", "MDS with Glimma"),
                      br(),
                      br(),
                      p("The following plot shows the genes that are differentially expressed in the current dataset via Volcano Plot."),
                      plotlyOutput("Volcano"),
                      br(),
                      p("The following table shows the genes that are differentially expressed in the current dataset."),
                      br(),
                      DT::dataTableOutput("Tables"),
                      br(),
                      br(),
                      p("The following table depicts how many differentially expressed genes are up-regulated or down-regulated for each normalization method and how many genes from the DEG list differ between the current normalization method and the others."),
                      textOutput("ResultsNameNorm"),
                      br(),
                      DT::dataTableOutput("DataUpDown"),
                      br(),
                      h5("Venn Diagram"),
                      h6("The following output exhibits the results from the DEGs lists"),
                      br(),
                      uiOutput("VennOptions1"),
                      br(),
                      plotOutput("Venn", width = "500px", height = "500px"),
                      br(),
                      plotOutput("DataGlimmaMDS", width = "1px", height = "1px"),
                      uiOutput("SelectGene"),
                      plotOutput("PlotlyCount")
                      
                      
                      #### Ui Differential expression analysis ####           
                      
             ),
             
             tabPanel("Differential Gene Expression Analysis with SEA", 
                      br(),
                      p("In this section the workflow of the Shannon Entropy Analysis (SEA) is exhibited."),
                      h6("It is crucial to note that this section uses the complete list of genes, rather than solely the differentially expressed genes (DEGs) list."),
                      br(),
                      p("Select at least 4 Normalization Methods to perform the SEA:"),
                      uiOutput("MatrizE"),
                      br(),
                      p("The following table displays the inicial E matrix that will be used to calculate the Shannon entropy weigths."),
                      h6("The values inputed in the matrix are the adjusted p-values obtained for each gene and normalization method."),
                      br(),
                      DT::dataTableOutput("MEE"), 
                      br(),
                      p("The following table shows the degree of importance (Wl) of each normalization method."),
                      br(),
                      DT::dataTableOutput("weigths"), 
                      br(),
                      p("The table shows the ranking of each gene."),
                      br(),
                      DT::dataTableOutput("ShannonRT"),
                      br(),
                      p("The following table allows to compare the ranks of each gene between methods."),
                      br(),
                      DT::dataTableOutput("FinalShannonRT"),
                      br(),
                      p("The following table shows the genes ranking and the corresponding score."),
                      h6("The score is determined by the frequency at which a specific gene appears among the top n genes across multiple gene lists derived from the selected normalization methods."),
                      br(),
                      sliderInput("integer", "Select the number of top ranked genes to be observed",
                                  min = 0, max = 100,
                                  value = 20),
                      DT::dataTableOutput("tenta"),
                      br(),
                      h5("Venn Diagram"),
                      h6("Comparison between the SEA and various normalization DEGs lists."),
                      br(),
                      uiOutput("VennOptions"),
                      br(),
                      plotOutput("Venn2", width = "500px", height = "500px"),
                      br())
             
             
           ))
         
),

### Ui Erichment Analysis ####

tabPanel(title = "Functional Enrichment Analysis",
         
         fluidPage(
           
           sidebarLayout(
             sidebarPanel(width = 3,
                          br(),
                          h6("Gene Ontology OPTIONS"),
                          uiOutput("Ontology"),
                          uiOutput("NormOnto"),
                          uiOutput("OrgDb"),
                          uiOutput("KeyType"),
                          
                          sliderInput(
                            inputId = "ShowCategory",
                            label = "Number of Categories to be shown",
                            min = 5,
                            max = 25,
                            value = 10,
                            step = 1),
                          br(),
                          h6("Kyoto Encyclopedia of Genes and Genomes OPTIONS"),
                          
                          uiOutput("KEGGkeytype"),
                          uiOutput("NormKegg"),
                          uiOutput("KEGGorganism"),
                          
                          sliderInput(
                            inputId = "ShowKEGGCategory",
                            label = "Number of Categories to be shown",
                            min = 1,
                            max = 25,
                            value = 5,
                            step = 1),
                          h6("To elaborate the KEEG enrichment analysis, the GeneOntology", em("OrgDB terms"), "and", em("keyType"), "options must also be selected."),
                          
                          
             ),
             
             mainPanel(width = 9,
                       
                       
                       tabPanel("Plot outputs", tabsetPanel(
                         
                         
                         tabPanel("Barplot", tabsetPanel(
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12,                         
                                               br(),
                                               textOutput("NameNorm13"),
                                               br(),
                                               plotOutput("fit"),
                                               br()),
                                        column(12,  
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("fitShannon"),
                                               br())))),
                           
                           
                           
                           tabPanel("KEGG", 
                                    
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm14"),
                                               br(),
                                               plotOutput("KEGGbarplot"),
                                               br()),
                                        column(12,           
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGbarplot_Shannon"),
                                               br()
                                        ))))
                           
                         )),
                         
                         
                         
                         tabPanel("Dotplot", tabsetPanel(
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm11"),
                                               br(),
                                               plotOutput("dotplot"),
                                               br()),
                                        column(12,  
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("dotplotShannon"),
                                               br(),
                                               br(),
                                               textOutput("NameNorm17"),
                                               br(),
                                               plotOutput("dotplotGSEA"),
                                               br())))),
                           
                           
                           tabPanel("KEGG", 
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm12"),
                                               br(),
                                               plotOutput("KEGGdotplot"),
                                               br()),
                                        column(12, 
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGdotplot_Shannon"),
                                               br()))))
                           
                         )),
                         
                         
                         
                         tabPanel("Goplot",tabsetPanel(
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12,   
                                               br(),
                                               textOutput("NameNorm8"),
                                               br(),
                                               plotOutput("goplot"),
                                               br()),
                                        column(12, 
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("goplotShannon"),
                                               br()))))
                           
                           # goplot com KEGG does not work
                           
                         )),
                         
                         
                         
                         tabPanel("Emapplot", tabsetPanel(
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12,  
                                               br(),
                                               textOutput("NameNorm5"),
                                               br(),
                                               plotOutput("emapplot"),
                                               br()),
                                        
                                        column(12,  
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("emapplotShannon"),
                                               br(),
                                               br(),
                                               textOutput("NameNorm19"),
                                               br(),
                                               plotOutput("emapplotGSEA"),
                                               br())))),
                           
                           
                           
                           tabPanel("KEGG", 
                                    
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm16"),
                                               br(),
                                               plotOutput("KEGGemapplot"),
                                               br()),
                                        column(12,           
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGemapplot_Shannon"),
                                               br()))))
                           
                           
                         )),
                         
                         
                         
                         tabPanel("cneplot", tabsetPanel(
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm"),
                                               br(),
                                               plotOutput("cnetplot"),
                                               br(),
                                               br(),
                                               plotOutput("cnetplot2"),
                                               br()),
                                        column(12,  
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("cnetplotShannon"),
                                               br(),
                                               br(),
                                               plotOutput("cnetplot2Shannon"),
                                               br(),
                                               br(),
                                               textOutput("NameNorm20"),
                                               br(),
                                               plotOutput("cnetplotGSEA"),
                                               br(),
                                               plotOutput("cnetplot2GSEA"),
                                               br())))),
                           tabPanel("KEGG", 
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm1"),
                                               br(),
                                               plotOutput("KEGGcnetplot"),
                                               br(),
                                               br(),
                                               plotOutput("KEGGcnetplot21"),
                                               br()),
                                        column(12, 
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGcnetplot_Shannon"),
                                               br(),
                                               br(),
                                               plotOutput("KEGGcnetplot21_Shannon"),
                                               br()))))
                           
                           
                           
                           
                         )),
                         
                         
                         tabPanel("Heatplot", tabsetPanel(
                           
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12,                         
                                               br(),
                                               textOutput("NameNorm3"),
                                               br(),
                                               plotOutput("heatplot2"),
                                               br()),
                                        column(12,  
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("heatplot2Shannon"),
                                               br(),
                                               br(),
                                               textOutput("NameNorm21"),
                                               br(),
                                               plotOutput("heatplot2GSEA"),
                                               br())))),
                           
                           
                           
                           tabPanel("KEGG", 
                                    
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm4"),
                                               br(),
                                               plotOutput("KEGGheatplot"),
                                               br()),
                                        column(12,           
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGheatplot_Shannon"),
                                               br()))))
                           
                         )),
                         
                         
                         tabPanel("Treeplot", tabsetPanel(
                           tabPanel("Gene Ontology",
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12,   
                                               br(),
                                               textOutput("NameNorm7"),
                                               br(),
                                               plotOutput("treeplot"),
                                               br()),
                                        
                                        column(12, 
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("treeplotShannon"),
                                               br(),
                                               br(),
                                               textOutput("NameNorm22"),
                                               br(),
                                               plotOutput("treeplotGSEA"),
                                               br())))),
                           
                           tabPanel("KEGG", 
                                    
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm15"),
                                               br(),
                                               plotOutput("KEGGtreeplot"),
                                               br()),
                                        column(12,           
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGtreeplot_Shannon"),
                                               br()))))
                           
                           
                           
                           
                           
                         )),
                         
                         
                         tabPanel("upsetplot",tabsetPanel(
                           tabPanel("Gene Ontology",
                                    
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12,                         
                                               br(),
                                               textOutput("NameNorm10"),
                                               br(),
                                               plotOutput("upsetplot"),
                                               br()),
                                        column(12,  
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("upsetplotShannon"),
                                               br(),
                                               br(),
                                               textOutput("NameNorm23"),
                                               br(),
                                               plotOutput("upsetplotGSEA"),
                                               br())))),
                           
                           tabPanel("KEGG",  
                                    
                                    fluidPage(
                                      tags$div(class="clearfix"),
                                      fluidRow(
                                        
                                        tags$div(style = "margin-top:-97em"),
                                        column(12, 
                                               br(),
                                               textOutput("NameNorm9"),
                                               br(),
                                               plotOutput("KEGGupsetplot"),
                                               br()),
                                        column(12,           
                                               br(),
                                               h5("ORA with SEA:"),
                                               br(),
                                               plotOutput("KEGGupsetplot_Shannon"),
                                               br()))))
                           
                           
                           
                         ))
                         
                         
                       )))))),



navbarMenu(title = "Links",
           nav_item(link_shiny)),

)





# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  ##### Metadata call ##############
  
  
  randomVals <- eventReactive(input$gogo, {
    read.table("exp.txt", header = TRUE)
    
  })
  
  
  Metadata <- reactive({
    
    
    infile <- input$metadata
    
    if (is.null(infile)) {
      
      
      return(randomVals())
      
      
    } else if (stringr::str_ends(infile$datapath, "(xlsx|xls)")) {
      
      return(read_xlsx(infile$datapath))
      
      
    } else if (stringr::str_ends(infile$datapath, "csv")) {
      
      return(read.csv(infile$datapath, header = TRUE))
      
      
    } else {
      return(read.table(infile$datapath, header = TRUE))
      
    }
    
    
  })
  
  
  
  output$Metadatatable <- renderDT({
    df <- Metadata()
    df
  })
  
  
  ##### coldata call ##############
  
  randomVals2 <- eventReactive(input$gogo, {
    read.table("meta.txt", header = TRUE)
    
  })
  
  
  
  Informationtable <- reactive({
    
    
    infile <- input$coldata
    
    
    if (is.null(infile)) {
      
      
      return(randomVals2())
      
      
    } else if (stringr::str_ends(infile$datapath, "(xlsx|xls)")) {
      
      read_xlsx(infile$datapath)
      
    } else if (stringr::str_ends(infile$datapath, "csv")) {
      
      read.csv(infile$datapath, header = TRUE)
      
    } else {
      read.table(infile$datapath, header = TRUE)
    }
    
    
  })
  
  
  output$Infotable<- renderDT({
    df1 <- Informationtable()
    df1
  })
  
  
  
  #### Design variable #####
  
  
  
  
  output$picker <-  renderUI ({
    
    infile <- Informationtable()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- pickerInput(inputId = 'pick', 
                          label = 'Variable to be studied', 
                          choices = colnames(Informationtable()),
                          selected = colnames(Informationtable()[2]), # to select which variable to start with
                          options = list(`actions-box` = TRUE),multiple = F)
      
      data
    }
    
    return(data)
    
  })
  
  
  group <- reactive({ 
    
    infile <- Informationtable()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- Informationtable1()[,input$pick]
      
    }
  })
  
  
  
  
  
  design <- reactive ({
    
    infile <- group()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- model.matrix(~ 0 + group()) 
      colnames(data) <- levels(group())
    }
    return(data)
    
  })
  
  
  
  
  
  #### Design Contrasts ####
  
  contrastList <- reactive({
    
    infile <- Informationtable1()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- as.factor(Informationtable1()[,input$pick])
    }
    return(data)
  })
  
  
  group2 <- reactive({ 
    
    infile <- contrastList()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- levels(contrastList()[input$pick])
    }
    return(data)
  })
  
  
  output$Cont <-  renderUI ({
    data <- pickerInput(inputId = 'cont', 
                        label = 'Variables to be compared', 
                        choices = group2(),
                        selected = group2()[1:2],
                        options = list(`actions-box` = FALSE,  "max-options" = 2), multiple = T)
    
    req(data)
    data
    
  })
  
  
  vetorContraste1 <- reactive({
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      data <- group2() %in{}% c(input$cont) 
      data <- as.integer(as.logical(data))
    }
    return(data)
  })
  
  
  vetorContraste <- reactive({
    data <- vetorContraste1()
    index <- which(vetorContraste1() == 1)[1]
    data[index] <- -1
    data
  })
  
  
  othervetorContraste <- reactive({
    data <- vetorContraste1()
    index <- which(vetorContraste1() == 1)[2]
    data[index] <- -1
    data
  })
  
  output$contt <- renderText({vetorContraste()})
  output$contt2 <- renderText({ paste(c(c1()), c(c2()), sep =",")})
  
  
  c1 <- reactive(vetorContraste())
  c2 <- reactive(othervetorContraste())
  
  output$SelectContrast <-  renderUI ({
    data <- pickerInput(inputId = 'SelectContr', 
                        label = 'Contrast', 
                        choices = list(paste(c1(), c(rep(c(","), times = (length(c1())-1)),""), collapse = " "), paste(c2(), c(rep(c(","), times = (length(c2())-1)),""), collapse = " ")),
                        options = list(`actions-box` = TRUE), multiple = F)
    
    req(data)
    data
    
  })
  
  
  ## Contrast phrase ##
  
  output$selected_contr <- renderText({ 
    
    infile <- input$SelectContr
    req(infile)
    
    
    if (is.null(infile)) {
      
      return(NULL)
      
    } else {
      num_vector <- as.numeric(strsplit(input$SelectContr, ",")[[1]]) # Gives me the numerical contrast
      data1 <- input$cont# Gives me the 2 variable names 
      name_vector <- group2() # Gives me all the variable names
      data3 <- group2() %in% input$cont # Gives me a True False vector
      
      # Initialize variables to store the names at 1 and -1 spots
      name_at_1 <- NULL
      name_at_neg1 <- NULL
      
      # Loop through the numerical vector to find the positions of 1 and -1
      for (i in 1:length(num_vector)) {
        if (num_vector[i] == 1) {
          name_at_1 <- name_vector[i]
        } else if (num_vector[i] == -1) {
          name_at_neg1 <- name_vector[i]
        }
      }
      
      f <- paste(name_at_1, "vs", name_at_neg1) 
    }
    return(f)
    
  })
  
  
  ###
  
  
  
  Metadata1 <- reactive({                   
    data <- Metadata()
    row.names(data) <- Metadata()[,1]
    data <- data[,-1]
    data <- as.matrix(data)
    data
  })
  
  
  
  output$Metadatatable1<- DT::renderDataTable({
    df1 <- Metadata1()
    DT::datatable(df1)
  })
  
  
  
  
  
  Informationtable1 <- reactive({    
    infile <- Informationtable()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      data <- as.data.frame(Informationtable())
      data <- mutate(data,  across(2,as.factor))
      data
    }
    return(data)
  })
  
  
  
  
  
  #### DESeq Analysis #####
  
  
  dds <- reactive({ 
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- DESeqDataSetFromMatrix(countData = Metadata1(),
                                     colData = Informationtable1(),
                                     design= ~ 0 + group()) 
      data
    }
    return(data)
  })
  
  
  
  cutoffvalueDeseq2 <- reactive ({ 
    
    L <-  as.numeric(min(colSums(Metadata1()))) 
    threshold <- as.numeric((10)/(L/10^6))
    threshold 
  })
  
  
  dds01 <- reactive({                   
    data <- dds()
    row.names(data) <- c(row.names(Metadata1()))
    data
  })
  
  ## MDS  Glimma PLOT ##
  
  deseqGlimmaMDS <- reactive({glMDSPlot(dds(), groups = as.data.frame(SummarizedExperiment::colData(dds())),
                                        labels = rownames(SummarizedExperiment::colData(dds())))}) 
  
  ####
  
  
  
  ### Deseq2 MDS
  
  MDS_Deseq2 = function(dds, group){
    
    
    Keep <- rowSums(cpm(counts(dds)) > cutoffvalueDeseq2() ) >= input$Replicates # number of replicates
    dds2 <- dds[Keep,]
    
    dds1 <- DESeq(dds2)
    res <- results(dds1, contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]]))
    
    top_genes <- reactive({head(order(res$padj, na.last=NA), length(row.names(dds1))) }) 
    
    dds_norm <- reactive({
      #dds_obj <- vst(dds2, blind = FALSE)
      #dds_obj <- rlog(dds, blind = FALSE)
      dds_obj <- normTransform(dds1)
      dds_obj <- assay(dds_obj)[top_genes(), ]
      dds_obj
    })
    
    
    y <- reactive ({ dds_norm()})
    
    
    col.group <- reactive ({ 
      data <- group
      data <- factor(group)
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      
      limma::plotMDS(y(), top = 100, labels = group, col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 100, labels = rownames(SummarizedExperiment::colData(dds())), cex = 1.4, col= col.group())
      title(main="Samples")
      
    })
    
    
    #plot2 <- reactive({  
    #  limma::plotMDS(y(), top = 500, col= col.group())
    #  title(main="Sample groups")
    #})
    
    return(plot1())
    
  }
  
  results_MDSDeseq2 <- reactive ({ MDS_Deseq2(dds(), group()) })
  
  
  
  
  
  ####
  
  
  dds1 <- reactive({ DESeq(dds01()) })
  
  
  output$ress <- renderText({ input$pick }) #resultsNames(dds3())
  
  
  
  Keep <- reactive({ rowSums(cpm(counts(dds1())) > cutoffvalueDeseq2()) >= input$Replicates })
  dds2 <- reactive({ dds1()[Keep(),] })
  
  
  dds3 <- reactive({ DESeq(dds2())})
  
  
  res7 <- reactive({results(dds3(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) })
  
  res8 <- reactive ({as.data.frame(res7())})
  res9 <- reactive({ res8()[complete.cases(res8()), ] })
  
  res10 <-reactive({ res9() %>% 
      mutate(regulated = case_when(
        padj < input$padj_variable & log2FoldChange < -input$lfc_variable ~ "down", # significantly down
        padj < input$padj_variable & log2FoldChange > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  res11 <- reactive({    
    
    data <- as.data.frame(res10())
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(res10()))
    data
  })
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resdeseq10 <- reactive({ res9() %>% 
      mutate(regulated = case_when(
        padj< input$FDR & log2FoldChange < -input$LFC ~ "down", # significantly down
        padj < input$FDR & log2FoldChange > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resdeseq11 <- reactive({    
    
    data <- resdeseq10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resdeseq10()))
    data$Regulated <- as.factor(data[,"regulated"])
    data
  })
  
  
  datDeseq<- reactive({
    data <- res9()
    data <- as.data.frame(res9())
    data$regulated <- with(data, ifelse(padj< input$FDR & log2FoldChange < -input$LFC, "down-regulated", ifelse(padj < input$FDR & log2FoldChange > input$LFC, "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$padj < input$FDR)
    data <- data %>% filter(abs(data$log2FoldChange) > input$LFC)
    data
  })
  
  Deseq2v2 <- reactive({                   
    data <- datDeseq()
    data
  })
  
  Deseq2v3 <- reactive({
    data <- res9()
    data <- as.data.frame(res9())
    data$regulated <- with(data, ifelse(padj< input$FDR & log2FoldChange < -input$LFC, "down-regulated", ifelse(padj < input$FDR & log2FoldChange > input$LFC, "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  ### Data for Volcano Home Page ###
  
  
  resdeseqa <- reactive({ res9() %>% 
      mutate(regulated = case_when(
        padj< input$padj_variable & log2FoldChange < -input$lfc_variable ~ "down", # significantly down
        padj < input$padj_variable & log2FoldChange > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resdeseqb <- reactive({    
    
    data <- resdeseqa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resdeseqa()))
    data$Regulated <- as.factor(data[,"regulated"])
    data
  })
  
  
  
  dataset<- reactive({
    data <- res7()
    data <- as.data.frame(res7())
    data <- data %>% filter(data$padj < input$padj_variable)
    data <- data %>% filter(abs(data$log2FoldChange) > input$lfc_variable)
  })
  
  output$DEnumberDeseq2 <- renderText({
    data <- length(c(row.names(dataset())))
    data
  }) 
  
  ###
  
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_Deseq2 <- reactive({
    data <- c(row.names(datDeseq()))
    data
  })
  
  ###
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownDeseq <- reactive({ 
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpDeseq <- reactive({ length(which(datDeseq()$regulated == "up-regulated")) })
      DownDeseq <- reactive({ length(which(datDeseq()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `DESeq2` = c(UpDeseq(), DownDeseq()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownDeseq <- DT::renderDataTable({
    df2 <- TableUpDownDeseq()
    DT::datatable(df2)
    
    
  })
  
  
  
  GenesUpDownDeseq2 <- reactive({ 
    
    
    ## NN
    equal2Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datNN()))) })
    dif2Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datNN()))))  })
    
    ## Median 
    equal1Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datMedian())))  })
    dif1Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datTC()))) })
    dif6Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datPoissonSeq()))) })
    dif3Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datDeseq()))) })
    dif11Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datRLE()))) })
    dif5Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datTMM()))) })
    dif7Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datTMMwsp()))) })
    dif8Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datSVA()))) })
    dif10Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datUQ()))) })
    dif9Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4Deseq2 <- reactive({ length(which(row.names(datDeseq()) %in% row.names(datQuantile()))) })
    dif4Deseq2 <- reactive({ length(which(!(row.names(datDeseq()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      #`DESeq2` = c(equal11Deseq2(), dif11Deseq2()),
      `Median` = c(equal1Deseq2(), dif1Deseq2()),
      `No Normalization` = c(equal2Deseq2(), dif2Deseq2()),
      `PoissonSeq` = c(equal3Deseq2(), dif3Deseq2()),
      `Quantile` = c(equal4Deseq2(), dif4Deseq2()),
      `RLE` = c(equal5Deseq2(), dif5Deseq2()),
      #`Total Count` = c(equal6Deseq2(), dif6Deseq2()),
      `TMM` = c(equal7Deseq2(), dif7Deseq2()),
      `TMMwsp` = c(equal8Deseq2(), dif8Deseq2()),
      `Upper-Quartile` = c(equal9Deseq2(), dif9Deseq2()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  
  
  
  #
  
  
  
  
  dataset1 <- reactive({                   
    data <- dataset()
    data
  })
  
  
  Deseq2 <- reactive({                   
    data <- dataset()
    data
  })
  
  
  output$DEgenes<- DT::renderDataTable({
    df2 <- dataset()
    DT::datatable(df2)
    
    
  })
  
  
  
  
  output$DE1<- DT::renderDataTable({
    dataset1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  results_deseq1 <- reactive ({
    data <- res7()
    data <- as.data.frame(res7())
    data$rank <- c(1:nrow(res7()))
    data$Method <- c(rep("Deseq2", times = nrow(res7())))
    data$Gene <- c(rownames(res7()))
    data
  })
  
  
  results_deseq11 <- reactive ({
    data <- as.data.frame(results_deseq1())
    data <- results_deseq1()[order(row.names(results_deseq1())),]
    data
  })
  
  
  
  ## Data for the PlotCount ##
  
  # create the list of genes
  
  
  listGene <- reactive({ Metadata1()[order(row.names(Metadata1())),] }) 
  
  
  output$SelectGene <-  renderUI ({
    
    infile <- input$NormalizationMethod
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  pickerInput(inputId = 'SGene', 
                           label = 'Select the Gene of interest: ', 
                           choices = rownames(listGene()),
                           options = list(`actions-box` = TRUE), multiple = F)
      fit3
      
    } else {
      
      fit3 <- NULL
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  PlotCount <- reactive({
    
    infile <- input$NormalizationMethod
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  plotCounts(dds1(), gene = input$SGene , intgroup= input$pick )
      fit3
      
    } else {
      
      fit3 <- NULL
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  output$PlotlyCount <-  renderPlot({
    data <- PlotCount()
    data
  })
  
  
  
  
  
  
  
  
  ### TMM ####
  
  DE_edgeR_TMM = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      cutoffvalue <- reactive ({ 
        
        L <- as.numeric(min(colSums(expression.matrix))) 
        threshold <- as.numeric((10)/(L/10^6))
        threshold 
      })
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue()) >= input$Replicates }) # Use a random column to decide the cutoff threshold
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      # TMM normalise
      expression.matrix.normalised <- reactive({
        data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "TMM")) ## method = "TMM"
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
    }
    return(res2())
    
  }
  
  results_TMM <- reactive ({ DE_edgeR_TMM(Metadata1(),Informationtable1(), design()) })
  
  
  ### For Results Page ###
  
  resTMM10 <- reactive({ results_TMM() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resTMM11 <- reactive({    
    
    data <- resTMM10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resTMM10()))
    data
  })
  
  
  ### 
  
  
  ### For Volcanos Home Page ###
  
  
  resTMMa <- reactive({ results_TMM() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resTMMb <- reactive({    
    
    data <- resTMMa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resTMMa()))
    data
  })
  
  
  
  
  datasetTMM<- reactive({
    data <- results_TMM()
    data <- as.data.frame(results_TMM())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  output$DEnumberTMM <- renderText({
    data <- length(c(row.names(datasetTMM())))
    data
  })  
  
  
  ###
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  results_TMM1 <- reactive ({
    data <- results_TMM()
    data <- as.data.frame( results_TMM())
    data$rank <- c(1:nrow(results_TMM()))
    data$Method <- c(rep("TMM", times = nrow(results_TMM())))
    data$Gene <- c(rownames(results_TMM()))
    data
    
  })
  
  
  
  results_TMM11 <- reactive ({
    
    data <- results_TMM1()[order(results_TMM1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
  })
  
  
  datTMM <- reactive({
    data <- results_TMM()
    data <- as.data.frame(results_TMM())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  TMMv2 <- reactive({                   
    data <- datTMM()
    data
  })
  
  TMMv3 <- reactive({
    data <- results_TMM()
    data <- as.data.frame(results_TMM())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  
  ### Code for Venn Diagram ###
  
  Venn_TMM <- reactive({
    data <- c(row.names(datTMM()))
    data
  })
  
  ###
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownTMM <- reactive({ 
    
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpTMM <- reactive({ length(which(datTMM()$regulated == "up-regulated")) })
      DownTMM <- reactive({ length(which(datTMM()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `TMM` = c(UpTMM(), DownTMM()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownTMM <- DT::renderDataTable({
    df2 <- TableUpDownTMM()
    DT::datatable(df2)
    
    
  })
  
  
  
  
  
  GenesUpDownTMM <- reactive({ 
    
    
    ## NN
    equal2TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datNN()))) })
    dif2TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datNN()))))  })
    
    ## Median 
    equal1TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datMedian())))  })
    dif1TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datTC()))) })
    dif6TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datPoissonSeq()))) })
    dif3TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datDeseq()))) })
    dif11TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datRLE()))) })
    dif5TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datTMM()))) })
    dif7TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datTMMwsp()))) })
    dif8TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datSVA()))) })
    dif10TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datUQ()))) })
    dif9TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4TMM <- reactive({ length(which(row.names(datTMM()) %in% row.names(datQuantile()))) })
    dif4TMM <- reactive({ length(which(!(row.names(datTMM()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11TMM(), dif11TMM()),
      `Median` = c(equal1TMM(), dif1TMM()),
      `No Normalization` = c(equal2TMM(), dif2TMM()),
      `PoissonSeq` = c(equal3TMM(), dif3TMM()),
      `Quantile` = c(equal4TMM(), dif4TMM()),
      `RLE` = c(equal5TMM(), dif5TMM()),
      #`Total Count` = c(equal6TMM(), dif6TMM()),
      #`TMM` = c(equal7TMM(), dif7TMM()),
      `TMMwsp` = c(equal8TMM(), dif8TMM()),
      `Upper-Quartile` = c(equal9TMM(), dif9TMM()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  
  
  
  
  #
  
  
  datasetTMM1 <- reactive({                   
    data <- datasetTMM()
    data
  })
  
  
  TMM <- reactive({                   
    data <- datasetTMM()
    data
  })
  
  
  
  output$DEgenesTMM<- DT::renderDataTable({
    df2 <- datasetTMM()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DETMM1<- DT::renderDataTable({
    datasetTMM1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  ### TMM MDS
  
  MDS_TMM = function(expression.matrix, metadata, design){
    
    cutoffvalue <- reactive ({ 
      
      L <-  as.numeric(min(colSums(expression.matrix))) 
      threshold <- as.numeric((10)/(L/10^6))
      threshold 
    })
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates }) # Use a random column to decide the cutoff threshold
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # TMM normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "TMM")) ## method = "TMM"
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col=  col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    
    return(plot1())
    
  }
  
  results_MDSTMM <- reactive ({ MDS_TMM(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  ### TMM Glimma
  
  MDSGlimma_TMM = function(expression.matrix, metadata, design){
    
    cutoffvalue <- reactive ({ 
      
      L <-  as.numeric(min(colSums(expression.matrix))) 
      threshold <- as.numeric((10)/(L/10^6))
      threshold 
    })
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates }) # Use a random column to decide the cutoff threshold
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # TMM normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "TMM")) ## method = "TMM"
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSTMM <- reactive ({ MDSGlimma_TMM(Metadata1(),Informationtable1(), design()) })
  
  
  DataGlimmaMDS <- renderPlot({
    data <- results_GlimmaMDSTMM()
    data
    
  })
  
  
  
  
  
  ### RLE ####
  
  DE_edgeR_RLE = function(expression.matrix, metadata, design){
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      cutoffvalue <- reactive ({               
        L <-  as.numeric(min(colSums(expression.matrix)))        
        threshold <- as.numeric((10)/(L/10^6))       
        threshold      })          
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      # RLE normalise
      expression.matrix.normalised <- reactive({
        data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "RLE")) ## method = "RLE"
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
    }
    return(res2())
    
  }
  
  results_RLE <- reactive ({ DE_edgeR_RLE(Metadata1(),Informationtable1(), design()) })
  
  
  
  ## DATA FOR THE VOLCANO PLOT  from results ##
  
  resRLE10 <- reactive({ results_RLE() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resRLE11 <- reactive({    
    
    data <- resRLE10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resRLE10()))
    data
  })
  
  
  
  datRLE <- reactive({
    data <- results_RLE()
    data <- as.data.frame(results_RLE())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  RLEv2 <- reactive({                   
    data <- datRLE()
    data
  })
  
  RLEv3 <- reactive({
    data <- results_RLE()
    data <- as.data.frame(results_RLE())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  
  ### For Volcanos Home Page ###
  
  
  resRLEa <- reactive({ results_RLE() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resRLEb <- reactive({   
    data <- resRLEa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resRLEa()))
    data
  })
  
  
  
  
  datasetRLE<- reactive({
    data <- results_RLE()
    data <- as.data.frame(results_RLE())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  output$DEnumberRLE <- renderText({
    data <- length(c(row.names(datasetRLE())))
    data
  })    
  
  
  
  ###
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_RLE <- reactive({
    data <- c(row.names(datRLE()))
    data
  })
  
  ###
  
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownRLE <- reactive({ 
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpRLE <- reactive({ length(which(datRLE()$regulated == "up-regulated")) })
      DownRLE <- reactive({ length(which(datRLE()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `RLE` = c(UpRLE(), DownRLE()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownRLE <- DT::renderDataTable({
    df2 <- TableUpDownRLE()
    DT::datatable(df2)
    
    
  })
  
  
  
  GenesUpDownRLE <- reactive({ 
    
    
    ## NN
    equal2RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datNN()))) })
    dif2RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datNN()))))  })
    
    ## Median 
    equal1RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datMedian())))  })
    dif1RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datTC()))) })
    dif6RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datPoissonSeq()))) })
    dif3RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datDeseq()))) })
    dif11RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datRLE()))) })
    dif5RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datTMM()))) })
    dif7RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datTMMwsp()))) })
    dif8RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datSVA()))) })
    dif10RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datUQ()))) })
    dif9RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4RLE <- reactive({ length(which(row.names(datRLE()) %in% row.names(datQuantile()))) })
    dif4RLE <- reactive({ length(which(!(row.names(datRLE()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11RLE(), dif11RLE()),
      `Median` = c(equal1RLE(), dif1RLE()),
      `No Normalization` = c(equal2RLE(), dif2RLE()),
      `PoissonSeq` = c(equal3RLE(), dif3RLE()),
      `Quantile` = c(equal4RLE(), dif4RLE()),
      #`RLE` = c(equal5RLE(), dif5RLE()),
      #`Total Count` = c(equal6RLE(), dif6RLE()),
      `TMM` = c(equal7RLE(), dif7RLE()),
      `TMMwsp` = c(equal8RLE(), dif8RLE()),
      `Upper-Quartile` = c(equal9RLE(), dif9RLE()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  
  
  
  #
  
  
  results_RLE1 <- reactive ({
    data <- results_RLE()
    data <- as.data.frame(results_RLE())
    data$rank <- c(1:nrow(results_RLE()))
    data$Method <- c(rep("RLE", times = nrow(results_RLE())))
    data$Gene <- c(rownames(results_RLE()))
    data
    
  })
  
  
  
  results_RLE11 <- reactive ({
    
    data <- results_RLE1()[order(results_RLE1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
  })
  
  
  
  datasetRLE1 <- reactive({                   
    data <- datasetRLE()
    data
  })
  
  
  RLE <- reactive({                   
    data <- datasetRLE()
    data
  })
  
  
  
  output$DEgenesRLE<- DT::renderDataTable({
    df2 <- datasetRLE()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DERLE1<- DT::renderDataTable({
    datasetRLE1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  ### RLE MDS
  
  MDS_RLE = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # RLE normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "RLE")) ## method = "RLE"
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    
    return(plot1())
    
  }
  
  results_MDSRLE <- reactive ({ MDS_RLE(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  ### RLE Glimma
  
  MDSGlimma_RLE = function(expression.matrix, metadata, design){
    
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # RLE normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "RLE")) ## method = "RLE"
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSRLE <- reactive ({ MDSGlimma_RLE(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  ### PoissonSeq ####
  
  DE_edgeR_PoissonSeq = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      cutoffvalue <- reactive ({ 
        
        L <-  as.numeric(min(colSums(expression.matrix))) 
        threshold <- as.numeric((10)/(L/10^6))
        threshold 
      })
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates }) # Use a random column to decide the cutoff threshold
      
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      invisible(capture.output( scaling.factor <- reactive ({  PoissonSeq::PS.Est.Depth(expression.matrix1()) }) ))
      dat.normed <- reactive ({ t(t(expression.matrix1())/scaling.factor()) })
      
      # PoissonSeq normalise
      expression.matrix.normalised <- reactive({
        data <- dat.normed()
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
    }
    return(res2())
    
  }
  
  results_PoissonSeq <- reactive ({ DE_edgeR_PoissonSeq(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resPoissonSeq10 <- reactive({ results_PoissonSeq() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resPoissonSeq11 <- reactive({    
    
    data <- resPoissonSeq10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resPoissonSeq10()))
    data
  })
  
  
  
  datPoissonSeq <- reactive({
    data <- results_PoissonSeq()
    data <- as.data.frame(results_PoissonSeq())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  PoissonSeqv2 <- reactive({                   
    data <- datPoissonSeq()
    data
  })
  
  
  PoissonSeqv3 <- reactive({
    data <- results_PoissonSeq()
    data <- as.data.frame(results_PoissonSeq())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  ##
  
  
  ### For Volcanos Home Page ###
  
  
  resPoissonSeqa <- reactive({ results_PoissonSeq() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resPoissonSeqb <- reactive({    
    data <- resPoissonSeqa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resPoissonSeqa()))
    data
  })
  
  
  
  
  
  datasetPoissonSeq <- reactive({
    data <- results_PoissonSeq()
    data <- as.data.frame(results_PoissonSeq())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  
  output$DEnumberPoissonSeq <- renderText({
    data <- length(c(row.names(datasetPoissonSeq())))
    data
  })    
  
  
  
  ###
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_PoissonSeq <- reactive({
    data <- c(row.names(datPoissonSeq()))
    data
  })
  
  ###
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownPoissonSeq <- reactive({ 
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpPoissonSeq <- reactive({ length(which(datPoissonSeq()$regulated == "up-regulated")) })
      DownPoissonSeq <- reactive({ length(which(datPoissonSeq()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `PoissonSeq` = c(UpPoissonSeq(), DownPoissonSeq()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownPoissonSeq <- DT::renderDataTable({
    df2 <- TableUpDownPoissonSeq()
    DT::datatable(df2)
    
    
  })
  
  
  GenesUpDownPoissonSeq <- reactive({ 
    
    
    ## NN
    equal2PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datNN()))) })
    dif2PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datNN()))))  })
    
    ## Median 
    equal1PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datMedian())))  })
    dif1PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datTC()))) })
    dif6PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datPoissonSeq()))) })
    dif3PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datDeseq()))) })
    dif11PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datPoissonSeq()))) })
    dif5PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datPoissonSeq())))) })
    
    ## TMM
    equal7PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datTMM()))) })
    dif7PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datTMMwsp()))) })
    dif8PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datSVA()))) })
    dif10PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datUQ()))) })
    dif9PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4PoissonSeq <- reactive({ length(which(row.names(datPoissonSeq()) %in% row.names(datQuantile()))) })
    dif4PoissonSeq <- reactive({ length(which(!(row.names(datPoissonSeq()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11PoissonSeq(), dif11PoissonSeq()),
      `Median` = c(equal1PoissonSeq(), dif1PoissonSeq()),
      `No Normalization` = c(equal2PoissonSeq(), dif2PoissonSeq()),
      #`PoissonSeq` = c(equal3PoissonSeq(), dif3PoissonSeq()),
      `Quantile` = c(equal4PoissonSeq(), dif4PoissonSeq()),
      `RLE` = c(equal5PoissonSeq(), dif5PoissonSeq()),
      #`Total Count` = c(equal6PoissonSeq(), dif6PoissonSeq()),
      `TMM` = c(equal7PoissonSeq(), dif7PoissonSeq()),
      `TMMwsp` = c(equal8PoissonSeq(), dif8PoissonSeq()),
      `Upper-Quartile` = c(equal9PoissonSeq(), dif9PoissonSeq()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  
  
  
  
  #
  
  results_PoissonSeq1 <- reactive ({
    data <- results_PoissonSeq()
    data <- as.data.frame(results_PoissonSeq())
    data$rank <- c(1:nrow(results_PoissonSeq()))
    data$Method <- c(rep("PoissonSeq", times = nrow(results_PoissonSeq())))
    data$Gene <- c(rownames(results_PoissonSeq()))
    data
  })
  
  
  results_PoissonSeq11 <- reactive ({
    
    data <- results_PoissonSeq1()[order(results_PoissonSeq1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
  })
  
  
  
  
  datasetPoissonSeq1 <- reactive({                   
    data <- datasetPoissonSeq()
    data
  })
  
  
  PoissonSeq <- reactive({                   
    data <- datasetPoissonSeq()
    data
  })
  
  
  
  output$DEgenesPoissonSeq<- DT::renderDataTable({
    df2 <- datasetPoissonSeq()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DEPoissonSeq1<- DT::renderDataTable({
    datasetPoissonSeq1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  ### Poisson MDS
  
  MDS_Poisson = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    invisible(capture.output( scaling.factor <- reactive ({  PoissonSeq::PS.Est.Depth(expression.matrix1()) }) ))
    dat.normed <- reactive ({ t(t(expression.matrix1())/scaling.factor()) })
    
    # PoissonSeq normalise
    expression.matrix.normalised <- reactive({
      data <- dat.normed()
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    
    
    return(plot1())
    
  }
  
  results_MDSPoisson <- reactive ({ MDS_Poisson(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  ### PoissonSeq Glimma
  
  MDSGlimma_PoissonSeq = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    invisible(capture.output( scaling.factor <- reactive ({  PoissonSeq::PS.Est.Depth(expression.matrix1()) }) ))
    dat.normed <- reactive ({ t(t(expression.matrix1())/scaling.factor()) })
    
    # PoissonSeq normalise
    expression.matrix.normalised <- reactive({
      data <- dat.normed()
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSPoissonSeq <- reactive ({ MDSGlimma_PoissonSeq(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  ### Median ####
  
  
  DE_edgeR_Median = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      cutoffvalue <- reactive ({               
        L <-  as.numeric(min(colSums(expression.matrix)))        
        threshold <- as.numeric((10)/(L/10^6))       
        threshold      })          
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
      
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      
      # Median normalise
      expression.matrix.normalised <- reactive({
        data <- limma::normalizeMedianValues(expression.matrix1())
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
    }
    return(res2())
    
  }
  
  results_Median <- reactive ({ DE_edgeR_Median(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resMedian10 <- reactive({ results_Median() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resMedian11 <- reactive({    
    
    data <- resMedian10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resMedian10()))
    data
  })
  
  
  
  
  datMedian <- reactive({
    data <- results_Median()
    data <- as.data.frame(results_Median())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  Medianv2 <- reactive({                   
    data <- datMedian()
    data
  })
  
  
  Medianv3 <- reactive({
    data <- results_Median()
    data <- as.data.frame(results_Median())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  
  ### For Volcanos Home Page ###
  
  
  resMediana <- reactive({ results_Median() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  
  resMedianb <- reactive({    
    data <- resMediana()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resMediana()))
    data
  })
  
  
  datasetMedian <- reactive({
    data <- results_Median()
    data <- as.data.frame(results_Median())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  
  output$DEnumberMedian <- renderText({
    data <- length(c(row.names(datasetMedian())))
    data
  })   
  
  ###
  
  
  ### Code for Venn Diagram ###
  
  Venn_Median <- reactive({
    data <- c(row.names(datMedian()))
    data
  })
  
  ###
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownMedian <- reactive({ 
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpMedian <- reactive({ length(which(datMedian()$regulated == "up-regulated")) })
      DownMedian <- reactive({ length(which(datMedian()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `Median` = c(UpMedian(), DownMedian()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data 
    }
    
    return(data)
  })
  
  
  output$TUpDownMedian <- DT::renderDataTable({
    df2 <- TableUpDownMedian()
    DT::datatable(df2)
    
    
  })
  
  
  
  
  GenesUpDownMedian <- reactive({ 
    
    
    ## NN
    equal2Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datNN()))) })
    dif2Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datNN()))))  })
    
    ## Media 
    equal1Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datMedian())))  })
    dif1Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datTC()))) })
    dif6Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datPoissonSeq()))) })
    dif3Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datDeseq()))) })
    dif11Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datRLE()))) })
    dif5Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datTMM()))) })
    dif7Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datTMMwsp()))) })
    dif8Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datSVA()))) })
    dif10Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datUQ()))) })
    dif9Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4Median <- reactive({ length(which(row.names(datMedian()) %in% row.names(datQuantile()))) })
    dif4Median <- reactive({ length(which(!(row.names(datMedian()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11Median(), dif11Median()),
      #`Median` = c(equal1Median(), dif1Median()),
      `No Normalization` = c(equal2Median(), dif2Median()),
      `PoissonSeq` = c(equal3Median(), dif3Median()),
      `Quantile` = c(equal4Median(), dif4Median()),
      `RLE` = c(equal5Median(), dif5Median()),
      #`Total Count` = c(equal6Median(), dif6Median()),
      `TMM` = c(equal7Median(), dif7Median()),
      `TMMwsp` = c(equal8Median(), dif8Median()),
      `Upper-Quartile` = c(equal9Median(), dif9Median()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  #
  
  
  
  results_Median1 <- reactive ({
    data <- results_Median()
    data <- as.data.frame(results_Median())
    data$rank <- c(1:nrow(results_Median()))
    data$Method <- c(rep("Median", times = nrow(results_Median())))
    data$Gene <- c(rownames(results_Median()))
    data
  })
  
  results_Median11 <- reactive ({
    
    data <- results_Median1()[order(results_Median1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
    
  })
  
  
  
  
  datasetMedian1 <- reactive({                   
    data <- datasetMedian()
    data
  })
  
  
  Median <- reactive({                   
    data <- datasetMedian()
    data
  })
  
  
  
  output$DEgenesMedian<- DT::renderDataTable({
    df2 <- datasetMedian()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DEMedian1<- DT::renderDataTable({
    datasetMedian1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  
  ### Median MDS
  
  MDS_Median = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    
    # Median normalise
    expression.matrix.normalised <- reactive({
      data <- limma::normalizeMedianValues(expression.matrix1())
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    
    return(plot1())
    
  }
  
  results_MDSMedian <- reactive ({ MDS_Median(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  ### Median Glimma
  
  MDSGlimma_Median = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    
    # Median normalise
    expression.matrix.normalised <- reactive({
      data <- limma::normalizeMedianValues(expression.matrix1())
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSMedian <- reactive ({ MDSGlimma_Median(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  ### TMMwsp ####
  
  DE_edgeR_TMMwsp = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      
      cutoffvalue <- reactive ({               
        L <-  as.numeric(min(colSums(expression.matrix)))        
        threshold <- as.numeric((10)/(L/10^6))       
        threshold      })          
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      # TMMwsp normalise
      expression.matrix.normalised <- reactive({
        data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "TMMwsp")) ## method = "TMMwsp"
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
      
      return(res2())
    }
  }
  
  results_TMMwsp <- reactive ({ DE_edgeR_TMMwsp(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resTMMwsp10 <- reactive({ results_TMMwsp() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resTMMwsp11 <- reactive({    
    
    data <- resTMMwsp10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resTMMwsp10()))
    data
  })
  
  datTMMwsp <- reactive({
    data <- results_TMMwsp()
    data <- as.data.frame(results_TMMwsp())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  TMMwspv2 <- reactive({                   
    data <- datTMMwsp()
    data
  })
  
  
  TMMwspv3 <- reactive({
    data <- results_TMMwsp()
    data <- as.data.frame(results_TMMwsp())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  
  
  ### For Volcanos Home Page ###
  
  
  resTMMwspa <- reactive({ results_TMMwsp() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resTMMwspb <- reactive({    
    data <- resTMMwspa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resTMMwspa()))
    data
  })
  
  
  
  datasetTMMwsp <- reactive({
    data <- results_TMMwsp()
    data <- as.data.frame(results_TMMwsp())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  output$DEnumberTMMwsp <- renderText({
    data <- length(c(row.names(datasetTMMwsp())))
    data
  })   
  
  
  
  ###
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_TMMwsp <- reactive({
    data <- c(row.names(datTMMwsp()))
    data
  })
  
  ###
  
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownTMMwsp <- reactive({ 
    
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpTMMwsp <- reactive({ length(which(datTMMwsp()$regulated == "up-regulated")) })
      DownTMMwsp <- reactive({ length(which(datTMMwsp()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `TMMwsp` = c(UpTMMwsp(), DownTMMwsp()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownTMMwsp <- DT::renderDataTable({
    df2 <- TableUpDownTMMwsp()
    DT::datatable(df2)
    
    
  })
  
  
  
  GenesUpDownTMMwsp <- reactive({ 
    
    
    ## NN
    equal2TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datNN()))) })
    dif2TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datNN()))))  })
    
    ## Media 
    equal1TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datMedian())))  })
    dif1TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datTC()))) })
    dif6TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datPoissonSeq()))) })
    dif3TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datDeseq()))) })
    dif11TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datRLE()))) })
    dif5TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datTMM()))) })
    dif7TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datTMMwsp()))) })
    dif8TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datSVA()))) })
    dif10TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datUQ()))) })
    dif9TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4TMMwsp <- reactive({ length(which(row.names(datTMMwsp()) %in% row.names(datQuantile()))) })
    dif4TMMwsp <- reactive({ length(which(!(row.names(datTMMwsp()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11TMMwsp(), dif11TMMwsp()),
      `Median` = c(equal1TMMwsp(), dif1TMMwsp()),
      `No Normalization` = c(equal2TMMwsp(), dif2TMMwsp()),
      `PoissonSeq` = c(equal3TMMwsp(), dif3TMMwsp()),
      `Quantile` = c(equal4TMMwsp(), dif4TMMwsp()),
      `RLE` = c(equal5TMMwsp(), dif5TMMwsp()),
      #`Total Count` = c(equal6TMMwsp(), dif6TMMwsp()),
      `TMM` = c(equal7TMMwsp(), dif7TMMwsp()),
      #`TMMwsp` = c(equal8TMMwsp(), dif8TMMwsp()),
      `Upper-Quartile` = c(equal9TMMwsp(), dif9TMMwsp()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  #
  
  
  
  results_TMMwsp1 <- reactive ({
    data <- results_TMMwsp()
    data <- as.data.frame( results_TMMwsp())
    data$rank <- c(1:nrow(results_TMMwsp()))
    data$Method <- c(rep("TMMwsp", times = nrow(results_TMMwsp())))
    data$Gene <- c(rownames(results_TMMwsp()))
    data
    
  })
  
  
  
  results_TMMwsp11 <- reactive ({
    
    data <- results_TMMwsp1()[order(results_TMMwsp1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
  })
  
  
  
  
  
  datasetTMMwsp1 <- reactive({                   
    data <- datasetTMMwsp()
    data
  })
  
  
  
  TMMwsp <- reactive({                   
    data <- datasetTMMwsp()
    data
  })
  
  
  
  output$DEgenesTMMwsp<- DT::renderDataTable({
    df2 <- datasetTMMwsp()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DETMMwsp1<- DT::renderDataTable({
    datasetTMMwsp1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  
  
  ### TMMwsp MDS
  
  MDS_TMMwsp = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # TMMwsp normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "TMMwsp")) ## method = "TMMwsp"
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    
    return(plot1())
    
  }
  
  results_MDSTMMwsp <- reactive ({ MDS_TMMwsp(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  ### TMMwsp Glimma
  
  MDSGlimma_TMMwsp = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # TMMwsp normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "TMMwsp")) ## method = "TMMwsp"
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSTMMwsp <- reactive ({ MDSGlimma_TMMwsp(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  ### Quantile #####
  
  
  DE_edgeR_Quantile = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      
      cutoffvalue <- reactive ({               
        L <-  as.numeric(min(colSums(expression.matrix)))        
        threshold <- as.numeric((10)/(L/10^6))       
        threshold      })          
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      # Quantile normalise
      expression.matrix.normalised <- reactive({
        data <- preprocessCore::normalize.quantiles(expression.matrix1())
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
      
      return(res2())
    }
  }
  
  results_Quantile <- reactive ({ DE_edgeR_Quantile(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resQuantile10 <- reactive({ results_Quantile() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resQuantile11 <- reactive({    
    
    data <- resQuantile10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resQuantile10()))
    data
  })
  
  
  datQuantile <- reactive({
    data <- results_Quantile()
    data <- as.data.frame(results_Quantile())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  Quantilev2 <- reactive({                   
    data <- datQuantile()
    data
  })
  
  
  Quantilev3 <- reactive({
    data <- results_Quantile()
    data <- as.data.frame(results_Quantile())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  
  ### For Volcanos Home Page ###
  
  
  resQuantilea <- reactive({ results_Quantile() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resQuantileb <- reactive({    
    data <- resQuantilea()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resQuantilea()))
    data
  })
  
  
  datasetQuantile <- reactive({
    data <- results_Quantile()
    data <- as.data.frame(results_Quantile())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  output$DEnumberQuantile <- renderText({
    data <- length(c(row.names(datasetQuantile())))
    data
  })   
  
  
  ###
  
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_Quantile <- reactive({
    data <- c(row.names(datQuantile()))
    data
  })
  
  ###
  
  
  
  ### Code for the Up Down Table
  
  
  TableUpDownQuantile <- reactive({ 
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpQuantile <- reactive({ length(which(datQuantile()$regulated == "up-regulated")) })
      DownQuantile <- reactive({ length(which(datQuantile()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `Quantile` = c(UpQuantile(), DownQuantile()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownQuantile <- DT::renderDataTable({
    df2 <- TableUpDownQuantile()
    DT::datatable(df2)
    
    
  })
  
  
  GenesUpDownQuantile <- reactive({ 
    
    
    ## NN
    equal2Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datNN()))) })
    dif2Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datNN()))))  })
    
    ## Media 
    equal1Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datMedian())))  })
    dif1Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datTC()))) })
    dif6Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datPoissonSeq()))) })
    dif3Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datDeseq()))) })
    dif11Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datRLE()))) })
    dif5Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datTMM()))) })
    dif7Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datTMMwsp()))) })
    dif8Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datSVA()))) })
    dif10Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datUQ()))) })
    dif9Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4Quantile <- reactive({ length(which(row.names(datQuantile()) %in% row.names(datQuantile()))) })
    dif4Quantile <- reactive({ length(which(!(row.names(datQuantile()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11Quantile(), dif11Quantile()),
      `Median` = c(equal1Quantile(), dif1Quantile()),
      `No Normalization` = c(equal2Quantile(), dif2Quantile()),
      `PoissonSeq` = c(equal3Quantile(), dif3Quantile()),
      #`Quantile` = c(equal4Quantile(), dif4Quantile()),
      `RLE` = c(equal5Quantile(), dif5Quantile()),
      #`Total Count` = c(equal6Quantile(), dif6Quantile()),
      `TMM` = c(equal7Quantile(), dif7Quantile()),
      `TMMwsp` = c(equal8Quantile(), dif8Quantile()),
      `Upper-Quartile` = c(equal9Quantile(), dif9Quantile()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  #
  
  
  
  results_Quantile1 <- reactive ({
    data <- results_Quantile()
    data <- as.data.frame(results_Quantile())
    data$rank <- c(1:nrow(results_Quantile()))
    data$Method <- c(rep("Quantile", times = nrow(results_Quantile())))
    data$Gene <- c(rownames(results_Quantile()))
    data <- as.data.frame(data)
    data
    
  })
  
  
  
  results_Quantile11 <- reactive ({
    data <- as.data.frame(results_Quantile1())
    data <- data[order(row.names(results_Quantile1())),]
    data
    
  })
  
  
  
  
  
  datasetQuantile1 <- reactive({                   
    data <- datasetQuantile()
    data
  })
  
  
  Quantile <- reactive({                   
    data <- datasetQuantile()
    data
  })
  
  
  output$DEgenesQuantile<- DT::renderDataTable({
    df2 <- datasetQuantile()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DEQuantile1<- DT::renderDataTable({
    datasetQuantile1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  ### Quantile MDS
  
  MDS_Quantile = function(expression.matrix, metadata, design){
    
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # Quantile normalise
    expression.matrix.normalised <- reactive({
      data <- preprocessCore::normalize.quantiles(expression.matrix1())
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <-  reactive ({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups") })
    
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    }) 
    
    
    
    return(plot1())
    
  }
  
  results_MDSQuantile <- reactive ({ MDS_Quantile(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  
  
  ### Quantile Glimma
  
  MDSGlimma_Quantile = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # Quantile normalise
    expression.matrix.normalised <- reactive({
      data <- preprocessCore::normalize.quantiles(expression.matrix1())
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSQuantile <- reactive ({ MDSGlimma_Quantile(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  ### Upper quartile ####
  
  
  
  DE_edgeR_UpperQuartile = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      
      cutoffvalue <- reactive ({               
        L <-  as.numeric(min(colSums(expression.matrix)))        
        threshold <- as.numeric((10)/(L/10^6))       
        threshold      })          
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      # Upper quartile normalise
      expression.matrix.normalised <- reactive({
        data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "upperquartile")) 
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
    }
    return(res2())
    
  }
  
  results_UpperQuartile <- reactive ({ DE_edgeR_UpperQuartile(Metadata1(),Informationtable1(), design()) })
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resUQ10 <- reactive({ results_UpperQuartile() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resUQ11 <- reactive({    
    
    data <- resUQ10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resUQ10()))
    data
  })
  
  
  datUQ <- reactive({
    data <- results_UpperQuartile()
    data <- as.data.frame(results_UpperQuartile())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  UQv2 <- reactive({                   
    data <- datUQ()
    data
  })
  
  UQv3 <- reactive({
    data <- results_UpperQuartile()
    data <- as.data.frame(results_UpperQuartile())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  ##
  
  
  ### For Volcanos Home Page ###
  
  
  resUQa <- reactive({ results_UpperQuartile() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resUQb <- reactive({ 
    data <- resUQa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resUQa()))
    data
  })
  
  
  
  datasetUpperQuartile <- reactive({
    data <- results_UpperQuartile()
    data <- as.data.frame(results_UpperQuartile())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  
  output$DEnumber_UpperQuartile <- renderText({
    data <- length(c(row.names(datasetUpperQuartile())))
    data
  })   
  
  
  ###
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_UQ <- reactive({
    data <- c(row.names(datUQ()))
    data
  })
  
  ###
  
  
  ### Code for the Up Down Table
  
  
  
  TableUpDownUQ <- reactive({ 
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpUQ <- reactive({ length(which(datUQ()$regulated == "up-regulated")) })
      DownUQ <- reactive({ length(which(datUQ()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `Upper Quartile` = c(UpUQ(), DownUQ()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownUQ <- DT::renderDataTable({
    df2 <- TableUpDownUQ()
    DT::datatable(df2)
    
    
  })
  
  
  
  GenesUpDownUQ <- reactive({ 
    
    
    ## NN
    equal2UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datNN()))) })
    dif2UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datNN()))))  })
    
    ## Media 
    equal1UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datMedian())))  })
    dif1UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datTC()))) })
    dif6UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datPoissonSeq()))) })
    dif3UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datDeseq()))) })
    dif11UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datRLE()))) })
    dif5UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datTMM()))) })
    dif7UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datTMMwsp()))) })
    dif8UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datSVA()))) })
    dif10UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datUQ()))) })
    dif9UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4UQ <- reactive({ length(which(row.names(datUQ()) %in% row.names(datQuantile()))) })
    dif4UQ <- reactive({ length(which(!(row.names(datUQ()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11UQ(), dif11UQ()),
      `Median` = c(equal1UQ(), dif1UQ()),
      `No Normalization` = c(equal2UQ(), dif2UQ()),
      `PoissonSeq` = c(equal3UQ(), dif3UQ()),
      `Quantile` = c(equal4UQ(), dif4UQ()),
      `RLE` = c(equal5UQ(), dif5UQ()),
      #`Total Count` = c(equal6UQ(),dif6UQ()),
      `TMM` = c(equal7UQ(), dif7UQ()),
      `TMMwsp` = c(equal8UQ(), dif8UQ()),
      #`Upper-Quartile` = c(equal9UQ(), dif9UQ()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  
  #
  
  
  
  
  results_UpperQuartile1 <- reactive ({
    data <- results_UpperQuartile()
    data <-as.data.frame(results_UpperQuartile())
    data$rank <- c(1:nrow(results_UpperQuartile()))
    data$Method <- c(rep("Upper Quartile", times = nrow(results_UpperQuartile())))
    data$Gene <- c(rownames(results_UpperQuartile()))
    data
  })
  
  
  results_UpperQuartile11 <- reactive ({
    
    data <- results_UpperQuartile1()[order(results_UpperQuartile1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
  })
  
  
  
  datasetUpperQuartile1 <- reactive({                   
    data <- datasetUpperQuartile()
    data
  })
  
  
  UpperQuartile <- reactive({                   
    data <- datasetUpperQuartile()
    data
  })
  
  
  
  output$DEgenes_UpperQuartile<- DT::renderDataTable({
    df2 <- datasetUpperQuartile()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DE_UpperQuartile1<- DT::renderDataTable({
    datasetUpperQuartile1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  
  ### UpperQuartile MDS
  
  MDS_UpperQuartile = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # Upper quartile normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "upperquartile")) 
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    plot2 <- reactive({
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    return(plot1())
    
  }
  
  results_MDSUpperQuartile <- reactive ({ MDS_UpperQuartile(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  
  ### UQ Glimma
  
  MDSGlimma_UQ  = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # Upper quartile normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()%*%diag(edgeR::normLibSizes(expression.matrix1(), method = "upperquartile")) 
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSUQ  <- reactive({ MDSGlimma_UQ(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  
  
  
  
  ### NN ####
  
  DE_edgeR_NN = function(expression.matrix, metadata, design){
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      res2() <- NULL
      
    } else {
      
      
      cutoffvalue <- reactive ({               
        L <-  as.numeric(min(colSums(expression.matrix)))        
        threshold <- as.numeric((10)/(L/10^6))       
        threshold      })          
      
      keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
      
      
      expression.matrix1 <- reactive ({
        data <- expression.matrix[keep(),]
        data
      })
      
      # NN normalise
      expression.matrix.normalised <- reactive({
        data <- expression.matrix1()
        rownames(data) <- base::rownames(expression.matrix1())
        colnames(data) <- base::colnames(expression.matrix1())
        data
      })
      
      # process using edgeR
      expression.matrix.for.de <- reactive({
        data <- round(expression.matrix.normalised())
        data <- data[apply(data, 1, sum) > 0, ]
        data
        
      })
      
      edger <- reactive ({
        data <- edgeR::DGEList(counts = expression.matrix.for.de())
        data <- edgeR::estimateDisp(data, design)
        data
      })
      
      edger.fit <- reactive ({ edgeR::glmQLFit(edger(), design) })
      edger.lrt <- reactive ({ edgeR::glmQLFTest(edger.fit(), contrast = as.numeric(strsplit(input$SelectContr, ",")[[1]])) }) # == B vs A 
      
      # extract results
      res1 <- reactive ({ edgeR::topTags(edger.lrt(), n = Inf)$table })
      
      res2 <- reactive ({
        data <- res1()
        data$DE <- res1()$FDR < input$FDR & abs(res1()$logFC) > input$LFC
        data
      })
    }
    return(res2())
    
  }
  
  results_NN <- reactive ({ DE_edgeR_NN(Metadata1(),Informationtable1(), design()) })
  
  
  
  ## DATA FOR THE VOLCANO PLOT ##
  
  resNN10 <- reactive({ results_NN() %>% 
      mutate(regulated = case_when(
        FDR < input$FDR & logFC < -input$LFC ~ "down", # significantly down
        FDR < input$FDR & logFC > input$LFC ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resNN11 <- reactive({    
    
    data <- resNN10()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resNN10()))
    data
  })
  
  
  datNN <- reactive({
    data <- results_NN()
    data <- as.data.frame(results_NN())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data <- data %>% filter(data$FDR < input$FDR)
    data <- data %>% filter(abs(data$logFC) > input$LFC)
    data
  })
  
  NNv2 <- reactive({                   
    data <- datNN()
    data
  })
  
  
  NNv3 <- reactive({
    data <- results_NN()
    data <- as.data.frame(results_NN())
    data$regulated <- with(data, ifelse(FDR < input$FDR & logFC < -input$LFC, "down-regulated", ifelse(FDR < input$FDR & logFC > input$LFC , "up-regulated", "Not differentially expressed")))
    data$regulated <- as.factor(data$regulated)
    data
  })
  
  ##
  
  
  ### For Volcanos Home Page ###
  
  
  resNNa <- reactive({ results_NN() %>% 
      mutate(regulated = case_when(
        FDR < input$padj_variable & logFC < -input$lfc_variable ~ "down", # significantly down
        FDR < input$padj_variable & logFC > input$lfc_variable ~ "up", # significantly up
        TRUE ~ "Not differentially expressed") # not significant
      ) })
  
  
  resNNb <- reactive({  
    data <- resNNa()
    req(data)
    data[,"regulated"] <- as.factor(data[,"regulated"])
    data$RefSeq <- c(row.names(resNNa()))
    data
  })
  
  
  datasetNN <- reactive({
    data <- results_NN()
    data <- as.data.frame(results_NN())
    data <- data %>% filter(data$FDR < input$padj_variable)
    data <- data %>% filter(abs(data$logFC) > input$lfc_variable)
  })
  
  output$DEnumberNN <- renderText({
    data <- length(c(row.names(datasetNN())))
    data
  })   
  
  
  ###
  
  
  
  
  ### Code for Venn Diagram ###
  
  Venn_NN <- reactive({
    data <- c(row.names(datNN()))
    data
  })
  
  ###
  
  
  
  ### Code for the Up Down Table
  
  
  
  TableUpDownNN <- reactive({ 
    
    
    infile <- input$SelectContr
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      UpNN <- reactive({ length(which(datNN()$regulated == "up-regulated")) })
      DownNN <- reactive({ length(which(datNN()$regulated == "down-regulated")) })
      
      
      data <- data.frame(
        `No Normalization` = c(UpNN(), DownNN()),
        row.names = c("Up-regulated genes", "Down-regulated genes"))
      
      data
    }
    return(data)
  })
  
  
  output$TUpDownNN<- DT::renderDataTable({
    df2 <- TableUpDownNN()
    DT::datatable(df2)
    
    
  })
  
  
  
  GenesUpDownNN <- reactive({ 
    
    
    ## NN
    equal2NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datNN()))) })
    dif2NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datNN()))))  })
    
    ## Media 
    equal1NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datMedian())))  })
    dif1NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datMedian())))) })
    
    ## Total count
    equal6NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datTC()))) })
    dif6NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datTC())))) })
    
    ## PoissonSeq
    equal3NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datPoissonSeq()))) })
    dif3NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datPoissonSeq())))) })
    
    ## Deseq2
    equal11NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datDeseq()))) })
    dif11NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datDeseq())))) })
    
    ## RLE
    equal5NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datRLE()))) })
    dif5NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datRLE())))) })
    
    ## TMM
    equal7NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datTMM()))) })
    dif7NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datTMM())))) })
    
    ## TMMwsp
    equal8NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datTMMwsp()))) })
    dif8NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datTMMwsp())))) })
    
    ## SVA
    equal10NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datSVA()))) })
    dif10NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datSVA())))) })
    
    ## UQ
    equal9NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datUQ()))) })
    dif9NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datUQ())))) })
    
    ## Quantile
    equal4NN <- reactive({ length(which(row.names(datNN()) %in% row.names(datQuantile()))) })
    dif4NN <- reactive({ length(which(!(row.names(datNN()) %in% row.names(datQuantile())))) })
    
    
    
    
    data <- data.frame(
      `DESeq2` = c(equal11NN(), dif11NN()),
      `Median` = c(equal1NN(), dif1NN()),
      #`No Normalization` = c(equal2NN(), dif2NN()),
      `PoissonSeq` = c(equal3NN(), dif3NN()),
      `Quantile` = c(equal4NN(), dif4NN()),
      `RLE` = c(equal5NN(), dif5NN()),
      #`Total Count` = c(equal6NN(), dif6NN()),
      `TMM` = c(equal7NN(), dif7NN()),
      `TMMwsp` = c(equal8NN(), dif8NN()),
      `Upper-Quartile` = c(equal9NN(), dif9NN()),
      row.names = c("Equal DE genes", "Different DE genes"))
    
    data 
    
  })
  
  
  
  #
  
  
  
  results_NN1 <- reactive ({
    data <- results_NN()
    data <- as.data.frame(results_NN())
    data$rank <- c(1:nrow(results_NN()))
    data$Method <- c(rep("No Normalization", times = nrow(results_NN())))
    data$Gene <- c(rownames(results_NN()))
    data
  })
  
  results_NN11 <- reactive ({
    
    data <- results_NN1()[order(results_NN1()[,"Gene"]),]
    data <- as.data.frame(data)
    data
    
  })
  
  
  
  
  
  datasetNN1 <- reactive({                   
    data <- datasetNN()
    data
  })
  
  
  NN <- reactive({                   
    data <- datasetNN()
    data
  })
  
  
  
  output$DEgenesNN<- DT::renderDataTable({
    df2 <- datasetNN()
    DT::datatable(df2)
    
    
  })
  
  
  
  output$DERNN1<- DT::renderDataTable({
    datasetNN1()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  ### NN MDS
  
  MDS_NN = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # NN normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data
    })
    
    
    
    
    y <- reactive ({ normLibSizes(edger(), method = "none") })
    
    
    col.group <- reactive ({ 
      data <- group()
      data <- factor(group())
      levels(data) <- paletteer_d("rcartocolor::Bold")
      data <- as.character(data)
      data
      
    })
    
    
    plot1 <- reactive({ 
      par(mfrow=c(1,2))
      limma::plotMDS(y(), top = 500, labels = group(), col= col.group(), cex = 1.4)
      title(main="Sample groups")
      
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Samples")
    })
    
    
    plot2 <- reactive({  
      limma::plotMDS(y(), top = 500, col= col.group(), cex = 1.4)
      title(main="Sample groups")
    })
    
    return(plot1())
    
  }
  
  results_MDSNN <- reactive ({ MDS_NN(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  
  ### NN Glimma
  
  MDSGlimma_NN  = function(expression.matrix, metadata, design){
    
    
    cutoffvalue <- reactive ({               
      L <-  as.numeric(min(colSums(expression.matrix)))        
      threshold <- as.numeric((10)/(L/10^6))       
      threshold      })          
    
    keep = reactive({ rowSums(cpm(expression.matrix) > cutoffvalue() ) >= input$Replicates })  
    
    expression.matrix1 <- reactive ({
      data <- expression.matrix[keep(),]
      data
    })
    
    # NN normalise
    expression.matrix.normalised <- reactive({
      data <- expression.matrix1()
      rownames(data) <- base::rownames(expression.matrix1())
      colnames(data) <- base::colnames(expression.matrix1())
      data
    })
    
    # process using edgeR
    expression.matrix.for.de <- reactive({
      data <- round(expression.matrix.normalised())
      data <- data[apply(data, 1, sum) > 0, ]
      data
      
    })
    
    
    groups.dff <- reactive({ 
      data <- as.data.frame(lapply(metadata, function(col) {
        if (is.factor(col)) {
          as.character(col)
        } else {
          col
        }
      }))
      
      data
    })
    
    
    edger <- reactive ({
      data <- edgeR::DGEList(counts = expression.matrix.for.de())
      data <- edgeR::estimateDisp(data, design)
      data
    })
    
    
    
    plotGlimma <- reactive({
      #Glimma::glimmaMDS(edger())
      Glimma::glMDSPlot(edger(), labels = group(), group = groups.dff() )
    })
    
    return(plotGlimma())
    
  }
  
  results_GlimmaMDSNN  <- reactive({ MDSGlimma_NN(Metadata1(),Informationtable1(), design()) })
  
  
  
  
  
  
  
  
  
  #### Gene Ontology Normal methods ####
  
  
  Go_results <-  reactive({ 
    infile <- input$KeyType
    req(infile)
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- enrichGO(gene = c(row.names(DataGeo())), OrgDb = GotermsOrgDb(), keyType = GotermsKeyType() , ont = GotermsOnto()) 
    }
    return(data) 
  })
  
  
  
  ## PLOTS FOR THE GENE ONTOLOGY
  
  ## BarPlot ##
  
  output$fit<- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- plot(barplot(Go_results(), showCategory = input$ShowCategory)) 
    }
    return(data)
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  ###
  
  
  ## Goplot ##
  
  output$goplot<- renderPlot({
    
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit2 <- NULL
      
    } else {
      
      fit2 <- goplot(Go_results(), showCategory = input$ShowCategory)
      fit2 }
    
    return(fit2)
  }, height = "auto", width = "auto" )
  
  ###
  
  
  
  ## Dotplot ##
  
  output$dotplot<- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      fit3 <- dotplot(Go_results(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
    
  },  height = "auto",  width = "auto" , res = 72  )
  
  ###
  
  
  
  ## cnetplot ##
  
  output$cnetplot <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- cnetplot1()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  fold_change_geneList <- reactive({ setNames(object = DataGeo()[,"logFC"], nm = row.names(DataGeo())) })
  
  fold_change_geneListDeseq2 <- reactive({ setNames(object = DataGeo()[,"log2FoldChange"], nm = row.names(DataGeo())) })
  
  
  
  
  ###
  
  cnetplot1 <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(Go_results(), showCategory = input$ShowCategory, color.params = list(foldChange=fold_change_geneListDeseq2())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(Go_results(), showCategory = input$ShowCategory, color.params = list(foldChange=fold_change_geneList())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  output$cnetplot2 <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- cnetplot2()
    }
    return(data)
    
  }, height = "auto", width = "auto"  )
  
  
  
  cnetplot2 <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(Go_results(), showCategory = input$ShowCategory, foldChange=fold_change_geneListDeseq2(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(Go_results(), showCategory = input$ShowCategory, foldChange=fold_change_geneList(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ## upsetplot ##
  
  output$upsetplot <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- upsetplot(Go_results(), n = input$ShowCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  ## heatplot ##
  
  output$heatplot <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      fit3 <- heatplot(Go_results(), showCategory = input$ShowCategory)
      fit3
      
    }
    return(fit3)
    
  }, height = "auto",  width = "auto" , res = 72 )
  
  
  
  ###
  
  ## heatplot2 ##
  
  
  output$heatplot2 <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- Heatplot2()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  
  Heatplot2 <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  heatplot(Go_results(), showCategory = input$ShowCategory, foldChange=fold_change_geneListDeseq2()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- heatplot(Go_results(), showCategory = input$ShowCategory, foldChange= fold_change_geneList()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ###
  
  
  
  # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. 
  # This makes it easier to identify functional modules:
  
  ora_analysis_bp <- reactive ({
    data <- pairwise_termsim(Go_results()) #, method = GotermsOnto()
    data
  })
  
  
  
  ## emapplot ##
  
  output$emapplot <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- emapplot(ora_analysis_bp(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  
  ## treeplot ##
  
  output$treeplot <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- treeplot(ora_analysis_bp(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72  )
  
  
  # p2 <- treeplot(edox2, hclust_method = "average")
  # (e.g., ‘average’, ‘complete’, ‘median’, ‘centroid’, etc., see also the document of the hclust() function). 
  
  ###
  
  
  
  #### Gene Ontology with GSEA ####
  
  EnrichGSEA <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      
      GSEA <- setNames(object = DataGeoGSEA()[,"log2FoldChange"], nm = row.names(DataGeoGSEA()))
      GSEA = sort(GSEA, decreasing = TRUE)
      data <- gseGO(gene = GSEA, OrgDb = GotermsOrgDb(), keyType = GotermsKeyType() , ont = GotermsOnto()) 
    } else {
      
      
      GSEA <- setNames(object = DataGeoGSEA()[,"logFC"], nm = row.names(DataGeoGSEA()))
      GSEA = sort(GSEA, decreasing = TRUE)
      data <- gseGO(gene = GSEA, OrgDb = GotermsOrgDb(), keyType = GotermsKeyType() , ont = GotermsOnto()) 
    }
    
    
    return(data)
    
  })
  
  
  ## 
  
  output$NameNorm17 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm18 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm19 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm20 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm21 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm22 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm23 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('GSEA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('GSEA with ',data, ':')
    }
    return(f)
    
  })
  
  
  
  ## PLOTS FOR THE GENE ONTOLOGY
  
  
  
  
  ## Dotplot ##
  
  output$dotplotGSEA<- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      fit3 <- dotplot(EnrichGSEA(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
    
  },  height = "auto",  width = "auto" , res = 72  )
  
  ###
  
  
  
  ## cnetplot ##
  
  output$cnetplotGSEA <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- cnetplot1GSEA()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  fold_change_geneListGSEA <- reactive({ setNames(object = DataGeoGSEA()[,"logFC"], nm = row.names(DataGeoGSEA())) })
  
  fold_change_geneListDeseq2GSEA <- reactive({ setNames(object = DataGeoGSEA()[,"log2FoldChange"], nm = row.names(DataGeoGSEA())) })
  
  
  
  
  ###
  
  cnetplot1GSEA <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichGSEA(), showCategory = input$ShowCategory, color.params = list(foldChange=fold_change_geneListDeseq2GSEA())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichGSEA(), showCategory = input$ShowCategory, color.params = list(foldChange=fold_change_geneListGSEA())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  output$cnetplot2GSEA <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- cnetplot2GSEA()
    }
    return(data)
    
  }, height = "auto", width = "auto"  )
  
  
  
  cnetplot2GSEA <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichGSEA(), showCategory = input$ShowCategory, foldChange=fold_change_geneListDeseq2GSEA(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichGSEA(), showCategory = input$ShowCategory, foldChange=fold_change_geneListGSEA(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ## upsetplot ##
  
  output$upsetplotGSEA <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- upsetplot(EnrichGSEA(), n = input$ShowCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## heatplot2 ##
  
  
  output$heatplot2GSEA <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- Heatplot2_GSEA()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  
  Heatplot2_GSEA <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  heatplot(EnrichGSEA(), showCategory = input$ShowCategory, foldChange=fold_change_geneListDeseq2GSEA()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- heatplot(EnrichGSEA(), showCategory = input$ShowCategory, foldChange= fold_change_geneListGSEA()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ###
  
  
  
  # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. 
  # This makes it easier to identify functional modules:
  
  ora_analysis_bp_GSEA <- reactive ({
    data <- pairwise_termsim(EnrichGSEA()) 
    data
  })
  
  
  
  ## emapplot ##
  
  output$emapplotGSEA <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- emapplot(ora_analysis_bp_GSEA(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  
  ## treeplot ##
  
  output$treeplotGSEA <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- treeplot(ora_analysis_bp_GSEA(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72  )
  
  
  
  
  ###
  
  
  
  
  
  
  #### Gene Ontology with Shannnon List ####
  
  
  
  DataGO <- reactive({
    
    
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == "Deseq2") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(Deseq2v2()))])
      DataGO <- as.data.frame(res9()[vetor,])
      
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(RLEv2()))])
      DataGO <- as.data.frame(results_RLE()[vetor,])
      
      
    } else if (infile == "TC") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(TCv2())),]))
      DataGO <- results_TC()[vetor,]
      
      
    } else if (infile == "Upper Quartile") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(UQv2())),]))
      DataGO <- results_UpperQuartile()[vetor,]
      
      
      
    } else if (infile == "TMM") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(TMMv2()))])
      DataGO <- as.data.frame(results_TMM()[vetor,])
      
    } else if (infile == "TMMwsp") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(TMMwspv2()))])
      DataGO <- as.data.frame(results_TMMwsp()[vetor,])
      
      
    } else if (infile == "PoissonSeq") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(PoissonSeqv2())),]))
      DataGO <- results_PoissonSeq()[vetor,]
      
    } else if (infile == "SVA") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(SVAv2())),]))
      DataGO <- results_SVA()[vetor,]
      
      
    } else if (infile == "Median") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(Medianv2())),]))
      DataGO <- results_Median()[vetor,]
      
      
    } else if (infile == "Quantile") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(Quantilev2())),]))
      DataGO <- results_Quantile()[vetor,]
      
      
    } else {
      vetor <- c(row.names(tentativa()[1:length(row.names(NNv2())),]))
      DataGO <- results_NN()[vetor,]
      
    }
    
    return(DataGO)
    
  })
  
  
  
  Go_results_Shannon <-  reactive({ 
    infile <- input$Ei
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      data <- enrichGO(gene = c(row.names(DataGO())), OrgDb = GotermsOrgDb(), keyType = GotermsKeyType() , ont = GotermsOnto()) 
    }
    
    return(data)
  })
  
  
  ###
  
  
  
  
  ## BarPlot Shannon ##
  
  output$fitShannon<- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      data <- plot(barplot(Go_results_Shannon(), showCategory = input$ShowCategory)) # slider option
    }
    return(data)
  })
  
  
  ###
  
  
  
  
  
  
  ## Goplot ##
  
  output$goplotShannon<- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- plot(goplot(Go_results_Shannon(), showCategory = input$ShowCategory))
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## Dotplot ##
  
  output$dotplotShannon<- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- dotplot(Go_results_Shannon(), showCategory = input$ShowCategory)
      fit3
    } 
    
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## cnetplot ##
  
  output$cnetplotShannon <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- cnetplot1Shannon()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  fold_change_geneList_Shannon <- reactive({ setNames(object = DataGO()[,"logFC"], nm = row.names(DataGO())) })
  
  fold_change_geneListDeseq2_Shannon <- reactive({ setNames(object = DataGO()[,"log2FoldChange"], nm = row.names(DataGO())) })
  
  
  
  
  ###
  
  cnetplot1Shannon <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(Go_results_Shannon(), showCategory = input$ShowCategory, color.params = list(foldChange=fold_change_geneListDeseq2_Shannon())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(Go_results_Shannon(), showCategory = input$ShowCategory, color.params = list(foldChange=fold_change_geneList_Shannon())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  output$cnetplot2Shannon <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      data <- cnetplot2Shannon()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  
  cnetplot2Shannon <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(Go_results_Shannon(), showCategory = input$ShowCategory, foldChange=fold_change_geneListDeseq2_Shannon(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(Go_results_Shannon(), showCategory = input$ShowCategory, foldChange=fold_change_geneList_Shannon(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  ## upsetplot ##
  
  output$upsetplotShannon <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- upsetplot(Go_results_Shannon(), n = input$ShowCategory)
      fit3
    }
    return(fit3)
  }, height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  ## heatplot ##
  
  output$heatplotShannon <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- heatplot(Go_results_Shannon(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  
  ###
  
  ## heatplot2 ##
  
  
  output$heatplot2Shannon <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      data <- NULL
      
    } else {
      
      data <- Heatplot2Shannon()
    }
    return(data)
  }, height = "auto",  width = "auto" , res = 72 )
  
  
  
  Heatplot2Shannon <- reactive({
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  heatplot(Go_results_Shannon(), showCategory = input$ShowCategory, foldChange=fold_change_geneListDeseq2_Shannon()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- heatplot(Go_results_Shannon(), showCategory = input$ShowCategory, foldChange= fold_change_geneList_Shannon()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ###
  
  
  
  # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. 
  # This makes it easier to identify functional modules:
  
  ora_analysis_bp_Shannon <- reactive ({
    data <- pairwise_termsim(Go_results_Shannon()) #, method = GotermsOnto()
    data
  })
  
  
  
  ## emapplot ##
  
  output$emapplotShannon <- renderPlot({
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- emapplot(ora_analysis_bp_Shannon(), showCategory = input$ShowCategory)
      fit3
    }
    return(fit3)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  
  ## treeplot ##
  
  output$treeplotShannon <- renderPlot({
    
    
    infile <- input$KeyType
    req(infile)
    
    
    if (infile == 'Select the KeyType') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- treeplot(ora_analysis_bp_Shannon(), showCategory = input$ShowCategory)
      fit3
    } 
    return(fit3)
  }, height = "auto",  width = "auto" , res = 72 )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### Gene Ontology options ####
  
  
  output$OrgDb <-  renderUI ({
    data <- pickerInput(inputId = 'OrgDb', 
                        label = 'OrgDb terms:', 
                        choices = list('Select the OrgDb term','org.Mm.eg.db', 'org.Rn.eg.db', 'org.Hs.eg.db'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  # In order to get the orh.Mm.eg.db with " "
  GotermsOrgDb <- reactive({ 
    data <- input$OrgDb 
    data
  })
  
  
  
  output$KeyType <-  renderUI ({
    data <- pickerInput(inputId = 'KeyType', 
                        label = 'keyType:', 
                        choices = list('Select the KeyType','ENSEMBL', 'ENTREZID', 'IPI', 'PROSITE', 'ALIAS', 'CHR', 'CHRLOC', 'CHRLOCEND', 'ENZYME', 'MAP', 'PATH', 
                                       'PMID', 'REFSEQ', 'UNIGENE', 'SYMBOL', 'ENSEMBLPROT', 'ENSEMBLTRANS', 'GENENAME', 'UNIPROT', 
                                       'GO', 'EVIDENCE', 'ONTOLOGY', 'GOALL', 'EVIDENCEALL', 'ONTOLOGYALL', 'OMIM', 'UCSCKG', 'ACCNUM','MGI','PFAM'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  
  GotermsKeyType <- reactive({ 
    data <- input$KeyType 
    data
  })
  
  
  output$Ontology <-  renderUI ({
    data <- pickerInput(inputId = 'ontology', 
                        label = 'GO terms:', 
                        choices = list('Select the GO term','BP', 'MF', 'CC', 'ALL'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  
  GotermsOnto <- reactive({ 
    data <- input$ontology 
    data
  })
  
  
  
  output$NormOnto <-  renderUI ({
    data <- pickerInput(inputId = 'NormForOntology', 
                        label = 'Normalization methods:', 
                        choices = list('Select the Normalization method','Deseq2','RLE','TMM','TMMwsp','PoissonSeq', 'Median', 'Quantile', 'Upper Quartile', 'No Normalization'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  
  
  
  ### TEXT TO have to show which normalization method is currently being used in the Gene Ontology and KEGG
  
  
  output$NameNorm <- renderText({ 
    
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm1 <- renderText({ 
    
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm2 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm3 <- renderText({ 
    
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm4 <- renderText({ 
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm5 <- renderText({ 
    infile <-  input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm6 <- renderText({ 
    infile <-  input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm7 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm8 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm9 <- renderText({ 
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm10 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm11 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm12 <- renderText({ 
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm13 <- renderText({ 
    infile <- input$NormForOntology
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormForOntology
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  
  output$NameNorm14 <- renderText({ 
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm15 <- renderText({ 
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  output$NameNorm16 <- renderText({ 
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == 'Select the Normalization method') {
      f <- paste('ORA with: ')
      
    } else {
      data <- input$NormKegg
      
      f <- paste('ORA with ',data, ':')
    }
    return(f)
    
  })
  
  ###
  
  
  
  
  DataGeo <- reactive({
    
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      DataGeo <- Deseq2()
      
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataGeo <- RLE()
      
      
    } else if (infile == "TC") {
      
      DataGeo <- TC()
      
      
    } else if (infile == "Upper Quartile") {
      
      DataGeo <- UpperQuartile()
      
      
      
    } else if (infile == "TMM") {
      
      DataGeo <- TMM()
      
    } else if (infile == "TMMwsp") {
      
      DataGeo <- TMMwsp()
      
      
    } else if (infile == "PoissonSeq") {
      
      DataGeo <- PoissonSeq()
      
    } else if (infile == "SVA") {
      
      DataGeo <- SVA()
      
      
    } else if (infile == "Median") {
      
      DataGeo <- Median()
      
    } else if (infile == "Cyclic loess") {
      
      DataGeo <- Cloess()
      
    } else if (infile == "loess") {
      
      DataGeo <- loess()
      
      
    } else if (infile == "Quantile") {
      
      DataGeo <- Quantile()
      
      
    } else {
      DataGeo <- NN()
      
    }
    
    
  })
  
  
  
  
  DataGeoGSEA <- reactive({
    
    
    infile <- input$NormForOntology
    req(infile)
    
    if (infile == "Deseq2") {
      
      DataGeo <- as.data.frame(res7())
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataGeo <- as.data.frame(results_RLE())
      
      
    } else if (infile == "Upper Quartile") {
      
      DataGeo <- as.data.frame(results_UpperQuartile())
      
      
    } else if (infile == "TMM") {
      
      DataGeo <- as.data.frame(results_TMM())
      
    } else if (infile == "TMMwsp") {
      
      DataGeo <- as.data.frame(results_TMMwsp())
      
      
    } else if (infile == "PoissonSeq") {
      
      DataGeo <- as.data.frame(results_PoissonSeq())
      
      
    } else if (infile == "Median") {
      
      DataGeo <- as.data.frame(results_Median())
      
      
    } else if (infile == "Quantile") {
      
      DataGeo <- as.data.frame(results_Quantile())
      
      
    } else {
      DataGeo <- as.data.frame(results_NN())
      
    }
    
    return(DataGeo)
  })
  
  
  
  #### KEGG options ####
  
  
  output$KEGGkeytype <-  renderUI ({
    data <- pickerInput(inputId = 'KEGGkeytype', 
                        label = 'KEGG GENES Database:', 
                        choices = list('Select the database','kegg', 'ncbi-geneid', 'ncib-proteinid', 'uniprot'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  
  KEGGterms <- reactive({ 
    data <- input$KEGGkeytype
    data
  })
  
  
  
  output$KEGGorganism <-  renderUI ({
    data <- pickerInput(inputId = 'KEGGorganism', 
                        label = 'KEGG organism:', 
                        choices = list('Select the organism','hsa', 'mmu'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  
  KEGGorg <- reactive({ 
    data <- input$KEGGorganism
    data
  })
  
  
  
  
  output$NormKegg <-  renderUI ({
    data <- pickerInput(inputId = 'NormKegg', 
                        label = 'Normalization methods:', 
                        choices = list('Select the Normalization method','Deseq2','RLE','TMM','TMMwsp','PoissonSeq', 'Median', 'Quantile', 'Upper Quartile', 'No Normalization'),
                        options = list(`actions-box` = TRUE),multiple = F)
    
    data
    
  })
  
  
  
  
  
  DataGeo1 <- reactive({
    
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      DataGeo1 <- Deseq2()
      
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataGeo1 <- RLE()
      
      
    } else if (infile == "TC") {
      
      DataGeo1 <- TC()
      
      
    } else if (infile == "Upper Quartile") {
      
      DataGeo1 <- UpperQuartile()
      
      
      
    } else if (infile == "TMM") {
      
      DataGeo1 <- TMM()
      
    } else if (infile == "TMMwsp") {
      
      DataGeo1 <- TMMwsp()
      
      
    } else if (infile == "PoissonSeq") {
      
      DataGeo1 <- PoissonSeq()
      
    } else if (infile == "SVA") {
      
      DataGeo1 <- SVA()
      
      
    } else if (infile == "Median") {
      
      DataGeo1 <- Median()
      
    } else if (infile == "Cyclic loess") {
      
      DataGeo1 <- Cloess()
      
    } else if (infile == "loess") {
      
      DataGeo1 <- loess()
      
      
    } else if (infile == "Quantile") {
      
      DataGeo1 <- Quantile()
      
      
    } else {
      DataGeo1 <- NN()
      
    }
    return(DataGeo1)
    
  })
  
  
  DataKegg <- reactive({
    
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      DataKegg <- setNames(object =  Deseq2()[,"log2FoldChange"], nm = row.names(Deseq2()))
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataKegg <- setNames(object =  RLE()[,"logFC"], nm = row.names(RLE()))
      
      
    } else if (infile == "TC") {
      
      DataKegg <- TC()$logFC
      names(DataKegg) <- row.names(TC())
      
      
    } else if (infile == "Upper Quartile") {
      
      DataKegg <- UpperQuartile()$logFC
      names(DataKegg) <- row.names(UpperQuartile())
      
      
      
    } else if (infile == "TMM") {
      
      DataKegg <- TMM()$logFC
      names(DataKegg) <- row.names(TMM())
      
    } else if (infile == "TMMwsp") {
      
      DataKegg <- TMMwsp()$logFC
      names(DataKegg) <- row.names(TMMwsp())
      
      
    } else if (infile == "PoissonSeq") {
      
      DataKegg <- PoissonSeq()$logFC
      names(DataKegg) <- row.names(PoissonSeq())
      
    } else if (infile == "SVA") {
      
      DataKegg <- SVA()$logFC
      names(DataKegg) <- row.names(SVA())
      
      
    } else if (infile == "Median") {
      
      DataKegg <- Median()$logFC
      names(DataKegg) <- row.names(Median())
      
    } else if (infile == "Cyclic loess") {
      
      DataKegg <- Cloess()$logFC
      names(DataKegg) <- row.names(Cloess())
      
      
    } else if (infile == "Quantile") {
      
      DataKegg <- setNames(object =  Quantile()[,"logFC"], nm = row.names(Quantile()))
      
      
    } else {
      DataKegg <- NN()$logFC
      names(DataKegg) <- row.names(NN())
      
    }
    
    return(DataKegg)
    
  })
  
  
  ## Enrich KEGG ####
  
  ids <- reactive({
    
    ids <- bitr(names(DataKegg()), fromType = GotermsKeyType(), toType = "ENTREZID", OrgDb = GotermsOrgDb())
    ids
    
  })
  
  
  df2 <- reactive({
    dedup_ids <- ids()[!duplicated(ids()[c(GotermsKeyType())]),]
    
    df2 = DataGeo1()[row.names(DataGeo1()) %in% dedup_ids[,1],]
    df2$names = dedup_ids[,2]
    df2
  })
  
  
  
  EnrichKegg <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      
      kegg <- setNames(object = df2()[,"log2FoldChange"], nm = df2()[,"names"])
      kegg = sort(kegg, decreasing = TRUE)
      kk1 <- enrichKEGG(gene =  names(kegg) ,  organism=KEGGorg(), pvalueCutoff = 0.05, keyType = KEGGterms())
      
    } else {
      
      
      kegg <- setNames(object = df2()[,"logFC"], nm = df2()[,"names"])
      kegg = sort(kegg, decreasing = TRUE)
      kk1 <- enrichKEGG(gene =  names(kegg),  organism=KEGGorg(), pvalueCutoff = 0.05, keyType = KEGGterms())
    }
    
    
    return(kk1)
    
  })
  
  
  
  
  ## PLOTS FOR THE KEGG
  
  ## BarPlot ## 
  
  output$KEGGbarplot<- renderPlot({
    
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      
      data <- plot(barplot(EnrichKegg(), showCategory = input$ShowKEGGCategory)) # slider option
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  ###
  
  
  
  ## Goplot ##
  
  output$KEGGgoplot <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      
      data <- plot(goplot(EnrichKegg(), showCategory = input$ShowKEGGCategory))
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## Dotplot ##
  
  output$KEGGdotplot<- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- dotplot(EnrichKegg(), showCategory = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## cnetplot ##
  
  output$KEGGcnetplot <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGcnetplot1()
    }
    return(data)
  }, height = "auto",  width = "auto" , res = 72 )
  
  
  fold_change_geneList_KEGG <- reactive({ 
    data <- setNames(object = df2()[,"logFC"], nm = df2()[,"names"])
    data
  })
  
  fold_change_geneListDeseq2_KEGG <- reactive({ 
    data <- setNames(object = df2()[,"log2FoldChange"], nm = df2()[,"names"])
    data })
  
  
  
  
  ###
  
  KEGGcnetplot1 <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichKegg(), showCategory = input$ShowKEGGCategory, color.params = list(foldChange=fold_change_geneListDeseq2_KEGG())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichKegg(), showCategory = input$ShowKEGGCategory, color.params = list(foldChange=fold_change_geneList_KEGG())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  output$KEGGcnetplot21 <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGcnetplot2()
    }
    return(data)
    
  }, height = "auto", width = "auto", res = 72 )
  
  
  
  KEGGcnetplot2 <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichKegg(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneListDeseq2_KEGG(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichKegg(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneList_KEGG(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ## upsetplot ##
  
  output$KEGGupsetplot <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- upsetplot(EnrichKegg(), n = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72)
  
  ###
  
  
  ## heatplot2 ##
  
  
  output$KEGGheatplot <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGHeatplot2()
    }
    return(data)
  },  height = "auto",  width = "auto" , res = 72)
  
  
  
  KEGGHeatplot2 <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  heatplot(EnrichKegg(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneListDeseq2_KEGG()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- heatplot(EnrichKegg(), showCategory = input$ShowKEGGCategory, foldChange= fold_change_geneList_KEGG()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  
  KEGG_ora_analysis_bp <- reactive ({
    data <- pairwise_termsim(EnrichKegg())
    data
  })
  
  
  
  ## emapplot ##
  
  output$KEGGemapplot <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- emapplot(KEGG_ora_analysis_bp(), showCategory = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  
  ## treeplot ##
  
  output$KEGGtreeplot <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- treeplot(KEGG_ora_analysis_bp(), showCategory = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72  )
  
  
  
  
  ###
  
  
  
  
  
  ## Enrich KEGG GSEA ####
  # This code is functional but given the low terms in several datasets I did not insert it in the application
  # If wanted just need to call the plots in the ui
  
  
  DataGeo1GSEA <- reactive({
    
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      DataGeo <- res7()
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataGeo <- results_RLE()
      
      
    } else if (infile == "Upper Quartile") {
      
      DataGeo <- results_UpperQuartile()
      
      
    } else if (infile == "TMM") {
      
      DataGeo <- as.data.frame(results_TMM())
      
    } else if (infile == "TMMwsp") {
      
      DataGeo <- as.data.frame(results_TMMwsp())
      
      
    } else if (infile == "PoissonSeq") {
      
      DataGeo <- as.data.frame(results_PoissonSeq())
      
      
    } else if (infile == "Median") {
      
      DataGeo <- as.data.frame(results_Median())
      
      
    } else if (infile == "Quantile") {
      
      DataGeo <- as.data.frame(results_Quantile())
      
      
    } else {
      DataGeo <- as.data.frame(results_NN())
      
    }
    
    return(DataGeo)
    
  })
  
  
  DataKeggGSEA <- reactive({
    
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      DataKegg <- setNames(object =  res7()[,"log2FoldChange"], nm = row.names(res7()))
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      
      DataKegg <- setNames(object = results_RLE[,"logFC"], nm = row.names(results_RLE))
      
      
    } else if (infile == "Upper Quartile") {
      
      DataKegg <- results_UpperQuartile()$logFC
      names(DataKegg) <- row.names(results_UpperQuartile())
      
      
      
    } else if (infile == "TMM") {
      
      DataKegg <- results_TMM()$logFC
      names(DataKegg) <- row.names(results_TMM())
      
    } else if (infile == "TMMwsp") {
      
      DataKegg <- results_TMMwsp()$logFC
      names(DataKegg) <- row.names( results_TMMwsp())
      
      
    } else if (infile == "PoissonSeq") {
      
      DataKegg <- results_PoissonSeq()$logFC
      names(DataKegg) <- row.names( results_PoissonSeq())
      
    } else if (infile == "Median") {
      
      DataKegg <- results_Median()$logFC
      names(DataKegg) <- row.names(results_Median())
      
      
    } else if (infile == "Quantile") {
      
      DataKegg <- setNames(object =  results_Quantile()[,"logFC"], nm = row.names(results_Quantile()))
      
      
    } else {
      DataKegg <- results_NN()$logFC
      names(DataKegg) <- row.names(results_NN())
      
    }
    
    return(DataKegg)
    
  })
  
  
  
  
  idsGSEA <- reactive({
    
    ids <- bitr(names(DataKeggGSEA()), fromType = GotermsKeyType(), toType = "ENTREZID", OrgDb = GotermsOrgDb())
    ids
    
  })
  
  
  df2GSEA <- reactive({
    dedup_ids <- idsGSEA()[!duplicated(idsGSEA()[c(GotermsKeyType())]),]
    
    df2 = DataGeo1GSEA()[row.names(DataGeo1GSEA()) %in% dedup_ids[,1],]
    df2$names = dedup_ids[,2]
    df2
  })
  
  
  
  EnrichKeggGSEA <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      
      kegg <- setNames(object = df2GSEA()[,"log2FoldChange"], nm = df2GSEA()[,"names"])
      kegg = sort(kegg, decreasing = TRUE)
      kk1 <- gseKEGG(geneList  =  kegg ,  organism=KEGGorg(), pvalueCutoff = 0.05, keyType = KEGGterms())
      
    } else {
      
      kegg <- setNames(object = df2GSEA()[,"logFC"], nm = df2GSEA()[,"names"])
      kegg = sort(kegg, decreasing = TRUE)
      kk1 <- gseKEGG(geneList = kegg,  organism=KEGGorg(), pvalueCutoff = 0.05, keyType = KEGGterms())
    }
    
    
    return(kk1)
    
  })
  
  
  
  
  ## PLOTS FOR THE KEGG GSEA
  
  
  ## Dotplot ##
  
  output$KEGGdotplotGSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- dotplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## cnetplot ##
  
  output$KEGGcnetplotGSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGcnetplot1GSEA()
    }
    return(data)
  }, height = "auto",  width = "auto" , res = 72 )
  
  
  fold_change_geneList_KEGGGSEA <- reactive({ 
    data <- setNames(object = df2GSEA()[,"logFC"], nm = df2GSEA()[,"names"])
    data
  })
  
  fold_change_geneListDeseq2_KEGGGSEA <- reactive({ 
    data <- setNames(object = df2GSEA()[,"log2FoldChange"], nm = df2GSEA()[,"names"])
    data })
  
  
  
  
  ###
  
  KEGGcnetplot1GSEA <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory, color.params = list(foldChange=fold_change_geneListDeseq2_KEGGGSEA())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory, color.params = list(foldChange=fold_change_geneList_KEGGGSEA())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  output$KEGGcnetplot21GSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGcnetplot2GSEA()
    }
    return(data)
    
  }, height = "auto", width = "auto", res = 72 )
  
  
  
  KEGGcnetplot2GSEA <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneListDeseq2_KEGGGSEA(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneList_KEGGGSEA(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ## upsetplot ##
  
  output$KEGGupsetplotGSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- upsetplot(EnrichKeggGSEA(), n = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72)
  
  ###
  
  
  ## heatplot2 ##
  
  
  output$KEGGheatplotGSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGHeatplot2GSEA()
    }
    return(data)
  },  height = "auto",  width = "auto" , res = 72)
  
  
  
  KEGGHeatplot2GSEA <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  heatplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneListDeseq2_KEGG()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- heatplot(EnrichKeggGSEA(), showCategory = input$ShowKEGGCategory, foldChange= fold_change_geneList_KEGG()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  
  KEGG_ora_analysis_bp_GSEA <- reactive ({
    data <- pairwise_termsim(EnrichKeggGSEA())
    data
  })
  
  
  
  ## emapplot ##
  
  output$KEGGemapplotGSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- emapplot(KEGG_ora_analysis_bp_GSEA(), showCategory = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  
  ## treeplot ##
  
  output$KEGGtreeplotGSEA <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- treeplot(KEGG_ora_analysis_bp_GSEA(), showCategory = input$ShowKEGGCategory)
      fit3
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72  )
  
  
  
  
  ###
  
  
  
  
  
  ## KEGG with Shannon ####
  
  
  
  
  DataKEGGShannon <- reactive({
    
    
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == "Deseq2") {
      
      
      DataKegg <- setNames(object =  DataShannonKEGG()[,"log2FoldChange"], nm = row.names(DataShannonKEGG()))
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataKegg <- setNames(object =   DataShannonKEGG()[,"logFC"], nm = row.names( DataShannonKEGG()))
      
    } else if (infile == "TC") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(TCv2())),]))
      DataGO <- results_TC()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "Upper Quartile") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(UQv2())),]))
      DataGO <- results_UpperQuartile()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "TMM") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(TMMv2()))])
      DataGO <- as.data.frame(results_TMM()[vetor,])
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "TMMwsp") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(TMMwspv2()))])
      DataGO <- as.data.frame(results_TMMwsp()[vetor,])
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "PoissonSeq") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(PoissonSeqv2())),]))
      DataGO <- results_PoissonSeq()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "SVA") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(SVAv2())),]))
      DataGO <- results_SVA()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "Median") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(Medianv2())),]))
      DataGO <- results_Median()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else if (infile == "Quantile") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(Quantilev2())),]))
      DataGO <- results_Quantile()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    } else {
      vetor <- c(row.names(tentativa()[1:length(row.names(NNv2())),]))
      DataGO <- results_NN()[vetor,]
      
      DataKegg <- setNames(object =  DataGO[,"logFC"], nm = row.names(DataGO))
      
    }
    
    return(DataKegg)
    
  })
  
  
  
  
  
  DataShannonKEGG <- reactive({
    
    
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == "Deseq2") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(Deseq2v2()))])
      DataGO <- as.data.frame(res9()[vetor,])
      
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(RLEv2()))])
      DataGO <- as.data.frame(results_RLE()[vetor,])
      
      
    } else if (infile == "TC") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(TCv2())),]))
      DataGO <- results_TC()[vetor,]
      
      
    } else if (infile == "Upper Quartile") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(UQv2())),]))
      DataGO <- results_UpperQuartile()[vetor,]
      
      
      
    } else if (infile == "TMM") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(TMMv2()))])
      DataGO <- as.data.frame(results_TMM()[vetor,])
      
    } else if (infile == "TMMwsp") {
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(TMMwspv2()))])
      DataGO <- as.data.frame(results_TMMwsp()[vetor,])
      
      
    } else if (infile == "PoissonSeq") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(PoissonSeqv2())),]))
      DataGO <- results_PoissonSeq()[vetor,]
      
    } else if (infile == "SVA") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(SVAv2())),]))
      DataGO <- results_SVA()[vetor,]
      
      
    } else if (infile == "Median") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(Medianv2())),]))
      DataGO <- results_Median()[vetor,]
      
      
    } else if (infile == "Quantile") {
      
      vetor <- c(row.names(tentativa()[1:length(row.names(Quantilev2())),]))
      DataGO <- results_Quantile()[vetor,]
      
      
    } else {
      vetor <- c(row.names(tentativa()[1:length(row.names(NNv2())),]))
      DataGO <- results_NN()[vetor,]
      
    }
    
    return(DataGO)
    
  })
  
  
  idsShannon <- reactive({
    
    ids <- bitr(names(DataKEGGShannon()), fromType = GotermsKeyType(), toType = "ENTREZID", OrgDb = GotermsOrgDb())
    ids
    
  })
  
  
  
  ddd <- reactive({
    
    dedup_ids <- idsShannon()[!duplicated(idsShannon()[c(GotermsKeyType())]),] 
    
    data <- DataShannonKEGG()[na.omit(match(row.names(DataShannonKEGG()), dedup_ids[,1])), ]
    data$names <- dedup_ids[,2]
    data
  })
  
  dedup_ids <- reactive({ idsShannon()[!duplicated(idsShannon()[c(GotermsKeyType())]),] })
  
  bb <- reactive({  dedup_ids()[,2] })
  
  
  output$Tr <- DT::renderDataTable({
    df2 <- as.data.frame(bb())
    DT::datatable(df2)
    
  })
  
  
  output$Trr <- DT::renderDataTable({
    df2 <- as.data.frame(fold_change_geneList_KEGG_Shannon())
    DT::datatable(df2)
    
  })
  
  fold_change_geneList_KEGG_Shannon <- reactive({ 
    df2 <- setNames(object = ddd()[,"logFC"], nm = dedup_ids()[,2]) 
    df2 = sort(df2, decreasing = TRUE)
    df2
  })
  
  fold_change_geneListDeseq2_KEGG_Shannon <- reactive({ 
    data <- setNames(object = ddd()[,"log2FoldChange"], nm =  dedup_ids()[,2]) 
    data = sort(data, decreasing = TRUE)
    data
  })
  
  
  EnrichKeggShannon <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    
    if (infile == "Deseq2") {
      
      kk1 <- enrichKEGG(gene =  ddd()[,"names"],  organism=KEGGorg(), pvalueCutoff = 0.05, keyType = KEGGterms())
      kk1
    } else {
      
      kk1 <- enrichKEGG(gene = ddd()[,"names"] ,  organism=KEGGorg(), pvalueCutoff = 0.05, keyType = KEGGterms())
      kk1
    }
    
    
    return(kk1)
    
  })
  
  
  
  ## PLOTS FOR THE KEGG
  
  ## BarPlot ##
  
  output$KEGGbarplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- plot(barplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory)) # slider option
    }
    return(data)
  },  height = "auto",  width = "auto" , res = 72 )
  
  
  ###
  
  ## Goplot ##
  
  output$KEGGgoplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- plot(goplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory))
    }
    return(data)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  ## Dotplot ##
  
  output$KEGGdotplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- dotplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory)
      
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  ## cnetplot ##
  
  output$KEGGcnetplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGGcnetplot1_Shannon()
    }
    return(data)
  }, height = "auto",  width = "auto" , res = 72 )
  
  
  
  
  ###
  
  KEGGcnetplot1_Shannon <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory, color.params = list(foldChange=fold_change_geneListDeseq2_KEGG_Shannon())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory, color.params = list(foldChange=fold_change_geneList_KEGG_Shannon())) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  output$KEGGcnetplot21_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGG_Shannoncnetplot2()
    }
    return(data)
    
  }, height = "auto", width = "auto", res = 72 )
  
  
  
  KEGG_Shannoncnetplot2 <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <- cnetplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneListDeseq2_KEGG_Shannon(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- cnetplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneList_KEGG_Shannon(), circular = TRUE, colorEdge = TRUE) + scale_color_gradientn(name = "Log2 Fold Change", colours = c("#64AC59", "#98D7B7","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  ## upsetplot ##
  
  output$KEGGupsetplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- upsetplot(EnrichKeggShannon(), n = input$ShowKEGGCategory)
      
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72)
  
  ###
  
  
  ## heatplot2 ##
  
  
  output$KEGGheatplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      data <- NULL
      
    } else {
      data <- KEGG_ShannonHeatplot2()
    }
    return(data)
    
  },  height = "auto",  width = "auto" , res = 72)
  
  
  
  KEGG_ShannonHeatplot2 <- reactive({
    
    infile <- input$NormKegg
    req(infile)
    
    if (infile == "Deseq2") {
      
      fit3 <-  heatplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory, foldChange=fold_change_geneListDeseq2_KEGG_Shannon()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    } else {
      
      fit3 <- heatplot(EnrichKeggShannon(), showCategory = input$ShowKEGGCategory, foldChange= fold_change_geneList_KEGG_Shannon()) + scale_fill_gradientn(name = "Log2 Fold Change", colours = c("#BE2641", "#F5B35E","#325F8C"))
      fit3
      
    }
    
    return(fit3)
    
  })
  
  
  
  
  
  
  KEGG_ora_analysis_bp_Shannon <- reactive ({
    data <- pairwise_termsim(EnrichKeggShannon())
    data
  })
  
  
  
  ## emapplot ##
  
  output$KEGGemapplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      fit3 <- emapplot(KEGG_ora_analysis_bp_Shannon(), showCategory = input$ShowKEGGCategory)
      
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72 )
  
  ###
  
  
  
  
  ## treeplot ##
  
  output$KEGGtreeplot_Shannon <- renderPlot({
    infile <- input$KEGGkeytype
    req(infile)
    
    
    if (infile == 'Select the database') {
      fit3 <- NULL
      
    } else {
      
      fit3 <- treeplot(KEGG_ora_analysis_bp_Shannon(), showCategory = input$ShowKEGGCategory)
    }
    return(fit3)
  },  height = "auto",  width = "auto" , res = 72  )
  
  
  
  
  
  
  
  
  
  ###
  
  
  
  
  
  
  ### MDS plots output ####
  
  
  
  
  DataForMDS <- reactive({
    
    infile <- input$NormalizationMethod
    req(infile)
    
    
    
    if (infile == "Deseq2") {
      
      
      Data <- MDS_Deseq2(dds(), group())
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
    } else if (infile == "RLE") {
      
      
      Data <-  MDS_RLE(Metadata1(),Informationtable1(), design()) 
      
    } else if (infile == "TC") {
      
      Data <- MDS_TC(Metadata1(),Informationtable1(), design()) 
      
    } else if (infile == "Upper Quartile") {
      
      
      Data<-  MDS_UpperQuartile(Metadata1(),Informationtable1(), design()) 
      
      
    } else if (infile == "TMM") {
      
      Data <-  MDS_TMM(Metadata1(),Informationtable1(), design()) 
      
    } else if (infile == "TMMwsp") {
      
      
      Data <- MDS_TMMwsp(Metadata1(),Informationtable1(), design()) 
      
      
    } else if (infile == "PoissonSeq") {
      
      Data<- MDS_Poisson(Metadata1(),Informationtable1(), design()) 
      
    } else if (infile == "SVA") {
      
      Data <- MDS_SVA(Metadata1(),Informationtable1(), design()) 
      
      
    } else if (infile == "Median") {
      
      Data <-  MDS_Median(Metadata1(),Informationtable1(), design()) 
      
      
    } else if (infile == "Quantile") {
      
      
      Data<- MDS_Quantile(Metadata1(),Informationtable1(), design()) 
      
    } else {
      
      Data <- MDS_NN(Metadata1(),Informationtable1(), design()) 
    }
    
    return(Data)
  })
  
  
  
  output$DataMDS1 <-  renderPlot({
    DataForMDS()
    
  })
  
  
  
  ### Glimma Plots MDS ####
  
  
  DataForGlimmaMDS <- reactive({
    
    infile <- input$NormalizationMethod
    req(infile)
    
    if (infile == "Deseq2") {
      
      DataForGlimmaMDS <- deseqGlimmaMDS() 
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
    } else if (infile == "RLE") {
      
      DataForGlimmaMDS <- results_GlimmaMDSRLE()
      
    } else if (infile == "Upper Quartile") {
      
      
      DataForGlimmaMDS <- results_GlimmaMDSUQ()
      
      
    } else if (infile == "TMM") {
      
      DataForGlimmaMDS <- results_GlimmaMDSTMM()
      
    } else if (infile == "TMMwsp") {
      
      
      DataForGlimmaMDS <- results_GlimmaMDSTMMwsp()
      
      
    } else if (infile == "PoissonSeq") {
      
      DataForGlimmaMDS <- results_GlimmaMDSPoissonSeq()
      
    } else if (infile == "SVA") {
      
      DataForGlimmaMDS <- results_GlimmaMDSSVA()
      
      
    } else if (infile == "Median") {
      
      DataForGlimmaMDS <- results_GlimmaMDSMedian()
      
      
    } else if (infile == "Quantile") {
      
      
      DataForGlimmaMDS <- results_GlimmaMDSQuantile()
      
    } else {
      # No Normalization
      DataForGlimmaMDS <- results_GlimmaMDSNN()
    }
    
    return(DataForGlimmaMDS)
  })
  
  
  
  
  
  GlimmaPlot <- eventReactive(input$glimma, {
    DataForGlimmaMDS()
  })
  
  
  
  output$DataGlimmaMDS <-  renderPlot({
    data <- GlimmaPlot()
    data
  })
  
  
  
  ### Outputs in the new table ####
  
  
  output$Normalization <-  renderUI ({
    data <- pickerInput(inputId = 'NormalizationMethod', 
                        label = 'Select the Normalization method', 
                        choices = list('Deseq2','RLE','TMM','TMMwsp','PoissonSeq', 'Median', 'Quantile', 'Upper Quartile', 'No Normalization'),
                        selected = 'No Normalization',
                        options = list(`actions-box` = TRUE),multiple = F)
    req(data)
    data
    
  })
  
  
  
  
  
  DataForVolcano <- reactive({
    
    infile <- input$NormalizationMethod
    req(infile)
    
    if (infile == "Deseq2") {

      Data <- plot_ly(data = resdeseq11(), 
                      x = ~log2FoldChange, 
                      y = ~(-log10(padj)), 
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
      
    } else if (infile == "RLE") {
      
      Data <- plot_ly(data = resRLE11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),  
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
      
    } else if (infile == "Upper Quartile") {
      
      Data <- plot_ly(data = resUQ11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),  
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
      
      
    } else if (infile == "TMM") {
      
      Data <- plot_ly(data = resTMM11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),  
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
    } else if (infile == "TMMwsp") {
      
      Data <- plot_ly(data = resTMMwsp11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),  
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
      
    } else if (infile == "PoissonSeq") {
      
      Data <- plot_ly(data = resPoissonSeq11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),  
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)

      
    } else if (infile == "Median") {
      
      Data <- plot_ly(data = resMedian11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),   
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
      
    } else if (infile == "Quantile") {

      Data <- plot_ly(data = resQuantile11(), 
                      x = ~logFC,
                      y = ~(-log10(FDR)),  
                      color = ~regulated, 
                      colors = c("down" = "#9370DB", 
                                 "Not differentially expressed" = "#E0EEEE", 
                                 "up" = "#FFA54F"),
                      text = ~RefSeq, 
                      type = 'scatter', 
                      mode = 'markers',
                      marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
      
    } else {
      # No Normalization
      Data <- plot_ly(data = resdeseq11(), 
                      x = ~log2FoldChange, 
                      y = ~(-log10(padj)), 
                              color = ~regulated, 
                              colors = c("down" = "#9370DB", 
                                         "Not differentially expressed" = "#E0EEEE", 
                                         "up" = "#FFA54F"),
                              text = ~RefSeq, 
                              type = 'scatter', 
                              mode = 'markers',
                              marker = list(size = 7)) %>%
        layout(title = "Volcano Plot",
               xaxis = list(title = "log2 Fold Change"),
               yaxis = list(title = "-log10(padj)"),
               showlegend = TRUE)
    }
    
    return(Data)
  })
  
  
  
  
  
  
  output$Volcano <-  renderPlotly({ 
    data <- DataForVolcano()
    data
  })
  
  
  
  
  ### Volcano Plot Home Page ####
  
  output$VolcanoTMM <-  renderPlot({
    
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      d <- ggplot(data=resTMMb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
    }
    return(data)
  })
  
  
  output$VolcanoDeseq2<-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resdeseqb(), aes(x=log2FoldChange, y=-log10(padj), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d + ylab("-log10(FDR)")
    }
    return(data)
  })
  
  output$VolcanoRLE <-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resRLEb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  
  output$VolcanoTMMwsp<-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resTMMwspb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  
  output$VolcanoUpper <-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resUQb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  output$VolcanoQuantile <-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resQuantileb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  output$VolcanoPoisson <-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resPoissonSeqb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  
  
  output$VolcanoMedian <-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resMedianb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  
  output$VolcanoNN <-  renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      d <- ggplot(data=resNNb(), aes(x=logFC, y=-log10(FDR), color = regulated, text = RefSeq)) + geom_point() + theme_minimal() + scale_color_manual(values = c("down" ="#9370DB", "Not differentially expressed"= "#E0EEEE", "up" = "#FFA54F")) + theme(legend.position = "none") 
      data <- d +  xlab("log2FoldChange")
      data
    }
    return(data)
  })
  
  
  ####
  
  DataResultsPage <- reactive({
    
    
    infile <- input$NormalizationMethod
    req(infile)
    
    
    if (infile == "Deseq2") {
      DataResultsPage <-Deseq2v3()
      
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
      
    } else if (infile == "RLE") {
      
      DataResultsPage <- RLEv3()[,-6]
      
      
    } else if (infile == "Upper Quartile") {
      
      DataResultsPage <- UQv3()[,-6]
      
      
      
    } else if (infile == "TMM") {
      
      DataResultsPage <- TMMv3()[,-6]
      
    } else if (infile == "TMMwsp") {
      
      DataResultsPage <- TMMwspv3()[,-6]
      
      
    } else if (infile == "PoissonSeq") {
      
      DataResultsPage <- PoissonSeqv3()[,-6]
      
      
    } else if (infile == "Median") {
      
      DataResultsPage <- Medianv3()[,-6]
      
      
    } else if (infile == "Quantile") {
      
      DataResultsPage <- Quantilev3()[,-6]
      
      
    } else {
      DataResultsPage <- NNv3()[,-6]
      
    }
    
    
  })
  
  
  
  
  output$Tables<- DT::renderDataTable({
    df <- DataResultsPage()
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  output$ResultsNameNorm <- renderText({ 
    data <- input$NormalizationMethod
    
    f <- paste('Current normalization method: ', data)
    f
    
  })
  
  
  ### Up Down regulated genes table ####
  
  
  
  DataForUpDown <- reactive({
    
    infile <- input$NormalizationMethod
    req(infile)
    
    
    
    if (infile == "Deseq2") {
      
      
      UpDownData <- cbind(TableUpDownMedian(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownTMMwsp(),TableUpDownUQ() )
      Data <- rbind(UpDownData,GenesUpDownDeseq2())
      
    } else if (infile == "Select the Normalization method") {  
      
      return(NULL)
      
    } else if (infile == "RLE") {
      
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownTMM(),TableUpDownTMMwsp(),TableUpDownUQ() )
      Data <- rbind(UpDownData, GenesUpDownRLE())
      
    } else if (infile == "Upper Quartile") {
      
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownTMMwsp())
      Data <- rbind(UpDownData,GenesUpDownUQ())
      
    } else if (infile == "TMM") {
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMMwsp(),TableUpDownUQ() )  
      Data <- rbind(UpDownData,GenesUpDownTMM())
      
    } else if (infile == "TMMwsp") {
      
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownUQ() )
      Data <-rbind(UpDownData,GenesUpDownTMMwsp())
      
    } else if (infile == "PoissonSeq") {
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownNN(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownTMMwsp(),TableUpDownUQ() )
      Data <- rbind(UpDownData,GenesUpDownPoissonSeq())
      
      
    } else if (infile == "Median") {
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownTMMwsp(),TableUpDownUQ() )
      Data <- rbind(UpDownData,GenesUpDownMedian())
      
    } else if (infile == "Quantile") {
      
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownNN(),TableUpDownPoissonSeq(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownTMMwsp(),TableUpDownUQ() )
      Data <- rbind(UpDownData,GenesUpDownQuantile())
      
    } else {
      ## No normalization
      UpDownData <- cbind(TableUpDownDeseq(),TableUpDownMedian(),TableUpDownPoissonSeq(),TableUpDownQuantile(),TableUpDownRLE(),TableUpDownTMM(),TableUpDownTMMwsp(),TableUpDownUQ() )
      Data <- rbind(UpDownData,GenesUpDownNN())
    }
    
    return(Data)
  })
  
  
  
  
  
  output$DataUpDown <- DT::renderDataTable({
    
    infile <- input$SelectContr
    
    req(infile)
    
    if (is.null(infile)) {
      df2 <- NULL
      
    } else {
      df2 <- DataForUpDown()
    }
    
    DT::datatable(df2)
    
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  #
  
  
  
  
  ### Shannon Values ####
  
  output$MatrizE <-  renderUI ({
    data <- pickerInput(inputId = 'Ei', 
                        label = 'Select the Normalization Methods of interest', 
                        choices = list('Median', 'No Normalization','PoissonSeq', 'Quantile','RLE','TMM','TMMwsp', 'Upper Quartile'),
                        selected = c('Median','PoissonSeq','RLE','TMM','Upper Quartile'),
                        options = list(`actions-box` = TRUE),multiple = TRUE)
    req(data)
    data
    
  })
  
  
  MatrizInicialE <- reactive({
    
    
    infile <- input$Ei
    req(infile)
    MatrizInicialE <- NULL
    
    
    if ('Median' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_Median11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
      
    } 
    
    
    if ('No Normalization' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_NN11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
      
    }
    
    
    
    if ('PoissonSeq' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_PoissonSeq11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
      
    }
    
    if ('Quantile' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_Quantile11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
    } 
    
    
    
    if ( 'RLE' %in% infile) {
      MatrizInicialE <- cbind(MatrizInicialE,results_RLE11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
      
    }
    
    if ('TMM' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_TMM11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_TMM11())
    }
    
    if ('TMMwsp' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_TMMwsp11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
      
    }  
    
    
    if ('Upper Quartile' %in% infile) {
      
      MatrizInicialE <- cbind(MatrizInicialE,results_UpperQuartile11()[,"FDR"])
      row.names(MatrizInicialE) <- row.names(results_RLE11())
      
    }  
    
    
    
    return(MatrizInicialE)
  })
  
  
  ## First table
  
  output$MEE<- DT::renderDataTable({
    data <- MatrizInicialE() 
    colnames(data) <- c(input$Ei)
    data
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  #### Calculate the Shannon weigths and Ranks ####
  
  TakeTheDeletedToNA <- reactive({ 
    data <- MatrizInicialE()
    data <- as.matrix(data)
    data[data=="Deleted"] <- "NA"
    data <- as.numeric(data)
    data <- matrix(data,nrow=length(row.names(results_RLE11())))
    data
  })
  
  
  Ejl <- reactive({ apply(TakeTheDeletedToNA(), 2, function(col) col / sum(col, na.rm = T) ) })
  
  e0 <- reactive({ (log2(nrow(TakeTheDeletedToNA())))^-1 })
  
  m1 <- reactive({ 
    data <- Ejl()*log2(Ejl())
    data <- as.matrix(data)
    data
    
  })
  
  el <- reactive({ 
    data <- -e0() * colSums(m1(), na.rm = T)
    data <- matrix(data)
    data
  })
  
  dl <- reactive({ 
    data <- 1 - el() 
    data <- matrix(data)
    data
  })
  
  wl <- reactive ({ 
    data <- dl()/sum(dl(), na.rm = T) 
    data <- matrix(data)
    data
  })
  
  
  
  output$weigths<- DT::renderDataTable({
    
    infile <- MatrizInicialE()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      data <- wl() 
      data <- matrix(wl(),nrow = 1)
      colnames(data) <- c(input$Ei)
      data
    }
    return(data)
  })
  
  ## the ranks
  
  Bj <- reactive ({
    
    
    data <- rowSums(c(wl())*Ejl())
    data <- matrix(as.numeric(data), ncol = 1)
    data 
    
  })
  
  Bj1 <- reactive({
    
    data <- cbind(Bj(), c(row.names(results_RLE11())))
    data <- matrix(data, ncol = 2)
    dimnames(data) <- list(c(row.names(results_RLE11())),  c("Bj","Genes"))
    data
  })
  
  Bj2 <- reactive ({ 
    data <- matrix(as.numeric(Bj1()), ncol = 2)
    dimnames(data) <- list(c(row.names(results_RLE11())),  c("Bj","Genes"))
    data
  })
  
  Bj25 <- reactive({
    d <- Bj2()[order(Bj2()[,"Bj"], decreasing=FALSE),, drop=FALSE]
    f <- d[,"Bj"]
    g <- row.names(Bj2()[order(Bj2()[,"Bj"], decreasing=FALSE),, drop=FALSE])
    data <- Bj2()[order(Bj2()[,"Bj"], decreasing=FALSE),, drop=FALSE]
    data <- matrix(as.numeric(data), ncol = 2)
    dimnames(data) <- list(c(g), c("Bj","Genes"))
    data
  })
  
  Bj3 <- reactive({
    data <- cbind(Bj25(), c(1:nrow(Bj2())))
    data <- matrix(data,ncol =3)
    dimnames(data) <- list(c(row.names(Bj25())), c("Bj","Genes","rank"))
    data
  })
  
  
  Bj4 <- reactive({
    data <- Bj3()[order(row.names(Bj3())),, drop=FALSE]
    data <- matrix(data, ncol = 3)
    dimnames(data) <- list(c(Bj1()[,"Genes"]), c("Bj","Genes","rank"))
    data
    
  })
  
  
  Bj5 <- reactive({
    data <- Bj4()[,-c(2)]
    data <- matrix(data, ncol = 2)
    dimnames(data) <- list(c(row.names(results_RLE11())), c("Bj","rank"))
    data
    
  })
  
  
  
  
  MatrizInicialRanks<- reactive({
    # MUITO IMPORTANTE / Very important
    # I have to have the if's in the same order of the UI list ! It really is needed, else some values will be modified randomly
    
    infile <- input$Ei
    req(infile)
    MatrizInicialRanks <- NULL
    
    if ('Median' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_Median11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
      
    } 
    
    if ('No Normalization' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_NN11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
      
    }
    
    if ('PoissonSeq' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_PoissonSeq11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
      
    } 
    
    if ('Quantile' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_Quantile11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
    } 
    
    if ('RLE' %in% infile) {
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_RLE11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
      
    }
    
    
    if ('TMM' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_TMM11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
    }
    
    
    if ('TMMwsp' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_TMMwsp11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
      
    }  
    
    
    if ('Upper Quartile' %in% infile) {
      
      MatrizInicialRanks <- cbind(MatrizInicialRanks,results_UpperQuartile11()[,"rank"])
      row.names(MatrizInicialRanks) <- row.names(results_RLE11())
      
    }  
    
    
    return(MatrizInicialRanks)
  })
  
  
  
  ShannonRanksTable <- reactive ({
    
    data <- cbind(MatrizInicialRanks(),Bj5()[,"rank"])
    table <- matrix(data,nrow = length(row.names(results_RLE11())))
    rownames(table) <- rownames(results_RLE11())
    colnames(table) <- c(input$Ei,"SEA")
    table
    
  })
  
  
  
  output$ShannonRT<- DT::renderDataTable({
    
    infile <- MatrizInicialRanks()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- Bj5()
    }
    return(data)
    
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  output$FinalShannonRT<- DT::renderDataTable({
    
    
    infile <- MatrizInicialRanks()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- ShannonRanksTable()
      
    }
    return(data)
    
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  
  
  
  
  # Ranks Genes plot ####
  
  
  PlotRanks<- reactive({
    
    
    infile <- input$Ei
    req(infile)
    PlotRanks <- NULL
    
    
    if ('Median' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_Median11()[,c("rank", "Method", "Gene")])
      
    } 
    
    
    if ('No Normalization' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_NN11()[,c("rank", "Method", "Gene")])
      
      
    }
    
    
    if ('PoissonSeq' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_PoissonSeq11()[,c("rank", "Method", "Gene")])
      
      
    }
    
    
    if ('Quantile' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_Quantile11()[,c("rank", "Method", "Gene")])
      
    } 
    
    
    if ( 'RLE' %in% infile) {
      PlotRanks <- rbind(PlotRanks,results_RLE11()[,c("rank", "Method", "Gene")])
      
    }
    
    
    
    if ('TMM' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_TMM11()[,c("rank", "Method", "Gene")])
      
    }
    
    if ('TMMwsp' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_TMMwsp11()[,c("rank", "Method", "Gene")])
      
      
    }  
    
    
    if ('Upper Quartile' %in% infile) {
      
      PlotRanks <- rbind(PlotRanks,results_UpperQuartile11()[,c("rank", "Method", "Gene")])
      
      
    }  
    
    
    return(PlotRanks)
  })
  
  
  
  ShannonRanksT <- reactive ({
    
    infile <- input$Ei
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- PlotRanks()
      row.names(data) <- NULL
      data
    }
    
    return(data)
  })
  
  
  
  
  
  ## order by rank
  
  
  Bump2 <- reactive ({
    
    data <-  PlotRanks()[order(PlotRanks()$rank), ]
    data
  })
  
  
  
  #### Compare the Top Shannnon Genes with the others list ####
  
  
  
  results_Shannon <- reactive ({
    data <- ShannonRanksTable()[,"SEA"]
    data <-as.data.frame(ShannonRanksTable()[,"SEA"])
    data$rank <- c(ShannonRanksTable()[,"SEA"])
    data$Method <- c(rep("SEA", times = nrow(ShannonRanksTable())))
    data$Gene <- c(row.names(ShannonRanksTable()))
    data
  })
  
  
  tentativa <- reactive({
    data <-  results_Shannon()[order(results_Shannon()$rank), ]
    data
  })
  
  tentativa1 <- reactive({
    data <- tentativa()[1:input$integer,"Gene"]
    data <- as.data.frame(tentativa()[1:input$integer,"Gene", drop = FALSE])
    colnames(data) <- c("Gene")
    data
  })
  
  
  
  
  tentativa3 <- reactive({
    data <- as.data.frame(Bump2()[1:((input$integer)*(length(c(input$Ei)))),"Gene", drop = FALSE])
    colnames(data) <- c("Gene")
    data 
  })
  
  tentativa5 <- reactive({
    data <- tentativa3()
    data <- tentativa3() %>% 
      group_by(Gene) %>%
      summarise(n=n(), score = n()/length(c(input$Ei)))
    data
  })
  
  tentativa7 <- reactive({
    data <- as.data.frame(tentativa5())
    data
  })
  
  
  tentativa4 <- reactive({
    vec1 <- c(tentativa7()[,"Gene"])
    vec2 <- c(tentativa1()[,"Gene"])
    data <- tentativa5()[tentativa5()$Gene %in% vec2, ]
    data
  })
  
  tentativa6 <- reactive({
    data <- tentativa4()
    order_idx <- match(tentativa1()$Gene, data$Gene)
    
    valid_idx <- !is.na(order_idx)
    data <- data[order_idx[valid_idx], , drop = FALSE]
    
    # Optionally, fill in any missing genes from tentativa1() if needed
    if (any(!valid_idx)) {
      missing_genes <- tentativa1()$Gene[!valid_idx]
      missing_data <- data.frame(Gene = missing_genes, n = 0, score = 0)
      data <- rbind(data, missing_data)
    }
    
    data
    
    
  })
  
  
  output$tenta<- DT::renderDataTable({
    
    infile <- ShannonRanksTable()
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
      
      data <- tentativa6()
      row.names(data) <- NULL
      data
    }
    return(data)
    
  }, 
  # Aqui estou a adicionar alguns pormenores extras na tabela
  options = list(
    #lengthChange = FALSE,
    #scrollX = TRUE,
    #pageLength = 15,
    dom = 'Blfrtip',
    # como por exemplo fazer download da tabela nos seguintes ficheiros:
    buttons = c('copy','csv','excel','pdf','print')),
  autoHideNavigation = FALSE,
  filter = "top",
  extensions = 'Buttons',
  rownames = TRUE,
  server = TRUE)
  
  
  
  
  ### Venn Diagram ####
  
  
  
  output$VennOptions1 <-  renderUI ({
    data <- pickerInput(inputId = 'VennOptions1', 
                        label = 'Select Methods', 
                        choices = list('DESeq2','Median', 'No Normalization','PoissonSeq', 'Quantile','RLE','TMM','TMMwsp', 'Upper Quartile'),
                        selected = c('DESeq2','No Normalization','TMM'),
                        options = list(`actions-box` = FALSE, "max-options" = 7), multiple = T)
    req(data)
    data
    
  })
  
  PlotVenn <- reactive({
    
    infile <- input$VennOptions1
    req(infile)
    PlotVenn <- list()
    
    if ('DESeq2' %in% infile) {
      PlotVenn$DESeq2 <-   Venn_Deseq2()
    } 
    
    
    if ('Median' %in% infile) {
      PlotVenn$Median <- Venn_Median()
    } 
    
    if ('No Normalization' %in% infile) {
      PlotVenn$NN <- Venn_NN()
    }
    
    if ('PoissonSeq' %in% infile) {
      PlotVenn$PoissonSeq <- Venn_PoissonSeq()
    }
    
    
    if ('Quantile' %in% infile) {
      PlotVenn$Quantile <- Venn_Quantile()
    } 
    
    if ('RLE' %in% infile) {
      PlotVenn$RLE <- Venn_RLE()
    }
    
    if ('TMM' %in% infile) {
      PlotVenn$TMM <- Venn_TMM()
    }
    
    if ('TMMwsp' %in% infile) {
      PlotVenn$TMMwsp <- Venn_TMMwsp()
    }  
    
    if ('Upper Quartile' %in% infile) {
      PlotVenn$UQ <- Venn_UQ()
    }  
    
    
    return(PlotVenn)
  })
  
  
  
  output$Venn <- renderPlot({
    
    infile <- input$cont
    req(infile)
    
    if (is.null(infile)) {
      data <- NULL
      
    } else {
    
    p <- ggVennDiagram(PlotVenn(), label = "count")
    
    data <- p +  scale_x_continuous(expand = expansion(mult = .1)) + scale_fill_distiller(palette = "RdBu") 
    data
    
  }
    return(data)}
  )
  
  
  
  
  output$VennOptions <-  renderUI ({
    data <- pickerInput(inputId = 'VennOptions', 
                        label = 'Select Method', 
                        choices = list('Select','Median', 'No Normalization','PoissonSeq', 'Quantile','RLE','TMM','TMMwsp', 'Upper Quartile'),
                        selected = 'TMM',
                        options = list(`actions-box` = TRUE),multiple = F)
    req(data)
    data
    
  })
  
  
  
  PlotVenn2 <- reactive({
    
    infile <- input$VennOptions
    req(infile)
    PlotVenn2 <- list()
    
    if ('Select' %in% infile) {
      
      PlotVenn2 <- NULL
      
    } else if ('Median' %in% infile) {
      
      
      PlotVenn2$Median <- Venn_Median()
      
      PlotVenn2$SEA <- c(tentativa()[1:length(row.names(datMedian())),"Gene"])
      
    } else if ('No Normalization' %in% infile) {
      PlotVenn2$NN <- Venn_NN()
      
      vetor <- c(row.names(tentativa()[1:length(row.names(datasetNN())),]))
      PlotVenn2$SEA <- c(row.names(results_NN()[vetor,]))
      
      
    } else if ('PoissonSeq' %in% infile) {
      PlotVenn2$PoissonSeq <- Venn_PoissonSeq()
      
      vetor <- c(row.names(tentativa()[1:length(row.names(datasetPoissonSeq())),]))
      PlotVenn2$SEA  <- c(row.names(results_PoissonSeq()[vetor,]))
      
      
    }else if ('Quantile' %in% infile) {
      PlotVenn2$Quantile <- Venn_Quantile()
      
      vetor <- c(row.names(tentativa()[1:length(row.names(datasetQuantile())),]))
      PlotVenn2$SEA <- c(row.names(results_Quantile()[vetor,]))
      
      
    } else if ('RLE' %in% infile) {
      PlotVenn2$RLE <- Venn_RLE()
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(datasetRLE()))])
      PlotVenn2$SEA  <- c(row.names(results_RLE()[vetor,]))
      
      
    } else if ('TMM' %in% infile) {
      PlotVenn2$TMM <- Venn_TMM()
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(datasetTMM()))])
      PlotVenn2$SEA  <- c(row.names(results_TMM()[vetor,]))
      
      
    } else if ('TMMwsp' %in% infile) {
      PlotVenn2$TMMwsp <- Venn_TMMwsp()
      
      data <- tentativa()[,"Gene"]
      vetor <- c(data[1:length(row.names(datasetTMMwsp()))])
      PlotVenn2$SEA <- c(row.names(results_TMMwsp()[vetor,]))
      
      
    }  else  {
      PlotVenn2$UQ <- Venn_UQ()
      
      
      vetor <- c(row.names(tentativa()[1:length(row.names(datasetUpperQuartile())),]))
      PlotVenn2$SEA <- c(row.names(results_UpperQuartile()[vetor,]))
      
      
    }  
    
    
    
    return(PlotVenn2)
  })
  
  
  output$Venn2 <- renderPlot({
    
    infile <- input$VennOptions
    req(infile)
    
    if ('Select' %in% infile) {
      
      f <- NULL
      
    } else {
      
      data <- ggVennDiagram(PlotVenn2(), label = "count")
      
      #data
      
      f <- ggvenn(PlotVenn2(), show_elements = F, label_sep = "\n", fill_color = c("#AE123A","#1C65A3"), text_size = 6, count_column = T)
      
      f }
    return(f)
  })
  
  
  

  
  ##### Logo Image of GeneSEA Explorer #####
  
  # add output in server
  output$home_img <- renderImage({
    
    list(src = "official_banner.png",
         width = "100%", height = "80%")
    
  }, deleteFile = F)
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)