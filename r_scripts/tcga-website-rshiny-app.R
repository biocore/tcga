# RShiny_Microbe_DA_app.R
# Gregory Poore
# Aug 31, 2019
# Purpose: Create an interactive website for microbial abundances

# # Load dependencies
# require(devtools)
# require(doMC)
# require(dplyr)
# require(reshape2)
# require(ggplot2)
# require(ggrepel)
# require(ggsci)
# require(pheatmap)
# require(plotly)
# 
# numCores <- detectCores()
# registerDoMC(cores=numCores)
# 
# library(shiny)
# runExample("01_hello")

library(shiny)
library(shinyauthr)
library(shinyjs)
library(shinydashboard)
library(ggpubr)
library(ggsci)
library(dplyr)
library(glue)
library(DT)
library(markdown) # for displaying text files

## Load data
# Load Kraken data
load("data/mergedSnmDataAndMetadataQC.RData")

# Load Shogun data
load("./data/Shogun_All_Quantile_Nov202019_Final.RData")

# Load Plasma data - raw
load("./data/vbMergedDataKrakenDPCF_Dec3.RData")
load("./data/shogunDataAndMetadataPlasmaRAW_Dec3.RData")

# Load Plasma data - normalized
load("./data/snmKrakenAndMetadataFiltered_Dec2_Final.RData")
load("./data/snmCfdnaShogunAndMetadata_Dec2_Final.RData")

#------------------------------------------------------------#
# One-time operations
autocompleteNames <- tail(names(mergedSnmDataAndMetadataQC),-49)

abbreviationTable <- read.csv(file = "./data/abbreviationTable.csv",
                              header = TRUE,
                              strip.white = TRUE,
                              stringsAsFactors = FALSE)
colnames(abbreviationTable) <- c("TCGA Abbreviation", "Cancer Type")

dataTypeList <- c("Full Data (i.e. everything detected in TCGA)",
                  "Likely Contaminants Removed (i.e. filtered pathogens re-allowed)",
                  "Decontamination by Sequencing Plate & Center",
                  "All Putative Contaminants Removed",
                  "Most Stringent Decontamination")

shogunDataTypeList <- c("Shogun Normalized Data",
                        "Kraken-Matched-to-Shogun Normalized Data (i.e. no viruses, fewer samples)")

plasmaDataTypeList <- c("Raw Kraken Data",
                        "Raw Shogun Data",
                        "Decontaminated & Age-Sex Normalized Kraken Data",
                        "Decontaminated & Age-Sex Normalized Shogun Data")

selectedSampleTypeOrContrastList <- c("Primary Tumor",
                                      "Blood Derived Normal",
                                      "Primary Tumor vs Solid Tissue Normal",
                                      "Stage I vs IV")

shogunAutocompleteNames <- colnames(snmDataShogunQCFilt)
shogunMergedSnmDataAndMetadata <- cbind(shogunMetadataQCFilt, snmDataShogunQCFilt)

metadataPSMatchedDP_RAW <- metadataPSMatchedDPShogun 
# ^ It can be the Kraken version if desired. They're equiv. but in different order.

metadataPSMatchedDP_RAWQC <- droplevels(metadataPSMatchedDP_RAW[!(metadataPSMatchedDP_RAW$hiv_status_clean == "HIV+"),])

metadataPSMatchedDP_RAWQC$disease_type_consol <- metadataPSMatchedDP_RAWQC$baseline_plasma_and_serum
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type %in% c("Lung Adenocarcinoma",
                                                                                           "Lung Squamous Cell Carcinoma",
                                                                                           "nsclc NOS",
                                                                                           "NSCLC Sarcomatoid") |
                                                metadataPSMatchedDP_RAWQC$diagnosis == "NSCLC- Large cell"] <- "Lung Cancer"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type %in% c("Skin Cutaneous Melanoma") |
                                                metadataPSMatchedDP_RAWQC$diagnosis %in% c("Melanoma of RU chest wall",
                                                                                   "Melanoma, unknown primary")] <- "Melanoma"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type %in% c("Prostate Cancer")] <- "Prostate Cancer"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type %in% c("healthy control")] <- "Healthy Control"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type_consol %in% c("control empty well")] <- "Empty Well Control"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type_consol %in% c("bacteria monoculture")] <- "Bacteria Monoculture"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type_consol %in% c("control blank DNA extraction")] <- "Blank (DNA Extraction) Control"
metadataPSMatchedDP_RAWQC$disease_type_consol[metadataPSMatchedDP_RAWQC$disease_type_consol %in% c("control blank library prep")] <- "Blank (Library Prep) Control"

metadataPSMatchedDP_RAWQC$disease_type_consol <- factor(metadataPSMatchedDP_RAWQC$disease_type_consol)
# dim(metadataPSMatchedDP_RAWQC) # 266 217

#------------------------------------------------------------#

user_base <- data.frame(
  user = c("knightlab"), # knightlab
  password = c("cancermicrobiome"), # cancermicrobiome
  permissions = c("standard"),
  name = c("Knight Lab"),
  stringsAsFactors = FALSE,
  row.names = NULL
)

ui <- dashboardPage(
  
  skin = "blue",
  
  dashboardHeader(title = "TCGA Cancer Microbiome",
                  titleWidth = 350,
                  tags$li(class = "dropdown",
                          tags$a(icon("database"), 
                                 href = "https://drive.google.com/drive/folders/18V2ON-Go5AeEtZLe1f9EeJToWOhg81ab?usp=sharing",
                                 title = "Link to download microbe data")),
                  tags$li(class = "dropdown",
                          tags$a(icon("github"), 
                                 href = "https://github.com/biocore/tcga",
                                 title = "See the code repository on Github"),
                  tags$li(class = "dropdown", style = "padding: 8px;",
                          shinyauthr::logoutUI("logout"))
                  )
  ),
  
  dashboardSidebar(width = 350,
                   collapsed = TRUE, 
                   div(textOutput("welcome"), 
                       style = "padding: 5px",
                       imageOutput("klLogo", width = "100%", inline = TRUE),
                       # id = "sidebar", # id important for updateTabItems
                       sidebarMenu(
                       menuItem("Kraken microbe abundances in TCGA cancer types", tabName = "home", icon = icon("chart-bar")),
                       menuItem("Shogun microbe abundances in TCGA cancer types", tabName = "tab2", icon = icon("chart-bar")),
                       menuItem("Kraken TCGA model performance and feature list", tabName = "tab1", icon = icon("table")),
                       menuItem("Shogun TCGA model performance and feature list", tabName = "tab3", icon = icon("table")),
                       menuItem("Plasma cell-free microbial abundances (validation)", tabName = "tab4", icon = icon("chart-bar"))),
                       tableOutput("abbrev"))
  ),
  
  dashboardBody(
    # tags$head(tags$style(HTML(' .main-sidebar{ width: 300px; } .main-header > .navbar { margin-left: 300px; } .main-header .logo { width: 300px; } .content-wrapper, .main-footer, .right-side { margin-left: 300px; } '))),
    # must turn shinyjs on
    shinyjs::useShinyjs(),
    # add logout button UI 
    div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
    # add login panel UI function
    shinyauthr::loginUI(id = "login"),
    
    tags$style(type="text/css",".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; }"),
    tabItems(
      tabItem("home", uiOutput("krakenMicrobialAbundances")),
      tabItem("tab1", uiOutput("krakenPerfUI")), 
      tabItem("tab2", uiOutput("shogunMicrobialAbundances")),
      tabItem("tab3", uiOutput("shogunPerfUI")),
      tabItem("tab4", uiOutput("plasmaMicrobialAbundances"))
    )
  )
)

# Define server logic 
server <- function(input, output,session) {
  
  # call the logout module with reactive trigger to hide/show
  logout_init <- callModule(shinyauthr::logout, 
                            id = "logout", 
                            active = reactive(credentials()$user_auth))
  
  # call login module supplying data frame, user and password cols
  # and reactive trigger
  credentials <- callModule(shinyauthr::login, 
                            id = "login", 
                            data = user_base,
                            user_col = user,
                            pwd_col = password,
                            log_out = reactive(logout_init()))
  
  observe({
    if(credentials()$user_auth) {
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    } else {
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
  })
  
  # pulls out the user information returned from login module
  user_data <- reactive({credentials()$info})
  
  user_info <- reactive({credentials()$info})
  
  output$welcome <- renderText({
    req(credentials()$user_auth)

    glue("Welcome {user_info()$name}")
  })
  
  output$klLogo <- renderImage({
    req(credentials()$user_auth)
    list(src = "./logos/kl-logo11.png",
         width = 350)
    }, deleteFile = FALSE)
  
  output$abbrev <- renderTable({
    req(credentials()$user_auth)
    abbreviationTable
    })
    
  output$krakenMicrobialAbundances <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      titlePanel("Kraken-Derived Normalized Microbial Abundances in TCGA"),
      
      helpText("Select microbe of interest to display its
               normalized abundance distribution across 
               TCGA cancer types."),
      
      p("Note: Kraken-derived data were generated using Kraken 1 against 59,974 quality filtered
        microbial genomes, originally downloaded from RepoPhlan on 14 June 2016 (original size was
        71,782 before quality filtering). This included bacteria, archaea, and viruses. This database is
        *different* than the 'Web of Life' database (accepted for publication; https://biocore.github.io/wol/) used for generating Shogun-derived data in TCGA. The Voom-SNM
        normalized data are plotted below.", style = "color:red"),
      
      selectizeInput("kraken_microbe_name", "Either (i) select from the drop-down list or
                    (ii) click on the name and hit backspace, then start typing the name of your microbe of interest
                     (it will autocomplete). After selecting your microbe of interest, please wait a second
                     for the plot to appear below. A basic Kruskal-Wallis test is performed for each sample type
                     to show if the microbe varies across cancer types. The genus for the HPV virus (Alphapapillomavirus) 
                     is pre-selected for you.",
                     choices = autocompleteNames,
                     selected = "k__Viruses.f__Papillomaviridae.g__Alphapapillomavirus",
                     width = '100%'),
      
      mainPanel(
        
      renderPlot(width = 1000, height = 750, expr = {
                   
                   microbeName <- input$kraken_microbe_name
                   mergedSnmDataAndMetadataQC$sample_type <- factor(mergedSnmDataAndMetadataQC$sample_type, 
                                                                    levels = c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))
                   mergedSnmDataAndMetadataQC %>%
                     select(investigation, sample_type, microbeName) %>%
                     filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
                     ggboxplot(x = "investigation", y = microbeName, 
                               color = "sample_type",
                               palette = "lancet",
                               # notch = TRUE,
                               facet.by = "sample_type",
                               nrow = 3,
                               xlab = "Disease Type", ylab = "Normalized Abundance (log2-cpm)", 
                               title = paste(microbeName, "\n\nAbundances Across TCGA Cancer Types"),
                               legend = "top",
                               legend.title = "Sample Type") +
                     rotate_x_text(angle = 30) +
                     theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
                     stat_compare_means(label.x = 5)
                   
                 })
        )
      )
    })
  
  output$krakenPerfUI <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      
      fluidRow(
        
        column(4,
               selectInput("dataTypeSelected", 
                           label = "Select Data Used for Model Training/Testing:", 
                           selected = "Full Data (i.e. everything detected in TCGA)",
                           choices = dataTypeList)
        ),
        
        column(4,
               selectInput("selectedCancerType",
                           label = "Select Cancer Type:", 
                           selected = "Adrenocortical Carcinoma",
                           choices = abbreviationTable$`Cancer Type`)
        ),
        
        column(4,
               selectInput("selectedSampleTypeOrContrast",
                           label = "Select Sample Type OR Comparison:",
                           selected = "Primary Tumor",
                           choices = selectedSampleTypeOrContrastList)
        )
        
      ),
      
      helpText("Color bars show probability cutoff thresholds. 
               Performance metrics are one-cancer-type-versus-all-others unless looking at
               tumor versus normal comparisons or tumor stage comparisons. A constant random number seed was used to ensure
               that the models are comparable in their performance. For training, 70% of the samples
               were randomly selected in a class-stratified manner; the remaining 30% of samples were 
               used as a holdout test set, on which performance metrics were generated. During training,
               4-fold cross validation procedures were employed to tune the model hyperparameters. Again,
               the performance metrics shown here (ROC and PR curves, confusion matrices) are based on the
               30% independent holdout test set for each comparison of interest."),
      
      strong(helpText("NOTE: If no plots appear, then the comparison was not modeled. This is
               due to having <20 samples in the minority class.")),
      
      
      fluidRow(class = "plotRow",
        
        column(6,
               renderImage({
                 validate(
                   need(file.exists(paste0("newimages/",
                                           input$dataTypeSelected,
                                           "/",
                                           input$selectedCancerType,
                                           " -- ",
                                           input$selectedSampleTypeOrContrast,
                                           " -- ROC.png")),
                        "Comparison not available. Nothing will be plotted."
                   ))
                 
                 list(src = paste0("newimages/",
                                   input$dataTypeSelected,
                                   "/",
                                   input$selectedCancerType,
                                   " -- ",
                                   input$selectedSampleTypeOrContrast,
                                   " -- ROC.png"),
                      width = 450)
               }, deleteFile = FALSE)
        ),
        
        column(6,
               renderImage({
                 validate(
                   need(file.exists(paste0("newimages/",
                                           input$dataTypeSelected,
                                           "/",
                                           input$selectedCancerType,
                                           " -- ",
                                           input$selectedSampleTypeOrContrast,
                                           " -- PR.png")),
                        "Comparison not available. Nothing will be plotted."
                   ))
                 
                 list(src = paste0("newimages/",
                                   input$dataTypeSelected,
                                   "/",
                                   input$selectedCancerType,
                                   " -- ",
                                   input$selectedSampleTypeOrContrast,
                                   " -- PR.png"),
                      width = 450)
               }, deleteFile = FALSE)
        ),
        tags$head(tags$style(".plotRow{margin-bottom:-400px;}"))
      ),
      
      fluidRow(
        column(12, 
               align = "center",
               strong("Confusion Matrix (using a 50% probability cutoff threshold on the 30% holdout test set):", style = "color:black"),
               p("Note: A 50% probability cutoff may *not* always be the best choice for class discrimination", style = "color:red"),
               pre(renderText({
                 
                 
                 
                 cmPath <- paste0("./newconfusionmatrices/",
                                  input$dataTypeSelected,
                                  "/",
                                  input$selectedCancerType,
                                  " -- ",
                                  input$selectedSampleTypeOrContrast,
                                  " -- CM.txt")
                 
                 includeText(cmPath)
                 }))
        )
      ),
      
      h3("Ranked list of model features:"),
      
      helpText("The Variable Importance Score is calculated based on the R gbm package."),
      
      p("Note: Feature importance scores are *not* directly compatible between models, as their calculation depends on 
        how the model was tuned during training. It is thus a *relative* feature importance score. Moreover, a high feature 
        importance score for a given taxon does *not* guarantee or imply an overabundance of that taxon for this comparison.
        A high feature importance score only means that the taxon was important for making predictions; whether the taxon is 
        less or more abundant in certain samples requires statistical testing. Lastly, the models shown in this paper were 
        *not* regularized, meaning that they could theoretically use all taxa features available in the data to make predictions. 
        All features with non-zero feature importance scores are shown here.", style = "color:red"),
      
      renderDataTable({
        
        df <- read.csv(paste0("./newfeatures/",
                              input$dataTypeSelected,
                              "/",
                              input$selectedCancerType,
                              " -- ",
                              input$selectedSampleTypeOrContrast,
                              " -- Features.csv"))
        colnames(df) <- c("Taxonomy", "Variable Importance Score")
        df
      })
      
      
      
      # fluidRow(
      #   column(6, renderImage("data/barnacle-likelyContam-pproc-PR-predicting-cancer-type-1-vs-all/Adrenocortical Carcinoma - Primary Tumor (CV k-fold of 4|Train proportion of 70) PRROC PR.png")),
      #   column(6, renderImage("data/barnacle-likelyContam-pproc-ROC-predicting-cancer-type-1-vs-all/Adrenocortical Carcinoma - Primary Tumor (CV k-fold of 4|Train proportion of 70) PRROC ROC.png"))
      # )
      
    )
    
  })
  
  output$shogunMicrobialAbundances <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      titlePanel("Shogun-Derived Normalized Microbial Abundances in TCGA"),
      
      helpText("Select microbe of interest to display its 
               normalized abundance distribution across 
               TCGA cancer types."),
      
      p("Note: Shogun-derived data were generated using SHallow shOtGUN (i.e. SHOGUN) alignments against the 'Web of Life' database (accepted
        for publication; https://biocore.github.io/wol/), which contains 10,575 bacterial and archaeal genomes. This database is *different* than
        the RepoPhlan database used for generating Kraken-derived data in TCGA, particularly in that it does *not* include viruses.
        Also, given the computational burden of doing alignment-based taxonomy assignments, the number of samples in the Shogun analysis is smaller than that in the
        Kraken-only data (n=13,517 vs 17,625); however, the total number of analyzed cancer types is the same (n=32) and the
        Voom-SNM pipeline was also used for normalization.", style = "color:red"),
      
      selectizeInput("shogun_microbe_name", "Either (i) select from the drop-down list or
                     (ii) click on the name and hit backspace, then start typing the name of your microbe of interest
                     (it will autocomplete). After selecting your microbe of interest, please wait a second
                     for the plot to appear below. A basic Kruskal-Wallis test is performed for each sample type
                     to show if the microbe varies across cancer types. The Fusobacterium genus is pre-selected for you.",
                     choices = shogunAutocompleteNames,
                     selected = "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Fusobacterium",
                     width = '100%'),
      
      mainPanel(
        
        renderPlot(width = 1000, height = 750, expr = {
          
          microbeName <- input$shogun_microbe_name
          # shogunMergedSnmDataAndMetadata$sample_type <- factor(mergedSnmDataAndMetadataQC$sample_type, 
          #                                                  levels = c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))
          shogunMergedSnmDataAndMetadata %>%
            select(investigation, sample_type, microbeName) %>%
            filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
            ggboxplot(x = "investigation", y = microbeName, 
                      color = "sample_type",
                      palette = "lancet",
                      # notch = TRUE,
                      facet.by = "sample_type",
                      nrow = 3,
                      xlab = "Disease Type", ylab = "Normalized Abundance (log2-cpm)", 
                      title = paste(microbeName, "\n\nAbundances Across TCGA Cancer Types"),
                      legend = "top",
                      legend.title = "Sample Type") +
            rotate_x_text(angle = 30) +
            theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
            stat_compare_means(label.x = 5)
          
        })
      )
      )
  })
  
  output$shogunPerfUI <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      
      p("Note: To have fair comparisons between models built on Kraken-derived data and Shogun-derived data,
        a version of the Kraken-derived raw data were subsetted with *no viruses* (as the 'Web of Life' database did not
        contain them; https://biocore.github.io/wol/) and were matched to only those TCGA samples processed by Shogun (n=13,517). Again, this is *less than* the total
        number of samples analyzed in the Kraken-only data (n=17,625). All other microbial assignments were kept in the Kraken-derived data. 
        Then, the Shogun-derived data and Kraken-derived subsetted data were normalized equally by Voom and SNM before inputting into machine learning pipelines.", style = "color:red"),
      
      fluidRow(
        
        column(4,
               selectInput("shogunDataTypeSelected", 
                           label = "Select Data Used for Model Training/Testing:", 
                           selected = "Shogun Normalized Data",
                           choices = shogunDataTypeList)
        ),
        
        column(4,
               selectInput("shogunSelectedCancerType",
                           label = "Select Cancer Type:", 
                           selected = "Adrenocortical Carcinoma",
                           choices = abbreviationTable$`Cancer Type`)
        ),
        
        column(4,
               selectInput("shogunSelectedSampleTypeOrContrast",
                           label = "Select Sample Type OR Comparison:",
                           selected = "Primary Tumor",
                           choices = selectedSampleTypeOrContrastList)
        )
        
      ),
      
      helpText("Color bars show probability cutoff thresholds. 
               Performance metrics are one-cancer-type-versus-all-others unless looking at
               tumor versus normal comparisons or tumor stage comparisons. A constant random number seed was used to ensure
               that the models are comparable in their performance. For training, 70% of the samples
               were randomly selected in a class-stratified manner; the remaining 30% of samples were 
               used as a holdout test set, on which performance metrics were generated. During training,
               4-fold cross validation procedures were employed to tune the model hyperparameters. Again,
               the performance metrics shown here (ROC and PR curves, confusion matrices) are based on the
               30% independent holdout test set for each comparison of interest."),
      
      strong(helpText("NOTE: If no plots appear, then the comparison was not modeled. This is
                      due to having <20 samples in the minority class.")),
      
      
      fluidRow(class = "plotRow",
               
               column(6,
                      renderImage({
                        validate(
                          need(file.exists(paste0("shogunimages/",
                                                  input$shogunDataTypeSelected,
                                                  "/",
                                                  input$shogunSelectedCancerType,
                                                  " -- ",
                                                  input$shogunSelectedSampleTypeOrContrast,
                                                  " -- ROC.png")),
                               "Comparison not available. Nothing will be plotted."
                          ))
                        
                        list(src = paste0("shogunimages/",
                                          input$shogunDataTypeSelected,
                                          "/",
                                          input$shogunSelectedCancerType,
                                          " -- ",
                                          input$shogunSelectedSampleTypeOrContrast,
                                          " -- ROC.png"),
                             width = 450)
                      }, deleteFile = FALSE)
               ),
               
               column(6,
                      renderImage({
                        validate(
                          need(file.exists(paste0("shogunimages/",
                                                  input$shogunDataTypeSelected,
                                                  "/",
                                                  input$shogunSelectedCancerType,
                                                  " -- ",
                                                  input$shogunSelectedSampleTypeOrContrast,
                                                  " -- PR.png")),
                               "Comparison not available. Nothing will be plotted."
                          ))
                        
                        list(src = paste0("shogunimages/",
                                          input$shogunDataTypeSelected,
                                          "/",
                                          input$shogunSelectedCancerType,
                                          " -- ",
                                          input$shogunSelectedSampleTypeOrContrast,
                                          " -- PR.png"),
                             width = 450)
                      }, deleteFile = FALSE)
               ),
               tags$head(tags$style(".plotRow{margin-bottom:-400px;}"))
      ),
      
      fluidRow(
        column(12, 
               align = "center",
               strong("Confusion Matrix (using a 50% probability cutoff threshold on the 30% holdout test set):", style = "color:black"),
               p("Note: A 50% probability cutoff may *not* always be the best choice for class discrimination", style = "color:red"),
               pre(renderText({
                 
                 
                 
                 cmPath <- paste0("./shogunconfusionmatrices/",
                                  input$shogunDataTypeSelected,
                                  "/",
                                  input$shogunSelectedCancerType,
                                  " -- ",
                                  input$shogunSelectedSampleTypeOrContrast,
                                  " -- CM.txt")
                 
                 includeText(cmPath)
               }))
        )
      ),
      
      h3("Ranked list of model features:"),
      
      helpText("The Variable Importance Score is calculated based on the R gbm package."),
      
      p("Note: Feature importance scores are *not* directly compatible between models, as their calculation depends on 
        how the model was tuned during training. It is thus a *relative* feature importance score. Moreover, a high feature 
        importance score for a given taxon does *not* guarantee or imply an overabundance of that taxon for this comparison.
        A high feature importance score only means that the taxon was important for making predictions; whether the taxon is 
        less or more abundant in certain samples requires statistical testing. Lastly, the models shown in this paper were 
        *not* regularized, meaning that they could theoretically use all taxa features available in the data to make predictions. 
        All features with non-zero feature importance scores are shown here.", style = "color:red"),
      
      renderDataTable({
        
        df <- read.csv(paste0("./shogunfeatures/",
                              input$shogunDataTypeSelected,
                              "/",
                              input$shogunSelectedCancerType,
                              " -- ",
                              input$shogunSelectedSampleTypeOrContrast,
                              " -- Features.csv"))
        colnames(df) <- c("Taxonomy", "Variable Importance Score")
        df
      })
      
      
      
      # fluidRow(
      #   column(6, renderImage("data/barnacle-likelyContam-pproc-PR-predicting-cancer-type-1-vs-all/Adrenocortical Carcinoma - Primary Tumor (CV k-fold of 4|Train proportion of 70) PRROC PR.png")),
      #   column(6, renderImage("data/barnacle-likelyContam-pproc-ROC-predicting-cancer-type-1-vs-all/Adrenocortical Carcinoma - Primary Tumor (CV k-fold of 4|Train proportion of 70) PRROC ROC.png"))
      # )
      
      )
    
  })
  
  output$plasmaMicrobialAbundances <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      titlePanel("Real-World Cell-free Microbial Abundances in Cancer and Non-Cancer Controls"),
      
      helpText("Select microbe of interest to display its 
               abundance distribution across 
               samples types."),
      
      p("Note: There are several options given below, but some are *not* directly comparable. Kraken-derived
        data were generated using Kraken 1 and the same RepoPhlan database used for TCGA (n=59,974 microbial genomes [
        bacteria, archaea, viruses]) while Shogun-derived data used the 'Web of Life' database (n=10,575 microbial genomes
        [bacteria, archaea]; https://biocore.github.io/wol/). Shogun data were additionally filtered at the 0.1% level, removing any microbes that did not
        pass this abundance threshold as a conservative measure (i.e. its data are sparser). Both raw datasets were independently decontaminated using 'prevalence'
        and 'frequency' methods of the decontam tool (https://github.com/benjjneb/decontam), using the stringent filtering
        hyperparameter of 0.5. The count data were then Voom-SNM normalized, analogous to that done for the TCGA data,
        before inputting into machine learning pipelines. The raw, pre-filtered, pre-decontaminated data for both Kraken
        and Shogun are shown below as well as their filtered, decontaminated, Voom-SNM normalized counterparts. Also note
        that technical controls (e.g. negative extraction blanks, positive monoculture wells) have been removed in the filtered, decontaminated, Voom-SNM normalized data.
        Lastly, note that some genera may not be available across datasets, either due to database differences or decontamination differences (done independently).", style = "color:red"),
      
      selectInput("plasmaDataTypeSelected", 
                  label = "Select Data Type:", 
                  selected = "Raw Kraken Data",
                  choices = plasmaDataTypeList),
      
      renderUI({
        
        if(input$plasmaDataTypeSelected == "Raw Kraken Data"){
          plasmaMetadata <- metadataPSMatchedDP_RAWQC
          plasmaData <- vbMergedDataKrakenDPCF[rownames(plasmaMetadata),]
          plasmaAutocompleteNames <- colnames(plasmaData)
        } else if(input$plasmaDataTypeSelected == "Raw Shogun Data"){
          plasmaMetadata <- metadataPSMatchedDP_RAWQC
          plasmaData <- dataPSUniqueDP[rownames(plasmaMetadata),]
          plasmaAutocompleteNames <- colnames(plasmaData)
        } else if(input$plasmaDataTypeSelected == "Decontaminated & Age-Sex Normalized Kraken Data"){
          plasmaData <- snmDataKrakenCFDecontamDPQC
          plasmaMetadata <- metadataPSMatchedDPQCFiltered
          plasmaAutocompleteNames <- colnames(plasmaData)
        } else if(input$plasmaDataTypeSelected == "Decontaminated & Age-Sex Normalized Shogun Data"){
          plasmaData <- snmShogundataPSUniqueDecontamDPQC
          plasmaMetadata <- metadataPSMatchedDPQCFiltered
          plasmaAutocompleteNames <- colnames(plasmaData)
        }
        
        selectizeInput("plasma_microbe_name", "Either (i) select from the drop-down list or
                     (ii) click on the name and hit backspace, then start typing the name of your microbe of interest
                       (it will autocomplete). After selecting your microbe of interest, please wait a second
                       for the plot to appear below. A basic Kruskal-Wallis test is performed
                       to show if the microbe varies across sample types. The Nocardia genus is pre-selected for you as
                       an example.",
                       choices = plasmaAutocompleteNames,
                       selected = "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Corynebacteriales.f__Nocardiaceae.g__Nocardia",
                       width = '100%')
        
      }),
      
      mainPanel(

        renderPlot(width = 1000, height = 750, expr = {
          
          if(input$plasmaDataTypeSelected == "Raw Kraken Data"){
            plasmaMetadata <- metadataPSMatchedDP_RAWQC
            plasmaData <- vbMergedDataKrakenDPCF[rownames(plasmaMetadata),]
          } else if(input$plasmaDataTypeSelected == "Raw Shogun Data"){
            plasmaMetadata <- metadataPSMatchedDP_RAWQC
            plasmaData <- dataPSUniqueDP[rownames(plasmaMetadata),]
          } else if(input$plasmaDataTypeSelected == "Decontaminated & Age-Sex Normalized Kraken Data"){
            plasmaData <- snmDataKrakenCFDecontamDPQC
            plasmaMetadata <- metadataPSMatchedDPQCFiltered
          } else if(input$plasmaDataTypeSelected == "Decontaminated & Age-Sex Normalized Shogun Data"){
            plasmaData <- snmShogundataPSUniqueDecontamDPQC
            plasmaMetadata <- metadataPSMatchedDPQCFiltered
          }

          microbeName <- input$plasma_microbe_name
          mergedPlasmaData <- cbind(plasmaMetadata, plasmaData)

          mergedPlasmaData %>%
            select(disease_type_consol, microbeName) %>%
            # filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
            ggboxplot(x = "disease_type_consol", y = microbeName,
                      color = "disease_type_consol",
                      palette = "lancet",
                      # notch = TRUE,
                      # facet.by = "sample_type",
                      # nrow = 3,
                      xlab = "Sample Type", ylab = "Abundance (Raw or Normalized [log2-cpm])",
                      title = paste(microbeName, "\n\nAbundances (cell-free DNA)"),
                      legend = "top",
                      legend.title = "Sample Type") +
            rotate_x_text(angle = 30) +
            theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
            stat_compare_means(label.x = 2)

        })
      )
      
      
      
      )
  })
  
}

shinyApp(ui, server)