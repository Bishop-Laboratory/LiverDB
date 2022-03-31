#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(shiny)
library(shinycssloaders)
library(prompter)
library(dplyr)
library(tidyr)
#library(ComplexHeatmap)
library(tibble)
library(futile.logger)
library(ggplot2)

#-------------------------------------------------------------------------------
# HELPER FUNCTIONS, CONSTANTS & INITIALIZATIONS
#-------------------------------------------------------------------------------

source("ui_globals.R")
source("utils.R")

GEO_BASE <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
GENECARDS_BASE <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="

app_data <- readRDS("app_data.rds")
exps <- app_data[["exps"]]
degs <- app_data[["degs"]]
deg_contrasts <- app_data[["contrasts"]]

#-------------------------------------------------------------------------------
# UI
#-------------------------------------------------------------------------------

ui <- function(request) {
  tagList(
    tags$head(
      tags$style(HTML(headerHTML())),
      tags$script(src="https://kit.fontawesome.com/5071e31d65.js", crossorigin="anonymous"),
      tags$link(rel="stylesheet", type="text/css", href="https://cdnjs.cloudflare.com/ajax/libs/cookieconsent2/3.1.1/cookieconsent.min.css")
    ),
    tags$body(
      tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/cookieconsent2/3.1.1/cookieconsent.min.js", `data-cfasync`="false"),
      # tags$script(src="cookie_consent.js")  # Uncomment for cookie consent form
    ),
    use_prompt(),
    navbarPage(
      title = "FibroDB",
      id = "fibrodb",
      theme = bslib::bs_theme(bootswatch = "cosmo"),
      tabPanel(
        title = "Home",
        id = "home-tab",
        value = "aboutTab",
        icon = icon("home"),
        fluidPage(br(), includeHTML("www/home.html"))
      ),
      tabPanel(
        title = "Explore",
        id = "explore-tab",
        icon = icon('table'),
        ExplorePageContents(deg_contrasts) #ExplorePageContents(results)
      ),
      tabPanel(
        title = "Download",
        id = "download-tab",
        icon = icon('download'),
        DownloadPageContents()
      ),
      tabPanel(
        title = "Documentation",
        id = "docs-tab",
        icon = icon('file-alt'),
        tags$iframe(
          src = './documentation.html',
          width = '100%', height = '800px',
          frameborder = 0,
          scrolling = 'auto'
        )
      )
    ), 
    tags$footer(HTML(footerHTML()))
  )
}

?DT::formatSignif
#-------------------------------------------------------------------------------
# SERVER FUNCTION
#-------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  ## DEG results table
  output$degTable <- DT::renderDT(server = TRUE, {
    #req(input$selectStudy)
    study <- input$selectStudy
    pair <- input$selectContrast
    pair <- strsplit(pair, " vs. ")[[1]]
    degs[[study]] %>%
      filter(numerator == pair[[1]] & denominator == pair[[2]]) %>% 
      select(gene_name, logFC, FDR) %>% 
      mutate(
        gene_name = paste0(
          "<a href='", paste0(GENECARDS_BASE, gene_name),
          "' target='_blank'>", gene_name, "</a>"
        )
      ) %>%
      DT::datatable(
        selection = list(mode = "single", selected = 1),
        rownames = FALSE, escape = FALSE,
        colnames = c("Gene", "Fold Change (log2)", "Adjusted p-value"),
        options = list(pageLength = 8, scrollX = TRUE)
      ) %>% 
      DT::formatSignif(2:3, digits = 5)
  })
  
  # Reactive contrast selector
  output$UIselectContrast <- renderUI({
    study <- input$selectStudy
    cont_labels <- deg_contrasts %>%
      filter(study_id == study) %>% 
      unite("contrast", c("numerator", "denominator"), sep = " vs. ") %>% 
      pull(contrast)
    
    selectInput(
      inputId = "selectContrast",
      label = "Contrast (Numerator vs. Denominator)",
      selected = cont_labels[1],
      choices = cont_labels
    )
  })
  
  # Reactive link to study
  output$studyLink <- renderUI({
    study <- input$selectStudy
    
    helpText(
      "View study ",
      a(study, href = paste0(GEO_BASE, study), target = "_blank"),
      " on GEO."
    )
  })
}

#-------------------------------------------------------------------------------
# LAUNCH APP
#-------------------------------------------------------------------------------

shinyApp(ui, server)



#-------------------------------------------------------------------------------
# SCRATCH ZONE
#-------------------------------------------------------------------------------




