#-------------------------------------------------------------------------------
# LIBRARIES
#-------------------------------------------------------------------------------
library(shiny)
library(shinycssloaders)
library(prompter)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
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
metadata <- app_data[["metadata"]]
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


#-------------------------------------------------------------------------------
# SERVER FUNCTION
#-------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Reactive link to study
  output$studyLink <- renderUI({
    study <- input$selectStudy
    
    helpText(
      "View study ",
      a(study, href = paste0(GEO_BASE, study), target = "_blank"),
      " on GEO."
    )
  })
  
  # Reactive contrast selector
  output$UIselectContrast <- renderUI({
    study <- input$selectStudy
    cont_labels <- deg_contrasts %>%
      dplyr::filter(study_id == study) %>% 
      unite("contrast", c("numerator", "denominator"), sep = " vs. ") %>% 
      pull(contrast)
    
    selectInput(
      inputId = "selectContrast",
      label = "Contrast (Numerator vs. Denominator)",
      selected = cont_labels[1],
      choices = cont_labels
    )
  })
  
  # DEG results table
  output$degTable <- DT::renderDT(server = TRUE, {
    req(input$selectStudy, input$selectContrast)
    study <- input$selectStudy
    pair <- strsplit(input$selectContrast, " vs. ")[[1]]

    degs[[study]] %>%
      dplyr::filter(numerator == pair[[1]] & denominator == pair[[2]]) %>% 
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
  
  # Get currently-selected gene from DEG results table
  current_gene <- reactive({
    req(input$selectStudy, input$selectContrast)
    study <- input$selectStudy
    pair <- strsplit(input$selectContrast, " vs. ")[[1]]
    
    # Get selected row from datatable
    selectedRow <- ifelse(
      is.null(input$degTable_rows_selected), 1, input$degTable_rows_selected
    )
    
    gene <- degs[[study]] %>%
      dplyr::filter(numerator == pair[[1]] & denominator == pair[[2]]) %>% 
      dplyr::filter(row_number() == selectedRow) %>%
      pull(gene_name)
  })
  
  # Expression plot
  output$countplot <- plotly::renderPlotly({
    req(input$selectStudy, input$selectCTS, current_gene())
    study <- input$selectStudy
    cts_sel <- input$selectCTS
    gene <- current_gene()

    plt <- exps[[study]] %>%
      dplyr::filter(gene_name == gene) %>% 
      rename(expression = tolower(cts_sel)) %>% 
      left_join(
        metadata %>% select(sample_id, condition),
        by = c("sample_id")
      ) %>% 
      ggplot(
        aes(x = condition, y = expression, fill = condition)
      ) +
      geom_boxplot(width = .65, alpha = .6, outlier.shape = NA) +
      geom_jitter(width = .15) +
      xlab("Conditions of Samples") +
      ylab(paste0("Expression (", cts_sel, ")")) +
      theme_gray(base_size = 13) +
      ggtitle(gene) +
      theme(legend.position = "none")
    
    plotly::ggplotly(plt)
  })
  
  # Volcano plot
  output$volcanoPlot <- renderPlot({
    study <- input$selectStudy
    pair <- strsplit(input$selectContrast, " vs. ")[[1]]
    gene <- current_gene()
    
    toplt <- degs[[study]] %>%
      dplyr::filter(numerator == pair[[1]] & denominator == pair[[2]]) %>% 
      rename(padj = FDR, fc = logFC)
      
    req(!is.na(toplt$padj[1]))
    
    plot_title <- paste0(pair[1], " vs. ", pair[2])
    pltdata <- toplt %>%
      mutate(
        hlight = gene_name == gene,
        padj = case_when(
          padj == 0 ~ .Machine$double.xmin, TRUE ~ padj
        ),
        padj = -log10(padj),
        sigcond = case_when(
          padj < 3 ~ "n.s.",
          abs(fc) < 1 ~ "sig-only",
          fc > 1 ~ "Over-expressed",
          fc < -1 ~ "Under-expressed"
        )
      ) %>%
      arrange(hlight, desc(padj)) 
    maxval <- max(pltdata$padj)
    pltdata %>%
      ggplot(
        aes(x = fc, y = padj, color = sigcond, size = hlight)
      ) +
      geom_vline(xintercept = 1, linetype = "dashed", alpha = .25) +
      geom_vline(xintercept = -1, linetype = "dashed", alpha = .25) +
      geom_hline(yintercept = 3, linetype = "dashed", alpha = .25) +
      geom_point() +
      xlab("Log2 Fold Change") +
      ylab("-log10(Adjusted p-value)") +
      guides(size = guide_none(),
             color = guide_legend(title = NULL)) +
      theme_bw(base_size = 18) + 
      scale_y_continuous(expand = c(0,0), limits = c(0, 1.05 * maxval)) +
      ggtitle(plot_title, subtitle = gene) +
      scale_color_manual(
        values = c(
          "n.s." = "#d6d6d6",
          "sig-only" = "#91bac4",
          "Over-expressed" = "#2fa4c2",
          "Under-expressed" = "#c24e2f"
        )
      ) + 
      theme(legend.position = "bottom", legend.text = element_text(size = 16)) 
  })
  
  # Heatmap
  output$heatmap <- renderPlot({
    req(input$selectStudy, input$selectContrast, input$selectCTS2)
    study <- input$selectStudy
    pair <- strsplit(input$selectContrast, " vs. ")[[1]]
    cts_sel <- input$selectCTS2
    
    toplt <- degs[[study]] %>% 
      rename(padj = FDR, fc = logFC)
    
    req(!is.na(toplt$padj[1]))
    
    plot_title <- paste0(pair[1], " vs. ", pair[2])
    g2plt <- toplt %>%
      mutate(sigcond = case_when(
        padj > 0.05 ~ "n.s.",
        abs(fc) < 1 ~ "sig-only",
        fc > 1 ~ "Over-expressed",
        fc < -1 ~ "Under-expressed"
      )) %>%
      dplyr::filter(sigcond %in% c("Over-expressed", "Under-expressed")) %>% 
      group_by(sigcond) %>% 
      slice_min(order_by = padj, n = 12) %>%
      pull(gene_name)
    
    topvt <- exps[[study]] %>%
      dplyr::filter(gene_name %in% g2plt) %>% 
      rename(counts = tolower(cts_sel)) %>% 
      left_join(
        metadata %>% select(sample_id, condition),
        by = c("sample_id")
      ) %>%
      dplyr::filter(grepl(pair[1], condition) | grepl(pair[2], condition))
    
    annot <- topvt %>% 
      dplyr::select(sample_id, condition) %>% 
      unique() %>% 
      column_to_rownames("sample_id")
    
    plt <- pivot_wider(  
      data = topvt,
      id_cols = gene_name, names_from = sample_id, values_from = counts
    ) %>% 
      column_to_rownames("gene_name") %>% 
      as.matrix() %>% 
      pheatmap(
        scale = "row",
        angle_col = "45",
        annotation_col = annot,
        name = cts_sel,
        main = plot_title
      )
    
    plt
  })

  
  
  
}


#-------------------------------------------------------------------------------
# SCRATCH ZONE
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# LAUNCH APP
#-------------------------------------------------------------------------------

shinyApp(ui, server)

