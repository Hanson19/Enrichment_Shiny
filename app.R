#cluster enrichment app
library(shiny)
library(tidyverse)
library(bslib)
library(rms)

enriched_terms <- read.csv("data/Term_Cluster_Count_FDR.csv", header = TRUE)
enriched_terms <- enriched_terms[-c(1)]
enriched_terms$FDR <- as.character(enriched_terms$FDR)

z_scores <- read.csv("data/z_scores_rep_curves.csv", header = TRUE)
z_scores <- z_scores[-c(1)]
cluster_order <- unique(z_scores[c(1,16)])
z_scores$Name <- factor(z_scores$Name, levels = cluster_order[[2]])

ui <- page_sidebar(
   title = "Cluster Enrichment",
   sidebar = sidebar(
     helpText(
       "Select a specific cluster and type of enrichment term to get a table of 
       enriched terms for the cluster. Table includes the FDR for the term, the
       number of clusters that are significant for the term, and lists those clusters.
       Terms are aranged by FDR values, smallest to largest. If
       the table remains blank it means that cluster had no significant terms under
       that gene set."
     ),
     checkboxGroupInput(
       "C",
       label = "Choose Clusters",
       choices = list(
         "All" = "All",
         "All Complex" =  "Complex",
         "All LinearUp" = "LinearUp",
         "All LinearDown" = "LinearDown",
         "Complex-1" = "Complex-1", "Complex-2" = "Complex-2", "Complex-3" = "Complex-3", "Complex-4" = "Complex-4",
         "Complex-5" = "Complex-5", "Complex-6" = "Complex-6", "Complex-7" = "Complex-7", "Complex-8" = "Complex-8",
         "Complex-9" = "Complex-9", "Complex-10" = "Complex-10", "Complex-11" = "Complex-11", "Complex-12" = "Complex-12",
         "Complex-13" = "Complex-13", "Complex-14" = "Complex-14", "LinearUp-1" = "LinearUp-1", "LinearUp-2" = "LinearUp-2",
         "LinearUp-3" = "LinearUp-3", "LinearUp-4" = "LinearUp-4", "LinearUp-5" = "LinearUp-5", "LinearUp-6" = "LinearUp-6",
         "LinearUp-7" = "LinearUp-7", "LinearUp-8" = "LinearUp-8", "LinearDown-1" = "LinearDown-1", "LinearDown-2" = "LinearDown-2",
         "LinearDown-3" = "LinearDown-3", "LinearDown-4" = "LinearDown-4", "LinearDown-5" = "LinearDown-5", "LinearDown-6" = "LinearDown-6"
       ),
       selected = "All"
     ),
     checkboxGroupInput(
       "set",
       "Select All Enrichment Term Types",
       choices = list(
             "SLIM2 GO BP" = "SLIM2 GO BP", "SLIM2 GO CC" = "SLIM2 GO CC", "SLIM2 GO MF" = "SLIM2 GO MF",
             "DRSC GLAD Gene Group" = "DRSC GLAD Gene Group", "FlyBase Gene Group" = "FlyBase Gene Group",
             "REACTOME pathway" = "REACTOME pathway"
       ),
       selected = "SLIM2 GO BP"
     ),
     radioButtons(
       "analysis",
       "Cluster Comparison",
       choices = list("Show All Terms" = "all", "Show Shared Terms between Clusters" = "share", "Show Unique Terms between Clusters" = "unique"),
       selected = "all"
     )
   ),
    page_fillable(
      layout_columns(
        card(plotOutput("Dendrogram")),
        card(plotOutput("Rep_Curves")),
        card(tableOutput("Enriched_Terms")),
        col_widths = c(12,12,12),
        row_heights = c(0.5,0.5,0.50)
      )
    )
   # card(plotOutput("Dendrogram")),
   # card(plotOutput("Rep_Curves")),
   # card(tableOutput("Enriched_Terms"))
)

server <- function(input, output){
  
  output$Dendrogram <- renderPlot({
    dendo <- readRDS("data/dendrogram_recolored.RDS")
    if(input$C == "All"){
      dendo
    }else if(input$C == "Complex"){
      df <- dendo[["layers"]][[1]][["data"]]
      df_selected <- df %>% filter(Cluster %in% c("Complex-1", "Complex-2", "Complex-3", "Complex-4",
                                                  "Complex-5", "Complex-6", "Complex-7", "Complex-8",
                                                  "Complex-9", "Complex-10", "Complex-11", "Complex-12","Complex-13","Complex-14") | col == "black")
      
      df_nonselected <- df %>% filter(!Cluster %in% c("Complex-1", "Complex-2", "Complex-3", "Complex-4",
                                                      "Complex-5", "Complex-6", "Complex-7", "Complex-8",
                                                      "Complex-9", "Complex-10", "Complex-11", "Complex-12","Complex-13","Complex-14") & col != "black")
      
      df_nonselected$col <- paste(df_nonselected$col, "1A", sep = "")
      df_alpha <- rbind(df_selected, df_nonselected)
      dendo_alpha <- dendo
      dendo_alpha[["layers"]][[1]][["data"]] <- df_alpha
      dendo_alpha
    }else if(input$C == "LinearDown"){
      df <- dendo[["layers"]][[1]][["data"]]
      df_selected <- df %>% filter(Cluster %in% c("LinearDown-1", "LinearDown-2",
                                                  "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6") | col == "black")
      
      df_nonselected <- df %>% filter(!Cluster %in% c("LinearDown-1", "LinearDown-2",
                                                      "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6") & col != "black")
      
      df_nonselected$col <- paste(df_nonselected$col, "1A", sep = "")
      df_alpha <- rbind(df_selected, df_nonselected)
      dendo_alpha <- dendo
      dendo_alpha[["layers"]][[1]][["data"]] <- df_alpha
      dendo_alpha
    }else if (input$C == "LinearUp"){
      df <- dendo[["layers"]][[1]][["data"]]
      df_selected <- df %>% filter(Cluster %in% c("LinearUp-1", "LinearUp-2",
                                                  "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                  "LinearUp-7", "LinearUp-8") | col == "black")
      
      df_nonselected <- df %>% filter(!Cluster %in% c("LinearUp-1", "LinearUp-2",
                                                      "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                      "LinearUp-7", "LinearUp-8") & col != "black")
      
      df_nonselected$col <- paste(df_nonselected$col, "1A", sep = "")
      df_alpha <- rbind(df_selected, df_nonselected)
      dendo_alpha <- dendo
      dendo_alpha[["layers"]][[1]][["data"]] <- df_alpha
      dendo_alpha
    }else if(input$C != "All"){
      df <- dendo[["layers"]][[1]][["data"]]
      df_selected <- df %>% filter(Cluster %in% input$C | col == "black")
      df_nonselected <- df %>% filter(!Cluster %in% input$C & col != "black")
      df_nonselected$col <- paste(df_nonselected$col, "1A", sep = "")
      df_alpha <- rbind(df_selected, df_nonselected)
      dendo_alpha <- dendo
      dendo_alpha[["layers"]][[1]][["data"]] <- df_alpha
      dendo_alpha
    }else{
      dendo
    }
  })
  
  output$Rep_Curves <- renderPlot({
    if(input$C == "All"){
      z_scores %>% ggplot(aes(x=day, y=Z_score))+
        ylab("Z Score")+
        xlab("Day")+
        facet_wrap(~Name, nrow=1)+
        geom_smooth(method = "lm",aes(color=Name),formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))),show.legend = FALSE)+
        scale_color_manual(values = unique(z_scores[c(1,20)]) %>% pull(colors))+
        theme_bw()+
        theme(text = element_text(size=10))
    }else if(input$C == "Complex"){
      z_ploted <- z_scores %>% filter(Official_ID %in% c("Complex-1", "Complex-2", "Complex-3", "Complex-4",
                                                         "Complex-5", "Complex-6", "Complex-7", "Complex-8",
                                                         "Complex-9", "Complex-10", "Complex-11", "Complex-12","Complex-13","Complex-14"))
      z_ploted %>% ggplot(aes(x=day, y=Z_score))+
        ylab("Z Score")+
        xlab("Day")+
        facet_wrap(~Name, nrow=1)+
        geom_smooth(method = "lm",aes(color=Name),formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))),show.legend = FALSE)+
        scale_color_manual(values = unique(z_ploted[c(1,20)]) %>% pull(colors))+
        theme_bw()+
        theme(text = element_text(size=15))
    }else if(input$C == "LinearDown"){
      z_ploted <- z_scores %>% filter(Official_ID %in% c("LinearDown-1", "LinearDown-2",
                                             "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6")) 
      z_ploted %>% ggplot(aes(x=day, y=Z_score))+
        ylab("Z Score")+
        xlab("Day")+
        facet_wrap(~Name, nrow=1)+
        geom_smooth(method = "lm",aes(color=Name),formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))),show.legend = FALSE)+
        scale_color_manual(values = unique(z_ploted[c(1,20)]) %>% pull(colors))+
        theme_bw()+
        theme(text = element_text(size=15))
    }else if(input$C == "LinearUp"){
      z_ploted <- z_scores %>% filter(Official_ID %in% c("LinearUp-1", "LinearUp-2",
                                             "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                             "LinearUp-7", "LinearUp-8")) 
      z_ploted%>% ggplot(aes(x=day, y=Z_score))+
        ylab("Z Score")+
        xlab("Day")+
        facet_wrap(~Name, nrow=1)+
        geom_smooth(method = "lm",aes(color=Name),formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))),show.legend = FALSE)+
        scale_color_manual(values = unique(z_ploted[c(1,20)]) %>% pull(colors))+
        theme_bw()+
        theme(text = element_text(size=15))
    }else{
     z_ploted <-  z_scores %>% filter(Official_ID %in% input$C) 
     z_ploted %>% ggplot(aes(x=day, y=Z_score))+
        ylab("Z Score")+
        xlab("Day")+
        facet_wrap(~Name, nrow=1)+
        geom_smooth(method = "lm",aes(color=Name),formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))),show.legend = FALSE)+
        scale_color_manual(values = unique(z_ploted[c(1,20)]) %>% pull(colors))+
        theme_bw()+
       theme(text = element_text(size=15))
    }
  })
  
  output$Enriched_Terms <- renderTable({
    if(input$C == "All"){
      data_table <- enriched_terms
    }else if(input$C == "Complex"){
      data_table <- enriched_terms %>% filter(Cluster %in% c("Complex-1", "Complex-2", "Complex-3", "Complex-4",
                                               "Complex-5", "Complex-6", "Complex-7", "Complex-8",
                                               "Complex-9", "Complex-10", "Complex-11", "Complex-12","Complex-13","Complex-14"))
    }else if(input$C == "LinearDown"){
      data_table <- enriched_terms %>% filter(Cluster %in% c("LinearDown-1", "LinearDown-2",
                                               "LinearDown-3", "LinearDown-4", "LinearDown-5", "LinearDown-6"))
    }else if(input$C == "LinearUp"){
      data_table <- enriched_terms %>% filter(Cluster %in% c("LinearUp-1", "LinearUp-2",
                                                             "LinearUp-3", "LinearUp-4", "LinearUp-5", "LinearUp-6",
                                                             "LinearUp-7", "LinearUp-8"))
    }else{
      data_table <- enriched_terms %>% filter(Cluster %in% input$C)
    }

    if(input$analysis == "all"){
      data_table %>% filter(Gene.Set %in% input$set)
    }else if(input$analysis == "share"){
        terms_count <- data_table %>% count(Gene.Set.ID) %>% filter(n == length(unique(data_table$Cluster)))
        #terms_count <- data_table %>% count(Gene.Set.ID) %>% filter(n == length(input$C))
        data_table %>% filter(Gene.Set.ID %in% terms_count$Gene.Set.ID & Gene.Set %in% input$set) %>% arrange(Cluster)%>% arrange(Gene.Set.ID)
    }else if(input$analysis == "unique"){
      terms_count <- data_table %>% count(Gene.Set.ID) %>% filter(n == 1)
      data_table %>% filter(Gene.Set.ID %in% terms_count$Gene.Set.ID & Gene.Set %in% input$set)
    }
    
  })
  
}

shinyApp(ui = ui, server =server)



