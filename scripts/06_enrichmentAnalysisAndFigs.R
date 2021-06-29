# alorenzetti 202106

# description ####
# this script will perform
# enrichment analysis inside
# upregulated and downregulated groups
# and will also plot figures

# starting processing ####
# creating and parsing objx ####
# we are going to search for enriched
# upregulated and downregulated set of genes
# I will have to arbitrarily keep only one
# category for each gene
# I am going to keep only the first one
enrichdf = list()

for(i in contrasts){
  enrichdf[[i]][["downregulated"]] = left_join(x = finRes[[i]][["downregulated"]],
                                                                      y = funCat,
                                                                      by = "locus_tag") %>%
    dplyr::select(locus_tag,
                  geneName,
                  log2FoldChange,
                  padj,
                  cog_category) %>%
    dplyr::mutate(status = case_when(log2FoldChange >= 1 ~ "Upregulated",
                                     log2FoldChange <= -1 ~ "Downregulated"),
                  cog_category = str_replace(cog_category, "\\|.*$", ""),
                  cog_category = str_replace(cog_category, ";.*$", ""))
  
  enrichdf[[i]][["upregulated"]] = left_join(x = finRes[[i]][["upregulated"]],
                                                                    y = funCat,
                                                                    by = "locus_tag") %>%
    dplyr::select(locus_tag,
                  geneName,
                  log2FoldChange,
                  padj,
                  cog_category) %>%
    dplyr::mutate(status = case_when(log2FoldChange >= 1 ~ "Upregulated",
                                     log2FoldChange <= -1 ~ "Downregulated"),
                  cog_category = str_replace(cog_category, "\\|.*$", ""),
                  cog_category = str_replace(cog_category, ";.*$", ""))
  
  enrichdf[[i]][["all"]] = bind_rows(enrichdf[[i]][["downregulated"]],
                                                            enrichdf[[i]][["upregulated"]])
  
  enrichdf[[i]][["allfig"]] = enrichdf[[i]][["all"]] %>%
    group_by(cog_category, status) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(count = case_when(status == "Downregulated" ~ count * -1,
                             TRUE ~ as.numeric(count)))
  
  enrichdf[[i]][["funcat"]] = funCat %>%
    mutate(cog_category = str_replace(cog_category, "\\|.*$", ""))
  
  # # enrichment analysis ####
  # # regulation status as primary clusters
  # # hypergeometric enrichment test of
  # # cog_category var inside
  # # regulation status
  vars = "cog_category"
  
  # # creating list to store results
  enrich = list()
  inputobj = enrichdf[[i]][["all"]]
  funcatobj = enrichdf[[i]][["funcat"]]
  for(j in inputobj$status %>% unique()){
    curRegGroup = inputobj %>%
      filter(status == j)
    for(k in vars){
      curRegGroupVec = curRegGroup %>%
        dplyr::select(k) %>%
        unlist(use.names = F)
      
      curRegGroupLvs = curRegGroupVec %>%
        unique()
      
      for(l in curRegGroupLvs){
        wb = sum(curRegGroupVec == l)
        vecu = funcatobj %>%
          dplyr::select(k) %>%
          unlist(use.names = F)
        wu = sum(vecu == l)
        bu = sum(vecu != l)
        drawn = curRegGroupVec %>% length()
        
        pval = phyper(q= wb, m= wu, n= bu, k = drawn, lower.tail = F)
        
        tib = tibble(regRule = j,
                     criteria = k,
                     level = l,
                     pval = pval)
        enrich = bind_rows(enrich, tib)
      }
    }
  }
  # 
  # # correcting pvalues using BH method
  # # filtering by pval
  enrich$qval = p.adjust(enrich$pval, method = "BH")
  enrich = enrich[enrich$qval < qthr,]
  
  # adjusting dataset to include enrichment analysis
  # on differential expression category plot
  enrichdf[[i]][["allfig"]] = enrichdf[[i]][["allfig"]] %>%
    left_join(x = .,
              y = enrich,
              by = c("cog_category" = "level", "status" = "regRule")) %>%
    dplyr::mutate(enrichStatus = case_when(qval < 0.001 ~ "***",
                                           qval < 0.01 & qval >= 0.001 ~ "**",
                                           qval < 0.05 & qval >= 0.01 ~ "*",
                                           TRUE ~ NA_character_),
                  cog_category = factor(cog_category, cog_category %>% unique() %>% rev()))
  
  # plotting figure ####
  enrichdf[[i]][["plot"]] = enrichdf[[i]][["allfig"]] %>%
    ggplot(aes(x = count,
               fill = status,
               y = cog_category)) +
    geom_col() +
    geom_text(aes(label = enrichStatus,
                  group = status),
              vjust = 0.75) +
    scale_fill_manual(name = "",
                      values = c("Upregulated" = "#E15759",
                                 "Downregulated" = "#4E79A7")) +
    scale_x_continuous(name = "Count",
                       limits = c(-100, 100),
                       labels = abs) +
    theme(text = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          legend.position = "bottom") +
    ylab("COG")
  
  # saving
  ggsave(plot = enrichdf[[i]][["plot"]],
         filename = paste0("plots/", i, "_enrichmentOfCategories.png"),
         width = 7,
         height = 4,
         unit = "in",
         dpi = 600)
  
  ggsave(plot = enrichdf[[i]][["plot"]],
         filename = paste0("plots/", i, "_enrichmentOfCategories.pdf"),
         width = 7,
         height = 4,
         unit = "in",
         dpi = 600)
}
