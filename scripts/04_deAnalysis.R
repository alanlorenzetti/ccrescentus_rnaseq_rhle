# alorenzetti 202106

# description ####
# this script will prepare data
# generate the counts, and
# perform differential expression
# analysis

# preparing data ####
# providing info about data
samples = c("NA1_10C", "NA2_10C", "NA3_10C",
            "NA1_30C", "NA2_30C", "NA3_30C",
            "rhIE_1_10oC_S4", "rhIE_2_10oC_S5", "rhIE_3_10oC_S6",
            "rhIE_1_30oC_S1", "rhIE_2_30oC_S2", "rhIE_3_30oC_S3")

reps = rep(c("A", "B", "C"), 4)
strains = c(rep("NA1000", 6), rep("rhlE", 6))
conditions = c(rep(c(rep("10C", 3), rep("30C", 3)), 2))

# info datatable
colData = data.frame(row.names = samples,
                     replicate = reps,
                     strain = strains,
                     condition = conditions)

# loading raw counts ####
# loading raw counts obtained by kallisto
# to manually compute TPMs, follow the instructions on the following
# page: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# reading totalrna counts
totrna = read_tsv("data/tableEstCounts25.tsv")

# formatting deseq2 input object
assay = totrna %>% 
  dplyr::select(starts_with("NA"),
                starts_with("rhIE")) %>%
  as.matrix() %>% 
  round(digits = 0)

colnames(assay) = colnames(assay) %>%
  str_replace(., "^NA", "NA1000") %>%
  str_replace(., "rhIE_", "rhlE") %>%
  str_replace(., "oC", "C") %>%
  str_replace(., "_S[0-9]", "") %>%
  str_replace(., "([0-9])_", "_\\1_") %>%
  str_replace(., "(.*)_(.*)_(.*)", "\\1_\\3_\\2") %>%
  str_replace(., "1$", "A") %>%
  str_replace(., "2$", "B") %>%
  str_replace(., "3$", "C")

assay = assay[,c(1,3,5,2,4,6,7,9,11,8,10,12)]

rownames(colData) = colnames(assay)

totrnaSE = SummarizedExperiment(assay = list(counts=assay),
                                rowData = totrna[,c(1:3)],
                                colData = colData)

# giving gene names thinking about compatibility
# issues. I am gonna add some features available
# in the functional categorization object in order
# to make it work in the next script
rownames(totrnaSE) = rowData(totrnaSE) %>% 
  as_tibble() %>% 
  dplyr::select(target_id) %>% 
  mutate(target_id = str_replace(string = target_id,
                                 pattern = "\\|.*$",
                                 replacement =  "")) %>% 
  left_join(x = .,
            y = funCat,
            by = c("target_id" = "locus_tag")) %>% 
  unite(rn, c("ID",
              "target_id",
              "entrezID",
              "gene_symbol"),
        sep = "|") %>% 
  dplyr::select(rn) %>% 
  unlist(use.names = F)

# removing rRNA instances
totrnaSE = totrnaSE[str_detect(string = rownames(totrnaSE),
                               pattern = "CCNA_R0069|CCNA_R0066",
                               negate = T),]

# creating DESeq object for entire experiment ####
# to do exploratory analysis and differential expression analysis
# creating deseq2 objects
totrnadds = totrnaSE

totrnadds$group = factor(paste0(totrnadds$strain, "_", totrnadds$condition))
totrnadds = DESeqDataSet(totrnadds, design = ~ group)

# doing rlog transformation for distance and PCA
# before eliminating genes with zero counts
rld = rlog(totrnadds)
rldNonBlind = rlog(totrnadds, blind=F)

# removing genes with zero counts and performing DESeq2 analysis
totrnadds = totrnadds[rowSums(counts(totrnadds)) > 1, ]
totrnadds = DESeq(totrnadds)
#resultsNames(dds)

# setting distance matrix for dds
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste( rld$strain, rld$condition, rld$replicate, sep="_" )
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(rev(brewer.pal(9, "Reds")) )(255)

# result tables for contrasts
results = list()
results[["rhlE10C_vs_rhlE30C"]] = results(totrnadds, contrast= c("group", "rhlE_10C", "rhlE_30C"), alpha = padjthreshold)
results[["rhlE10C_vs_NA100010C"]] = results(totrnadds, contrast= c("group", "rhlE_10C", "NA1000_10C"), alpha = padjthreshold)
results[["rhlE30C_vs_NA100030C"]] = results(totrnadds, contrast= c("group", "rhlE_30C", "NA1000_30C"), alpha = padjthreshold)
