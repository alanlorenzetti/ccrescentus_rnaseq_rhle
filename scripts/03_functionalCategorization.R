# alorenzetti 202106

# description ####
# this script will get and wrangle annotation
# and functional data for genes

# setting up annotation ####
# loading annotation file
gffFile = file.path("data/Ccrescentus.gff")
annot = rtracklayer::import(gffFile)

# filtering genes
genes = subset(annot, type == "gene")
genes$Dbxref = genes$Dbxref[lapply(genes$Dbxref, grepl, pattern = "GeneID:")] %>%
  sub(pattern = "GeneID:", replacement = "") %>%
  as.numeric()

# filtering CDS to obtain protein products
CDS = subset(annot, type == "CDS")
CDS$Dbxref = CDS$Dbxref[lapply(CDS$Dbxref, grepl, pattern = "GeneID:")] %>%
  sub(pattern = "GeneID:", replacement = "") %>%
  as.character()

CDS = CDS[CDS$Dbxref %>% duplicated == FALSE,]
CDS = as_tibble(CDS)
CDS = CDS %>%
  dplyr::select(Dbxref, product, ID)

# filtering ncRNAs to obtain products
types = c("ncRNA",
          "tRNA",
          "SRP_RNA",
          "antisense_RNA",
          "tmRNA",
          "RNase_P_RNA",
          "rRNA")

ncrna = subset(annot, type %in% types)
ncrna$Dbxref = ncrna$Dbxref[lapply(ncrna$Dbxref, grepl, pattern = "GeneID:")] %>%
  sub(pattern = "GeneID:", replacement = "") %>%
  as.character()

ncrna = ncrna[ncrna$Dbxref %>% duplicated == FALSE,]
ncrna = as_tibble(ncrna)
ncrna = ncrna %>%
  dplyr::select(Dbxref, product, ID)

names(genes) = paste(genes$ID, genes$locus_tag, genes$Dbxref, genes$Name, sep = "|")

# getting cog info ####
if(file.exists("data/cog.tsv")){
  cogdict = read_tsv("data/cog.tsv")
} else {
  # this chunk will download
  # COG 2020 data to perform
  # functional caracterization
  
  # downloading files directly from COG ftp ####
  cog = list()
  
  # description of each file is available at
  # ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-11-25.txt
  
  # cog obj
  # it is going to throw a warning, since number of
  # cols is different depending on the entry
  cog$cog = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv",
                     col_names = c("gene_id",
                                   "refseq_id",
                                   "protein_id",
                                   "protein_length",
                                   "cog_coord_prot",
                                   "cog_length_prot",
                                   "cog_id",
                                   "reserved",
                                   "cog_memb_class",
                                   "PSI_BLAST_bit_score",
                                   "PSI_BLAST_evalue",
                                   "cog_prof_length",
                                   "prot_coord_cog"))
  
  # cog memb class definition:
  # 0: footprint covers most of the protein and most of the COG profile;
  # 1: footprint covers most of the COG profile and part of the protein;
  # 2: footprint covers most of the protein and part of the COG profile;
  # 3: partial match on both protein and COG profile)
  
  # cog definitions
  cog$def = read_tsv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab",
                     col_names = c("cog_id",
                                   "cog_funcat",
                                   "cog_name",
                                   "gene_symbol",
                                   "functional_pathway",
                                   "pubmed_id",
                                   "pdb_id")) #%>% 
  # mutate(cog_funcat = str_replace_all(string = cog_funcat,
  #                                     pattern = "([A-Z])",
  #                                     replacement = "\\1,"),
  #        cog_funcat = str_replace(string = cog_funcat,
  #                                 pattern = ",$",
  #                                 replacement = "")) %>% 
  # separate_rows(cog_funcat, sep = ",")
  
  # some cogs have two categories, but that makes it
  # complex since a gene can have more than one COG
  # a simple solution is to keep only the first
  # category, and that is going to be applied here
  cog$def = cog$def %>% 
    mutate(cog_funcat = str_replace(cog_funcat,
                                    "([A-Z]).*$",
                                    "\\1"))
  
  # cog patt
  # it is going to throw a warning, since number of
  # cols is different depending on the entry
  # cog$patt = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.patt.txt",
  #                    col_names = F)
  # 
  # # cog tax
  # cog$tax = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.tax.csv",
  #                     col_names = F)
  
  # cog function
  cog$fun = read_tsv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab",
                     col_names = c("cog_funcat",
                                   "cog_color_hex",
                                   "cog_category"))
  
  # list of organisms
  cog$org = read_csv(file = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.org.csv",
                     col_names = c("refseq_id", "name", "tx_id", "taxon")) %>% 
    filter(name == "Caulobacter_vibrioides_CB15")
  
  # parsing ####
  # getting Halobacterium salinarum NRC-1
  # refseq assembly id
  rsid = cog$org$refseq_id
  
  # filtering cog object to include
  # only halo
  cog$cog = cog$cog %>%
    filter(refseq_id == rsid)
  
  # unifying datasets
  # keeping only member class 0 or 1
  # according to what have been described above
  cog$final = cog$cog %>%
    left_join(x = ., y = cog$def,
              by = c("cog_id" = "cog_id")) %>% 
    left_join(x = ., y = cog$fun,
              by = c("cog_funcat" = "cog_funcat")) %>% 
    filter(cog_memb_class == 0 | cog_memb_class == 1) %>% 
    dplyr::select(gene_id,
                  cog_id,
                  cog_funcat,
                  cog_name,
                  gene_symbol,
                  functional_pathway,
                  cog_category,
                  cog_memb_class) %>% 
    mutate(functional_pathway = case_when(is.na(functional_pathway) ~ "Undefined",
                                          TRUE ~ as.character(functional_pathway)),
           gene_symbol = case_when(is.na(gene_symbol) ~ "Undefined",
                                   TRUE ~ as.character(gene_symbol))) %>% 
    group_by(gene_id) %>% 
    summarise(across(.cols = everything(),
                     .fns = ~ paste0(.x, collapse = "|"))) %>% 
    ungroup()
  
  # removing redundancy in case
  # a gene has two different COGs
  # of the same category
  removRed = function(x){
    noRed = str_split(string = x,
                      pattern = "\\|") %>%
      unlist() %>% 
      unique() %>% 
      paste0(collapse = "|")
    
    return(noRed)
  }
  
  # performing operations
  cog$final = cog$final %>% 
    rowwise() %>% 
    mutate(across(.cols = everything(),
                  .fns = ~ removRed(.x)))
  
  # cog only exists for ccrescentus cb15
  # we are working with na1000
  # ive got corresponding cb15 gis to na1000 locus_tag
  # from ortholuge db
  # to infer cog of na1000 based on cb15
  # reading ortholuge
  gi = read_delim("http://www.pathogenomics.sfu.ca/ortholugedb/paired/download/plain?strain1=Caulobacter%20crescentus%20CB15&strain2=Caulobacter%20crescentus%20NA1000", delim="\t")
  
  dict = gi %>% 
    dplyr::select(ltcb15 = "Locus Tag (Strain 1)",
                  ltna1000 = "Locus Tag (Strain 2)") %>% 
    dplyr::mutate(ltcb15 = as.character(ltcb15),
                  ltna1000 = as.character(ltna1000)) %>% 
    dplyr::group_by(ltna1000) %>% 
    filter(row_number()==1)
  
  # wrangling dict obj and merging with cog obj
  cogdict = dict %>% 
    left_join(x = .,
              y = cog$final,
              by = c("ltcb15" = "gene_id")) %>% 
    dplyr::select(ltna1000,
                  gene_symbol,
                  cog_id,
                  cog_name,
                  cog_category,
                  functional_pathway) %>% 
    drop_na() %>% 
    distinct()
  
  # writing cog file
  write_tsv(x = cogdict,
            file = "data/cog.tsv")
}

# creating a functional category table for all genes
funCat = genes %>%
  as_tibble() %>%
  dplyr::select(locus_tag, Dbxref) %>%
  dplyr::mutate(Dbxref = as.character(Dbxref)) %>% 
  dplyr::left_join(CDS, by=c("Dbxref" = "Dbxref")) %>% 
  dplyr::left_join(ncrna, by=c("Dbxref" = "Dbxref")) %>% 
  tidyr::unite(col = "product", c("product.x", "product.y"), na.rm = T) %>% 
  tidyr::unite(col = "ID", c("ID.x", "ID.y"), na.rm = T) %>% 
  dplyr::left_join(cogdict, by=c("locus_tag" = "ltna1000")) %>% 
  dplyr::rename(entrezID = "Dbxref") %>% 
  dplyr::mutate(across(.cols = everything(),
                       .fns = ~ case_when(is.na(.x) ~ "Undefined",
                                          TRUE ~ as.character(.x)))) %>% 
  dplyr::distinct()

# saving tables to store funCat object
write_tsv(x = funCat,
          file = paste0("results/", "functionalCategorization", ".tsv"))

# declaring functional categorization function ####
functCat = function(sigTable){
  funcat = sigTable[,c(1:4,6,10)] %>%
    dplyr::mutate(upORdown = case_when(log2FoldChange <= -log2fcthreshold ~ "DownRegulated",
                                       log2FoldChange >= log2fcthreshold ~ "Upregulated")) %>%
    dplyr::left_join(x = .,
                     y = funCat %>% 
                       dplyr::select(c(locus_tag,
                                       product,
                                       cog_id,
                                       cog_name,
                                       cog_category,
                                       functional_pathway)),
                     by="locus_tag") %>%
    dplyr::arrange(log2FoldChange) %>% 
    dplyr::mutate(across(.cols = everything(),
                         .fns = ~ case_when(is.na(.x) ~ "Undefined",
                                            TRUE ~ as.character(.x))))
  
  #  funcat$UNIPROTKB = paste0("<a href='https://www.uniprot.org/uniprot/",funcat$UNIPROTKB,"'>", funcat$UNIPROTKB,"</a>")
  return(funcat)
}
