#EWAS catalogue summary associations

library(tidyverse)
library(data.table)


#EWAS catalogue search

#1. load ewas catalogue download
#2. load bayesr results
#3. create df of paired CpG and trait from bayesr results 
    #For ewas catalogue:
        #Filter to epigenome-wide significance
        #Filter to whole blood only
        #Have not filtered to N participants as so few protein EWASs
#4. Filter ewas catalogue to CpGs of interest from BayesR
#5. Manually inspect filtered ewas catalogue for traits -> look for different ways of describing proteins
#6. If required:  manually filter ewas catalogue results by trait
#7. And/or create replicated association df using R

#load libraries
library(data.table)
library(R.utils)
library(stringr)

#1 & 2. load bayesr results & ewas catalogue
all_anno <- fread("cis_trans_anno_trial_1604.csv")
catalogue <- fread("ewascatalog-results_050825.txt.gz")
traits <- fread("ewascatalog-studies_050825.txt.gz")

#3 & 4
#Subset to CpGs from bayesR results & to results with epigenome-wide significance
br_cat <- catalogue %>%
  dplyr::filter(CpG %in% all_anno$CpG & P < 3.6e-8)
dim(br_cat) #1983 entries
length(unique(br_cat$CpG)) #280
length(unique(br_cat$StudyID)) #533

#subset traits file to whole blood and studyIDs in the filtered Cpg results
br_traits <- traits %>%
  dplyr::filter(StudyID %in% br_cat$StudyID & Tissue == "Whole blood")
br_traits <- br_traits %>% filter(N > 1000)
dim(br_traits) #217

#merge to common results fulfilling all criteria
both <- merge(br_cat, br_traits, by = "StudyID") 
dim(both) #1082 entries
length(unique(both$StudyID)) #217
length(unique(both$CpG)) #231 CpGs 
fwrite(fwrite(both, file = "ewascat_bayescpgs.csv"))

#summarise traits by CpG
result <- both %>%
  group_by(CpG) %>%
  dplyr::summarise(
    Traits = paste(unique(Trait), collapse = ", "),  # List traits as a comma-separated string
    Num_Traits = n_distinct(Trait))


#df of just associations for merging
associations <- all_anno %>%
  dplyr::select(CpGs, probe_gene, protein, protein_gene)
associations <- associations %>%
  mutate(protein_gene = stringr::str_split_i(protein_gene, " ", i = 1),
        probe_gene = stringr::str_split_i(probe_gene, ";", i = 1)) %>%
  rename(CpG = CpGs)

assoc <- associations %>%
  group_by(CpG) %>%
  summarise(
    probe_gene = probe_gene,
    proteins = paste(unique(protein), collapse = ", "),
    protein_gene = paste(unique(protein_gene), collapse = ", "),
    number_assoc_proteins = length(unique(protein))) %>%
    ungroup() %>%
    distinct()

sumdf <- left_join(assoc, result, by = "CpG")
sumdf <- sumdf %>% arrange(CpG)
sumdf <- sumdf %>% select(CpG, probe_gene, everything())
head(sumdf)
sumdf <- sumdf %>% arrange(desc(Num_Traits))

#write out both sets of results
fwrite(result, file = "catsearch_summary.csv")
fwrite(sumdf, file = "pgs_traits_summary_050825.csv")

#5 & 6 On visual inspection of traits associated with CpGs from ewas catalogue: 
   #proteins are almost exclusively listed by gene name
   #no UniProt IDs used
   #any others manually looked up and UniProt IDs checked against list of proteins from bayesr results
   #manually filtered list created

#7. Filter "both" dataframe by protein gene to find replicated associations
protein_gene_list <- unique(all_anno$protein_gene)
head(protein_gene_list)
search_traits <- protein_gene_list
pattern <- paste(search_traits, collapse = "|")
pattern #manually edit to clean up any special characters etc
pattern <- "IGKV3-7|IGLV4-69|IGLV8-61|IGLV4-60|IGLV2-18|IGHV3OR16-12|IGKV3D-15|IGHV1-45|IGHV3-49|IGKJ1|IGHV6-1|IGHV3-15|IGHV2-26|IGKV3D-20|IGHV1-3|IGHV4-28|IGHV3-35|IGHV3-38|IGHV5-51|IGHA2|IGHV3-64D|IGLC7|FBLN1|TF|GC|APOL1|CD5L|FCN3|ATRN|APOM|CP|F2|F9|F10|PLG|F12|SERPINC1|SERPINA1|AGT|A2M|C3|C5|KNG1|IGF2|IGKV1-5|IGKV3-20|IGLV1-47|IGLV1-51|IGLV1-40|IGLV2-23|IGLV2-11|IGLV2-8|IGLV3-19|IGLV3-27|IGHV1-69|IGHV1-46|IGHV3-11|IGHV3-13|IGHV3-7|IGKC|HBD|APOA1|APOE|FGA|APCS|C1QC|C9|APOH|LRG1|FN1|AMBP|ORM1|AHSG|PPBP|TF|HPX|C4BPA|VTN|APOB|LCAT|HRG|IGLV7-43|A1BG|VWF|IGKV1-16|F13B|SERPINA7|SERPIND1|IGKV2-30|IGKV4-1|GSN|APOA4|C8A|C8G|THBS1|SERPINA6|CD14|CFH|C7|C6|SELL|CPN1|LBP|ORM2|ITIH1|PZP|C4BPB|CPN2|FBLN1|AZGP1|PON1|SERPINA4|IGFALS|SERPINF1|BTD|AFM|LUM|APOC4|HBB|HBA1|ADARB1|GPLD1|LGALS3BP|APOF|HABP2|ADIPOQ|ECM1|PGLYRP2|CDC5L|CFHR5|LYVE1"

#filter CpG df to any rows containing a trait matching the search term
trait_filtered <- both %>% dplyr::filter(
  grepl(pattern, Trait))

#after visual inspection of trait confirming protein genes are listed first
#create new column with just the gene for matching with associations
trait_filtered <- trait_filtered %>%
  mutate(protein_gene = stringr::str_split_i(Trait, " ", i = 1))

#df of just associations for merging
associations <- all_anno %>%
  dplyr::select(CpGs,protein_gene)
associations <- associations %>%
  mutate(protein_gene = stringr::str_split_i(protein_gene, " ", i = 1))

#subset ewas catalogue results to replicated associations
associations <- associations %>% dplyr::rename(CpG = CpGs)
associations_repl <- left_join(trait_filtered, associations, by = c("CpG", "protein_gene"))  
dim(associations_repl)

fwrite(associations_repl, file ="replicated_associations.csv")

#extract details of replication
length(unique(associations_repl$CpG))
length(unique(associations_repl$StudyID)) 
length(unique(associations_repl$protein_gene)) 

associations_repl %>% dplyr::select(P) %>% summary()


#Checking annotations
head(all_anno)
sum(is.na(all_anno$CIS_TRANS))


#Extract dg results from ewas catalogue
dgr <- catalogue %>% filter(grepl("somascan_protein_measurement", StudyID))

#Extract dg traits from traits catalogue
dgt <- traits %>% filter(PMID == "35945220")
dgt <- dgt %>% mutate(protein_gene = str_split_i(Trait, " ", 1))

#add annotations for filtering
dt <- merge(dgr, dgt, by = "StudyID")
anno <- fread("133protein_annos_forpub_nodup.csv")
anno <- anno %>% rename(
            uniprot = protein,
            protein_gene = gene_names_simpl)

dt2 <- left_join(dt, anno, by = "protein_gene")

#Look at results for overlapping proteins
soma3 <- fread("GS+soma+QC+and+study+samples+normalized.csv")
somaids <- fread("somascan_ids.txt")

protein_list <- colnames(protein_df)

#Filter to overlapping proteins
soma_list <- colnames(soma3)
somaf <- soma_list[soma_list %in% protein_list] 
overlap <- unique(somaf)
length(unique(somaf))

#look at results by threshold for overlap proteins
dto <- dt2 %>% filter(uniprot %in% overlap)

#summarise results by p-value threshold
thresholds <- c(4.5e-10, 3.6e-8, 0.05)

summarise_p <- function(dt, thresholds) {
  summary_list <- lapply(thresholds, function(x) {
    subset_dt <- subset(dto, P <= x)
    summary <- data.frame(
      threshold = x,
      total_results = nrow(subset_dt),
      CpG_count = length(unique(subset_dt$CpG)),
      protein_count = length(unique(subset_dt$uniprot)))
    
    return(summary)
  })
  
    # Combine all summaries into one data frame
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)

}

# Call the function with the data frame and thresholds
summary_results <- summarise_p(dto, thresholds)
print(summary_results)

#check results 
dte <- dto %>% filter(P < 3.6e-8)

#___________________________________________________________
#Add JR results into mix - look at results that overlap 
head(all_anno)
jr <- all_anno %>% select(protein, CpGs, pip,probe_gene, protein_gene, effect_size) %>%
  rename(CpG = CpGs)

comb <- left_join(dt, jr, by = c("protein_gene", "CpG"))

#summarise results by p-value threshold
thresholds <- c(4.5e-10, 3.6e-8, 0.05)

summarise_p <- function(dt, thresholds) {
  summary_list <- lapply(thresholds, function(x) {
    subset_dt <- subset(dt, P <= x)
    summary <- data.frame(
      threshold = x,
      total_results = nrow(subset_dt),
      CpG_count = length(unique(subset_dt$CpG)),
      protein_count = length(unique(subset_dt$protein_gene)),
      n_assoc_gmrm =  sum(!is.na(subset_dt$pip)),
      cpg_gmrm = length(unique(na.omit(subset_dt$probe_gene))),
      protein_gmrm = length(unique(na.omit(subset_dt$protein))))   #contains NA values which get counted as 1
     

    
    return(summary)
  })
  
    # Combine all summaries into one data frame
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)

}

# Call the function with the data frame and thresholds
summary_results2 <- summarise_p(comb, thresholds)
print(summary_results2)
summary_results2 <- summary_results2 %>% mutate(prop_assoc_in_gmrm = n_assoc_gmrm/total_results)
fwrite(summary_results2, file = "summary_somascan_comp.csv")

com2 <- comb %>% mutate(p_thresh = ifelse(P < 4.5e-10, "bonf",
                                     ifelse(P < 3.6e-08 & P > 4.5e-10, "ewas",
                                       ifelse(P < 0.05 & P > 3.6e-08, "nominal", "non"))))  

com2 <- com2 %>% mutate(gmrm = ifelse(!is.na(pip), "yes", "no"))
com3 <- com2 %>% filter(gmrm == "yes")
fwrite(com3, file =  "somascan_comp_070825.csv")

