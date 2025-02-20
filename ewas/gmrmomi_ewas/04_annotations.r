#GMRMomi results: CpG and protein annotation
library(tidyverse)
library(data.table)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(biomaRt)

#Annotate CpGs_________________________________________________________________
highPIP <- fread("/data_dir/ewas/bayesr/highpipcpgs_0310.csv")
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- data.frame(epic_anno)
setDT(anno, keep.rownames = "CpG")
anno <- anno %>% filter(., CpG %in% highPIP$CpGs)
anno <- anno %>% dplyr::select(Name, chr, pos, strand, UCSC_RefGene_Name, Relation_to_Island)
head(anno)


anno$chr <- gsub("chr", "", anno$chr) # remove chr from chromosome names
anno$chr <- gsub("X", "23", anno$chr) # replace chr X and Y with 23 and 24
anno$chr <- gsub("Y", "24", anno$chr)
anno$chr <- as.numeric(anno$chr) # set chr as numeric 

anno$UCSC_RefGene_Name	 <- as.factor(anno$UCSC_RefGene_Name) # set genes and strands as factors 
anno$strand <- as.factor(anno$strand)
anno[which(anno$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA #set missing gene names to NA 

#Add annotations into highPIP dataframe
anno$CpGs <- anno$Name
anno <- anno %>% dplyr::select(-Name)
cpg_anno <- left_join(highPIP, anno, by = "CpGs")
chr_length <- read.csv("/ext_dir/chromosomes/CHR_Length.csv")
chr_length$CHR <- as.character(chr_length$CHR)
chr_length <- chr_length %>% dplyr::rename(chr = "CHR")
cpg_anno$chr <- as.character(cpg_anno$chr)
cpg_anno <- left_join(cpg_anno, chr_length, by = "chr")
cpg_anno <- cpg_anno %>% dplyr::select(-c(n, PostIP))
cpg_anno <- cpg_anno %>% dplyr::rename(probe_chr = chr, probe_gene = UCSC_RefGene_Name)

#save out results
fwrite(cpg_anno, "/data_dir/ewas/bayesr/output/new_step2/annotated_highPIP_cpgs2.csv")

#check annotations for completeness
anno <- fread("/data_dir/ewas/bayesr/output/new_step2/annotated_highPIP_cpgs2.csv")
head(anno)
sum(is.na(anno))

#Annotate proteins_________________________________________________________________
protein_anno <- fread("/data_dir/ewas/comparison_data/mass_spec_proteins_annotations.csv")
protein_anno <- protein_anno %>% dplyr::select(From, Gene.Names)
protein_anno <- protein_anno %>% dplyr::rename(c(protein = From, protein_gene = Gene.Names)) #rename so clear
protein_anno <- protein_anno %>% dplyr::filter(., protein %in% highpip$protein)
gene_list <- protein_anno$protein_gene
protein_anno$protein_gene <- gsub(" .*", "", protein_anno$protein_gene)
gene_list <- protein_anno$protein_gene #some have alternative names for each gene (manually checked with uniprot - see below)
gene_list <- gsub(" .*", "", gene_list) #remove name alternates


#set up BioMart database query
#below specifies to use build 37 (newer build now avaialble)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl, page = "structure") #lists outputs of a biomaRt query, page selects which page we want to see the outputs from
head(attributes)


#Extract required information for the genes
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"), #select attributes to be retrieved
            filters = "external_gene_name",  #filter to be used in query
            values = c(gene_list), #values of the filter (list of Protein Gene Names)
             mart = ensembl) #prepared Mart object setting up query

#_______________________________________________________________________________________________________
#Prepare data frame with information from BioMart on the protein genes

#create biomart data frame, with gene, start and end columns
biomart_df <- tss %>%
  dplyr::rename( c(protein_gene = external_gene_name,
            t_start = transcript_start,
            t_end = transcript_end, 
            chromosome_ensembl = chromosome_name))
head(biomart_df)

#check chromosome names 
biomart_df$chromosome_ensembl 


#separate into positive and negative stands
biomart_df_pos <- biomart_df1[biomart_df1$strand == 1,] 
biomart_df_neg <- biomart_df1[biomart_df1$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present
length(unique(biomart_df_pos$protein_gene))  
length(unique(biomart_df_neg$protein_gene)) 
length(unique(biomart_df1$protein_gene)) 

## create output data frames for genes on positive strand and those on negative strand to be r-bound later 
#out_pos == genes on positive strand 
out_pos <- matrix(nrow=length(unique(biomart_df_pos$protein_gene)), ncol=3)
colnames(out_pos) <- c("start", "end","chromosome")
rownames(out_pos) <- unique(biomart_df_pos$protein_gene) 
head(out_pos) #data frame with protein gene names on positive strand as row names 

##out_neg == genes on negative strand 
out_neg <- matrix(nrow=length(unique(biomart_df_neg$protein_gene)), ncol=3)
colnames(out_neg) <- c("start", "end","chromosome")
rownames(out_neg) <- unique(biomart_df_neg$protein_gene) 
head(out_neg) #data frame with protein gene names on negative strand as row names 


#_______________________________________________________________________________________________________
#Add data to prepared df 

#Add transcription start and end sites to prepared data frames
## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end 

for(protein_gene in unique(biomart_df_pos$protein_gene)) {
  tmp <- biomart_df_pos[which(biomart_df_pos$protein_gene==protein_gene), ] #puts data back into relevant gene's row
  out_pos[protein_gene,"start"] <-  min(tmp$t_start) #selects the minimum start position
  out_pos[protein_gene,"end"] <-  max(tmp$t_end) #selects the maximum end position
  out_pos[protein_gene,"chromosome"] <-  unique(tmp$chromosome_ensembl)[1] #inputs the chromosome for that gene name
}
head(out_pos)


## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end 
for(protein_gene in unique(biomart_df_neg$protein_gene)) {
  tmp <- biomart_df_neg[which(biomart_df_neg$protein_gene==protein_gene), ]
  out_neg[protein_gene,"start"] <-  max(tmp$t_end)                            
  out_neg[protein_gene,"end"] <-  min(tmp$t_start)
  out_neg[protein_gene,"chromosome"] <-  unique(tmp$chromosome_ensembl)[1]
}
head(out_neg)

## row bind the positive and negative strand dataframes to give one out_df 
out_df <- rbind(out_pos, out_neg)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference 
## set as data frames
out_df <- as.data.frame(out_df) 
out_df$protein_gene <- row.names(out_df) 
rownames(out_df) <- NULL
head(out_df) #this data frame now contains start, end, chromosome and GeneNames for the CpGs
fwrite(out_df, file = "/data_dir/ewas/bayesr/output/new_step2/annotated_protein_genes2.csv")

#Add to highPIP dataframe
out_df <- out_df %>% dplyr::rename(protein_chr = "chromosome")
out_df <- left_join(protein_anno, out_df, by = "protein_gene")
all_anno <- left_join(cpg_anno, out_df, by = "protein")
head(all_anno)


#Several proteins have no gene annotation
#Several protein genes have no tss etc
#Get list of proteins with no gene annotation
all_anno %>% is.na()
all_anno_2 <- all_anno %>% mutate(protein_gene = na_if(protein_gene, ""))
list_prot <- all_anno_2 %>% dplyr::filter(is.na(protein_gene))
unique(list_prot$protein) #2 proteins have not been allocated a gene name

#___________________________________________________________________________________________________________________
#Get list of protein genes with no transcript start site info 
list_gene <- all_anno_2 %>% dplyr::filter(is.na(start))
length(unique(list_gene$protein)) #affects 3 proteins
#From manual look up (uniprot.org) problem genes have additional labels


#___________________________________________________________________________________________________________________

#save out unannotated proteins for annotation by hand
unannotated <- all_anno %>%
  dplyr::filter(is.na(start))
fwrite(unannotated, file = "/data_dir/ewas/bayesr/output/new_step2/unannotated_protein_genes.csv")  


#_______________________________________________________________________________________________________
#Add Cis/Trans annotations for CpG sites

#add in chromosome lengths for the protein-related genes
chr_length <- read.csv("/ext_dir/chromosomes/CHR_Length.csv")
chr_length <- chr_length %>% dplyr::rename(protein_chr_length = "CHR_Length")
chr_length <- chr_length %>% dplyr::rename(protein_chr = "CHR")
chr_length$protein_chr <- as.character(chr_length$protein_chr)
all_anno <- left_join(all_anno, chr_length, by = "protein_chr")


#create new columns for cis/trans labels
#get distances from transcription start site (NB. those on different chromosomes are not relevant)
all_anno <- all_anno %>%
  mutate(., distance = as.numeric(start) - as.numeric(pos)) 

all_anno <- all_anno %>% #add column labelling as CIS or TRANS
  mutate(., CIS_TRANS = ifelse(protein_chr == probe_chr & abs(distance) <= 1e6, "Cis", "Trans"))
setDT(all_anno) 

#Count cis vs trans associations
table_cis_trans <- all_anno %>% dplyr::count(protein, CIS_TRANS, sort = FALSE)
table_cis_trans <- table_cis_trans %>% mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
all_anno2 <-  all_anno %>% mutate(across(CIS_TRANS, ~ifelse(.=="", NA, as.character(.))))
all_anno2 %>% dplyr::count(CIS_TRANS, sort = FALSE) 
#   CIS_TRANS     n
#       <char> <int>
# 1:       Cis   286
# 2:     Trans   387
# 3:      <NA>    24

#Save out final annotation file
fwrite(all_anno2, file = "/data_dir/ewas/bayesr/output/new_step2/cis_trans_anno_trial_1604.csv")
