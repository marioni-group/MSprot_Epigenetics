#Setting up required files for OSCA EWAS

#Load required packages
library(tidyverse)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(limma)
library(data.table)
library(here)
library(wrapr)
library(dplyr)
library(missMethods)

###--------------------------------------------------------------------------------
## Step 1 - read in DNA methylation file, prepare for OSCA and save as .txt file 
# Use file created from methylation regression (ewas_full_methylation_files_prep) 

data <- readRDS(here("input/methylation_data/GS_allchrom_resid_mvals_scaled.rds"))

#Set up for osca input file
data[1:5, 1:5]
osca_dat <- data.frame(FID=data$rn, IID=data$rn) #Create new data frame, adding IID and FID
osca_dat <- cbind(osca_dat, data) #This keeps Sample_Name and SSID
osca_dat <- osca_dat %>%
  select(-rn)
osca_dat[1:5, 1:5] #check

write.table(osca_dat, file=here("input/methylation_data/osca_data_1911.txt", row.names=F, quote=F, sep=' '))  

#--------------------------------------------------------------------------------
## Step 2 - in the command line (keep R open in a separate window) 
# Create binary methylation data (--methylation-beta if beta values, --methylation-m if M-values)

cd /data_dir/ewas/input_files #in terminal, check in correct folder
 osca_Linux \
 --efile methylation_data/osca_data_1911.txt \
 --methylation-m \
 --make-bod \
 --out osca_full 

#--------------------------------------------------------------------------------
## Step 3 - rewrite a blank .opi file created by OSCA in above step

# load EPIC annotation file 
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(anno)

#extract CpGs from methylation object 
opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]

# remove chr from chromosome names
opi$chr <- gsub("chr", "", opi$chr)

# replace chr X and Y with 23 and 24
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)

# set chr as numeric 
opi$chr <- as.numeric(opi$chr)

# set genes and strands as factors 
opi$UCSC_RefGene_Name	 <- as.factor(opi$UCSC_RefGene_Name	)
opi$strand <- as.factor(opi$strand)

# set missing gene names to NA 
opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA
head(opi)

# save file, overwriting .opi file generated in Step 2 
here()
write.table(opi, file="ewas/input/osca_full.opi",
	             col.names=F,
	             row.names=F,
	             quote=F, sep='\t') 

###--------------------------------------------------------------------------------
#Create phenotype files in batches (for parallel processing)

#set location for original phenotype files
ext_data <- "/data_dir/"

#First check order of IDs in methylation and protein files are the same
order <- read.csv(here("/input/methylation_data/order_resid_scaled_methylation_data.csv"))
pheno <- read.csv(here("ewas/input/protein_data/resid_proteins_basic_scaled.csv", check.names = F))

#Need to add SSID to phenotype file
target <- readRDS(paste0(ext_data, "GS/GS_methylation/GS20k/GS20k_Targets_18869.rds"))
target <- target %>%
  select(Sample_Name, Sample_Sentrix_ID) #select IDs in target file
rownames(target) <- c() #remove named rownames

target <- target %>%
filter(Sample_Name %in% pheno$id) #filter to IDs in phenotype(protein) file
identical(pheno$id, target$Sample_Name) #FALSE

names(target)[names(target) == "Sample_Name"] <- "id" #rename Sample_Name to id
pheno2 <- merge(pheno, target, by = "id") #merge files to add SSID to phenotype file

#move phenotype file to front of df
pheno2 <- pheno2 %>%
  mutate(FID = Sample_Sentrix_ID)
pheno2 <- pheno2 %>%
  select(Sample_Sentrix_ID, FID, everything())

#Match order of files
trial2 <- pheno2
trial2 <- trial2[match(osca_dat$FID, trial2$FID),] 
trial2 <- trial2 %>%
  select(-Sample_Sentrix_ID, -id, FID, everything()) #trial2 now has order-matched phenotype data
pheno <- trial2 #re-assign to pheno

# Check orders are now matched
identical(pheno$FID, order$order) # TRUE
pheno <- pheno %>%
  mutate(IID = FID) 
pheno <- pheno %>% 
  select(FID, IID, everything()) 

#-----------------------------------------------------------------------------------
#Write out protein data in batches for osca EWAS

#Test
location <- "protein_data/"

#Batch1
for(i in 3:113){ 
name <- as.character(names(pheno)[i])
write.table(select(pheno,c("FID","IID",name)), file = paste0(location, "batch1/",name,".phen"), col.names = T, row.names=F, quote = F, sep=' ')
}

#Batch2
for(i in 114:234){ 
name <- as.character(names(pheno)[i])
write.table(select(pheno,c("FID","IID",name)), file = paste0(location, "batch2/",name,".phen"), col.names = T, row.names=F, quote = F, sep=' ')
}

#Batch3
for(i in 235:335){ 
name <- as.character(names(pheno)[i])
write.table(select(pheno,c("FID","IID",name)), file = paste0(location, "batch3/",name,".phen"), col.names = T, row.names=F, quote = F, sep=' ')
}

#Batch4
for(i in 336:441){ 
name <- as.character(names(pheno)[i])
write.table(select(pheno,c("FID","IID",name)), file = paste0(location, "batch4/",name,".phen"), col.names = T, row.names=F, quote = F, sep=' ')
}

