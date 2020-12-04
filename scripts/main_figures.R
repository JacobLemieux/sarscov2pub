### Generate Main Figures 
# May 18 2020
# Updated November 2020
# lemieux@broadinstitute.org
# figure references correspond to original medrxiv submission (lemieux et al 2020, medrxiv)

# load packages
library(ggplot2)
library(cowplot)
library(ggtree)
library(treeio)
library(ape)
library(phytools)
library(tidytree)
library(reshape2)
library(MASS)
library(beastio)
library(lubridate)
library(tidyverse)
library(Biostrings)
library(gplots)
library(DECIPHER)
library(PopGenome)
library(ips)
library(devtools)
library(choroplethrZip)
library(dplyr)
library(gt)
library(webshot)
library(choroplethr)
library(choroplethrMaps)

setwd("COVID/SARS-CoV-2/")
source("scripts/SARS-CoV-2_functions.R")

###################################
# Set local SARS-CoV-2 directory 
###################################

###################################
# Compute or load some metadata
###################################

#system("pangolin data/sequences_and_trees/concatenated_unaligned.fasta -t 10 -p -o data/metadata --outfile Pangolin_lineage_report.csv")
#system("ln -s -f COVID/SARS-CoV-2/data/sequences_and_trees/mafft_aligned.fasta data/popgenome/mafft_aligned.fasta")

###################################
# Read in datasets
###################################

clinical_data <- read_REDCap_data("data/metadata/COVID19ViralSequenci_DATA_2020-07-27_2018.csv") # clinical data from REDCap

sample_set_metadata <- convert_REDCap_to_metadata(clinical_data) # convert to metadata format
dph_data <- read_DPH2("data/metadata/SabetiDeIdentifiedList_070220.csv") # read in DPH metadata
sample_set_metadata <- merge_metadata(dph_data, sample_set_metadata) # combine MGH and DPH metadata
sample_set_metadata <- left_join(sample_set_metadata, read_csv("data/metadata/zip_code_by_city_state.csv"), by = "zip_code_of_residence")
sample_set_metadata <- left_join(sample_set_metadata, read_csv("data/metadata/DPH_sample_county.csv"), by = "sample_id")

nosocomial_annotations <- data.frame(sample_id = c("MA_MGH_00902", "MA_MGH_00903", "MA_MGH_00905", "MA_MGH_00908", "MA_MGH_00227", "MA_MGH_00316"), nosocomial_cluster = c(rep("Cluster 2", 4), rep("Cluster 1", 2)))
sample_set_metadata <- left_join(sample_set_metadata, nosocomial_annotations, by = "sample_id") # combine cluster (eventually this info should be added to REDCap and this step deprecated)
sample_set_metadata$county[grep("DPH", sample_set_metadata$sample_id)] <- sample_set_metadata$DPH_county[grep("DPH", sample_set_metadata$sample_id)]

viral_quants <- read_viral_ct("data/metadata/CovidSeq_SeqPrepSheets_combined_RNAInfo_2020-06-16.csv") # Ct data (compiled by Katie)

pangolin_lineages <- read_pangolin("data/metadata/Pangolin_lineage_report.csv") # pangolin lineages

GENOME.class <- readData("data/popgenome") # alignment for pop gen stats

sample_set <- clean_sample_set(read_tsv("data/metadata/sample_set_entity_2020_06_28.tsv")) # load sample set (from Terra), remove water controls and missing lines

sample_set <- left_join(sample_set, clinical_data, by = "sample_id") # merge with patient identifiers

sample_set <- left_join(sample_set,viral_quants, by="sample_id")  %>% distinct(sample_id, .keep_all = TRUE) # add viral quants

sample_set <- left_join(sample_set, pangolin_lineages, by = "sample_id") # add pangolin lineages

sample_set <- left_join(sample_set, dph_data, by = "sample_id")

sample_set$record_id[grep("DPH", sample_set$sample_id)] <- sample_set$DPH_acc[grep("DPH", sample_set$sample_id)] # merge unique identifier for DPH samples
sample_set$date[grep("DPH", sample_set$sample_id)] <- sample_set$DPH_date[grep("DPH", sample_set$sample_id)] # merge unique identifier for DPH samples

sample_set_hq <- filter_sample_set(sample_set, assembly_percent_complete = 0.98) # filter for high quality genomes, based on completeness 

hq_serial_counts <- sample_set_hq %>% group_by(record_id) %>% summarise(n = n(), earliest_date = min(date) , sample_id = sample_id[date == min(date)][1]) # create table of counts by patient ID

sample_set_unique_hq <- left_join(hq_serial_counts , sample_set_hq, by = "sample_id")[,c("sample_id", "date", "earliest_date", "record_id.x", "assembly_fasta")] %>% distinct(sample_id, record_id.x, .keep_all = TRUE) %>% filter(!is.na(sample_id)) # uniqueify by patient ID and remove duplicates
case_stats <- read_and_parse_DPH_cases("data/metadata/Cases_05_31_2020.csv", sample_set_metadata, sample_set_unique_hq) # case data from DPH
out_of_state <- sample_set_metadata[!is.na(sample_set_metadata$zip_code_of_residence) & is.na(sample_set_metadata$county),]

out_of_state_hq <- out_of_state[out_of_state$sample_id %in% sample_set_unique_hq$sample_id,]

#correct for mis-identification

sample_set_metadata$exposure[sample_set_metadata$sample_id == "MA_MGH_00115"] <- "NoKnown"

# output county data by sample:
#write_tsv(sample_set_metadata[sample_set_metadata$sample_id %in% sample_set_unique_hq$sample_id ,c("sample_id", "county")], "data/metadata/sample_county.tsv")

###################################
# Align and trim
###################################

### download fasta sequences, align with mafft, trim UTRs, and write trimmed alignment
#fasta_seqs <- pull_fastas_and_concatenate(sample_set_unique_hq, clean_dir=TRUE, do_not_download = FALSE)
#mafft_aln <- mafft(as.DNAbin(fasta_seqs), thread = 10) # align and trim fastas
#write.fas(mafft_aln, "data/sequences_and_trees/mafft_aligned.fasta") # write alignment
mafft_alnDNAss <-readDNAMultipleAlignment("data/sequences_and_trees/mafft_aligned.fasta", format="fasta")
DNAss <- as(substr(mafft_alnDNAss, 268, ncol(mafft_alnDNAss) - 230), "DNAStringSet")

DNAss <- relabel_DNAss(DNAss, sample_set_metadata)
#writeXStringSet(DNAss, "data/sequences_and_trees/trimmed_alignment.fasta", format = "fasta")

########################
# Compute and load trees
########################

## Compute ML trees
#system("fasttree -nt -gtr -quote data/sequences_and_trees/trimmed_alignment.fasta > data/sequences_and_trees/trimmed_alignment.tree") # fasttree
#system("trimal -in data/sequences_and_trees/trimmed_alignment.fasta -out data/sequences_and_trees/trimmed_alignment.phy -phylip")
#system("phyml --no_memory_check -i data/sequences_and_trees/trimmed_alignment.phy") # phyml tree
#system("rm data/sequences_and_trees/trimmed_alignment.phy")
#system("iqtree -s data/sequences_and_trees/trimmed_alignment.fasta -bb 10000 -nt 10") # iqtree with 10000 UF bootstraps

# BEAST has to be run outside of R (at the moment).

# Construct metadata for iqtree
metadata_iqtree <- sample_set_metadata
metadata_iqtree$label <- gsub("\\|", "_", metadata_iqtree$label)
metadata_iqtree$exposure[metadata_iqtree$exposure == "NoKnown"] <- NA


# Read in root-to-tip distance (requires tempest, and usually takes as input PhyML tree)
root_to_tip <- read_tsv("data/sequences_and_trees/root-to-tip_trimmed_alignment_phyml.txt") # data file from tempEst for root-to-tip regression

# Read in nextstrain tree and metadata
#nextstrain_tree <- read.newick("data/sequences_and_trees/all_samples_aligned.fasta_aligned.filtered.masked_timetree.nwk")
nextstrain_tree <- read.newick("COVID/2020_08_07/nextstrain/all_samples_aligned.fasta_aligned.filtered.masked_timetree.nwk")
#nextstrain_mltree <- read.newick("data/sequences_and_trees/all_samples_aligned.fasta_aligned.filtered.masked_iqtree.nwk")
nextstrain_mltree <- read.newick("COVID/2020_08_07/nextstrain/all_samples_aligned.fasta_aligned.filtered.masked_iqtree.nwk")
#nextstrain_metadata <- read_tsv("data/metadata/output_gisaid_v2.tsv", col_types = cols(.default = "c")) # nextstrain metadata
nextstrain_metadata <- read_tsv("COVID/2020_08_07/nextstrain/output_gisaid_v3.tsv", col_types = cols(.default = "c")) # nextstrain metadata

# clean up nextstrain tree metadata
nextstrain_metadata$exposure <- rep(NA, nrow(nextstrain_metadata))
nextstrain_metadata$exposure[nextstrain_metadata$CONF_A_EXPOSURE == "YES"] <- "ConfA"
nextstrain_metadata$exposure[nextstrain_metadata$SNF_A_EXPOSURE == "YES"] <- "SNFA"
nextstrain_metadata$MA <- rep(NA, nrow(nextstrain_metadata))
nextstrain_metadata$MA[nextstrain_metadata$geoloc_cat == "Massachusetts"] <- "MA"
nextstrain_metadata$MGH_DPH <- rep(NA, nrow(nextstrain_metadata))
nextstrain_metadata$MGH_DPH[grep("MA_DPH", nextstrain_metadata$strain)] <- "DPH"
nextstrain_metadata$MGH_DPH[grep("MA_MGH", nextstrain_metadata$strain)] <- "MGH"

######################################
### identify and annotate clades
######################################

#ma_local_nodes <- c(8228, 8418, 7346, 7265, 5888)
#ma_local_nodes <- c(8125,8731,  7279, 7198, 6734) 
ma_local_nodes <- c( 8203,8787, 7354, 7012, 5875)  # set these from visual inspection of the tree
ggtree(nextstrain_tree) %<+% nextstrain_metadata +
  geom_highlight(node = ma_local_nodes[1], fill = "firebrick") +
  geom_highlight(node = ma_local_nodes[2], fill = "lightblue") +
  geom_highlight(node = ma_local_nodes[3], fill = "darkgreen") +
  geom_highlight(node = ma_local_nodes[4], fill = "purple") +
  geom_highlight(node = ma_local_nodes[5], fill = "yellow") +
  geom_tippoint(aes(color=MA, subset = !is.na(MA)), size = 6) +
  geom_text2(aes(subset=!isTip, label=node), size = 8)
ggsave("data/convenience_plots/full_nextstrain_tree_with_node_labels.pdf", height = 1000, width = 50, limitsize = FALSE)

ggtree(nextstrain_mltree) %<+% nextstrain_metadata +
  geom_tippoint(aes(color=MA, subset = !is.na(MA)), size = 6) +
  geom_text2(aes(subset=!isTip, label=node), size = 8) +
  geom_text2(aes(subset=isTip, label = label), size = 4, nudge_x = 0.001)
ggsave("data/convenience_plots/full_nextstrain_MLtree_with_node_labels.pdf", height = 1000, width = 50, limitsize = FALSE)

sample_set_unique_hq$clade <- rep(NA, nrow(sample_set_unique_hq))
for(i in 1:length(ma_local_nodes)){
  descendants_tmp <- nextstrain_tree$tip.label[getDescendants(nextstrain_tree, node = ma_local_nodes[i])]
  descendants_tmp <- descendants_tmp[!is.na(descendants_tmp)]
  descendants_tmp <- descendants_tmp[grep("MA_", descendants_tmp)]
  sample_set_unique_hq$clade[sample_set_unique_hq$sample_id %in% descendants_tmp] <- paste("Clade_", i, sep="")
}
clade_set <- sample_set_unique_hq[,c("sample_id", "clade")]
sample_set_metadata <- left_join(sample_set_metadata, clade_set, by = "sample_id")
sample_set_metadata$clade[is.na(sample_set_metadata$clade)] <- "None"
sample_set_metadata <- unite(sample_set_metadata, label, c(sample_id, exposure, clade, date), sep="|", remove=FALSE)
sample_set_metadata <- sample_set_metadata[,c("sample_id", "label", "date", "exposure", "clade", "nosocomial_cluster", "zip_code_of_residence", "clinic", "county")]

DNAss2 <- as(substr(mafft_alnDNAss, 268, ncol(mafft_alnDNAss) - 230), "DNAStringSet")
DNAss2 <- relabel_DNAss(DNAss2, sample_set_metadata)
#writeXStringSet(DNAss2, "data/sequences_and_trees/trimmed_alignment_with_clades.fasta", format = "fasta")

# Read in iqtree trees and modify metadata to have iqtree name format
iqtree <- read.iqtree("data/sequences_and_trees/trimmed_alignment.fasta.contree")
iqtree@phylo <- midpoint.root(iqtree@phylo) # midpoint root ML tree
metadata_iqtree_clade_labels <- unite(sample_set_metadata, label, c(sample_id, exposure, clade, date), sep="_", remove = FALSE)


###################################
# Main Figures
###################################

# choose color palette:
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
# To use for fills, add
#scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

## Figure 1
case_stats$Date <- as_date(case_stats$Date)
cases_by_site <- case_stats[!case_stats$Site == "Proportion",]
case_stats$Site <- gsub("Proportion", "Cumulative Sampling Proportion", case_stats$Site)
p1 <- ggplot(data = subset(cases_by_site, cases_by_site$Date < "2020-05-01"), aes(x = Date, y = Cases, col = Site)) + 
  geom_point() + 
  theme_bw() + 
  stat_smooth(span = 0.3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))+ 
  xlab("Date") + 
  ylab("Cumulative Cases") +
  theme(legend.text = element_text(size=14)) + 
  scale_y_log10() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) +
  scale_colour_manual(values=cbPalette)
p1
ggsave("data/figures/F1A.pdf", width = 5, height = 5)
ggsave("data/figures/F1A.jpg", width = 5, height = 5)

case_props <- case_stats[case_stats$Site == "Cumulative Sampling Proportion",]

p2 <- ggplot(data = subset(case_props, case_props$Date < "2020-05-01"), aes(x = Date, y = Cases, col = Site)) + 
  geom_point() + 
  stat_smooth() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))+ 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  xlab("Date") + 
  ylab("Cumulative Sampling Fraction") +
  theme(legend.text = element_text(size=14)) + 
  scale_y_log10() + 
  theme(legend.position = c(0.6, 0.915), legend.title=element_blank()) + 
  scale_colour_manual(values=cbPalette)
p2
ggsave("data/figures/F1B.pdf", width = 5, height = 5)
ggsave("data/figures/F1B.jpg", width = 5, height = 5)

# new figure 1C
test_data <- read_csv("COVID/COVID_Cumulative_Results_8_4_2020.csv", skip = 1)
test_data$`Collect Date` <- as_date(mdy(test_data$`Collect Date`))
test_data$Birthdate <- as_date(mdy(test_data$`Birthdate`))
test_data <- test_data %>% filter(!is.na(Birthdate))
test_data$Age <- age_calc(test_data$Birthdate, unit = "years")
test_data <- test_data %>% filter(`Collect Date` > "2020-01-01")



test_pos <- test_data %>% filter(`Result 1` == "Detected" | `Result 1` == "Positive" | `Result 1` == "Presumptive Positive" | `Result 1` == "Positive for 2019-novel Coronavirus (2019-nCoV) by PCR."| `Result 1` =="Presumptive positive for 2019-novel Coronavirus (2019-nCoV) by PCR." | `Result 1` == "Positive for SARS-Cov-2 (COVID-19) by PCR")
test_pos$res <- "Positive"
test_neg <- test_data %>% filter(`Result 1` == "NOT DETECTED" | `Result 1` == "Negative"| `Result 1` == "Not Detected"  | `Result 1` == "Negative for SARS-Cov-2 (COVID-19) by PCR" | `Result 1` == "Negative for 2019-novel Coronavirus (2019-nCoV) by PCR."| `Result 1` =="Negative for SARS-Cov-2 (COVID-19) by PCR")
test_neg$res <- "Negative"

test_patients <- rbind(test_pos, test_neg)

test_pos_patients <-  test_pos %>% group_by(PtNumber) %>% summarise(n = n(), earliest = min(`Collect Date`), latest = max(`Collect Date`), duration = max(`Collect Date`) -  min(`Collect Date`))
test_pos_patients_by_date <- test_pos_patients %>% group_by(earliest) %>% summarise(n = sum(n))

ma_cases <- read_csv("data/metadata/Cases_05_31_2020.csv")
ma_june_cases <- read_csv("data/metadata/Cases_06_17_2020.csv")
ma_june_cases <- ma_june_cases[,c("Date", "Positive New")]
names(ma_june_cases) <- c("Date", "New")
ma_cases <- ma_cases[,c("Date", "New")]
ma_cases <- rbind(ma_cases, ma_june_cases)

pcr_pos_date <- test_pos_patients_by_date[,c("earliest", "n")]
names(pcr_pos_date) <- c("Date", "New")

all_cases <- rbind(ma_cases, pcr_pos_date)
all_cases$Type <- c(rep("Reported", nrow(ma_cases)), rep("MGH", nrow(pcr_pos_date)))
#data = subset(all_cases, all_cases$Date > "2020-03-20")
p1 <- ggplot(subset(all_cases, all_cases$Date < "2020-06-15"), aes(x = Date, y = New, col = Type)) + 
  geom_point() + 
  geom_smooth() + 
  theme_bw() + 
  xlab("Collection Date") + ylab("N") + 
  #  theme(axis.text=element_text(size=20), axis.title=element_text(size = 20))+ 
  theme(legend.text = element_text(size=14)) + 
  theme(legend.position = c(0.82, 0.915), legend.title=element_blank()) + 
  scale_y_log10() +  
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  scale_colour_manual(values=cbPalette)
p1
ggsave("data/figures/F1C.pdf",width=5,height=5)
ggsave("data/figures/F1C.jpg",width=5,height=5)


# Plot zip code choropleths

dem_data <- read_delim("COVID/061220130343720223.txt", delim="|", col_types = cols(.default = "c"))

# create df of zip code count data
#zip_df <- data.frame(table(dem_data$Zip_code))
zip_df <- data.frame(table(sample_set_metadata$zip_code_of_residence[sample_set_metadata$sample_id %in% sample_set_unique_hq$sample_id & sample_set_metadata$exposure != "SNFA"]))
zip_df <- zip_df[!zip_df$Var1 == "@",]
write_tsv(zip_df,"data/metadata/zip_table.tsv")
zip_df <- read_tsv("data/metadata/zip_table.tsv")
names(zip_df) <- c("region", "value1")
zip_df$region <- as.character(zip_df$region)
zip_df$value1 <- as.numeric(zip_df$value1)

# plot map with zip code info
data(df_pop_zip)
data(df_zip_demographics)
ma_pop_zip <- df_pop_zip
ma_pop_zip <- left_join(ma_pop_zip, zip_df, by = "region")
ma_pop_zip2 <- df_pop_zip
ma_pop_zip2$value <- 0
ma_pop_zip2 <- left_join(ma_pop_zip2, zip_df, by = "region")
ma_pop_zip2$value1[is.na(ma_pop_zip2$value1)] <- 0
ma_pop_zip2 <- ma_pop_zip2[,c("region", "value1")]
ma_pop_ziplog <- data.frame(region = ma_pop_zip2$region, value = log10(ma_pop_zip2$value + 1))
names(ma_pop_zip2)[2] <- "value"

p1 <- zip_choropleth(ma_pop_ziplog, state_zoom = "massachusetts", num_colors = 1)
p1
ggsave("data/figures/F1E.pdf")
ggsave("data/figures/F1E.jpg")

p2 <- zip_choropleth(ma_pop_ziplog, county_zoom=c(25017,25025), num_colors = 1)
p2
ggsave("data/figures/F1F.pdf")
ggsave("data/figures/F1F.jpg")

# county choropleth
data(df_pop_county)
data(county.regions)
county_choropleth(df_pop_county, state_zoom = "massachusetts")

# merge with 2013 ACS survey data
zip_dem <- left_join(zip_df, df_zip_demographics, by = "region")

p5 <- ggplot(subset(zip_dem, zip_dem$value1 > 5), aes(x = value1/total_population, y = percent_white)) + geom_point() + geom_smooth(method = "lm")
p6 <- ggplot(subset(zip_dem, zip_dem$value1 > 5), aes(x = value1/total_population, y = percent_hispanic)) + geom_point() + geom_smooth(method = "lm")
p7 <- ggplot(subset(zip_dem, zip_dem$value1 > 5), aes(x = value1/total_population, y = percent_black)) + geom_point() + geom_smooth(method = "lm")
p8 <- ggplot(subset(zip_dem, zip_dem$value1 > 5), aes(x = value1/total_population, y = per_capita_income)) + geom_point() + geom_smooth(method = "lm")
plot_grid(p5,p6, p7, p8)
ggsave("data/figures/demographics_by_zip.pdf")
write_tsv(zip_dem, "data/metadata/demographics_by_zip_2013_ACS.tsv")

# aggregate by county

county_df <- data.frame(table(sample_set_metadata$county[sample_set_metadata$sample_id %in% sample_set_unique_hq$sample_id]))
county_df$Var1 <- tolower(county_df$Var1)
names(county_df) <- c("county.name", "value1")
county_df <- left_join(county_df, county.regions %>% filter(state.name == "massachusetts"))
county_df$value1 <- as.numeric(county_df$value1)
county_df <- county_df[,c("region", "value1")]

# plot map with country info

ma_pop_county <- left_join(df_pop_county, county_df, by = "region")
ma_pop_county$value1[is.na(ma_pop_county$value1)] <- 0
ma_pop_county <- ma_pop_county[,c("region", "value1")]
names(ma_pop_county)[2] <- "value"
ma_pop_county_log <- data.frame(region = ma_pop_county$region, value = log10(ma_pop_county$value + 1))

p1 <- county_choropleth(ma_pop_county_log, state_zoom = "massachusetts", num_colors = 1) + 
  scale_fill_gradient(labels = round(10^c(seq(0, 2.5, by = 0.5))-1, 0), high = "#0047ab", low = "white") + 
  theme(legend.title = element_blank())
p1

ggsave("data/figures/F1E.pdf", height =5, width = 5)
ggsave("data/figures/F1E.jpg", height =5, width = 5)


p2 <- zip_choropleth(ma_pop_ziplog, county_zoom=c(25017,25025), num_colors = 1) + 
  scale_fill_gradient(labels = round(10^c(seq(0, 2, by = 0.5))-1, 0), high = "#0047ab", low = "white") + 
  theme(legend.title = element_blank())
p2
ggsave("data/figures/F1F.pdf", height =5, width = 5)
ggsave("data/figures/F1F.jpg", height =5, width = 5)

# compare to total DPH counts

dph_county <- read_csv("data/metadata/County_8_11_2020.csv")
dph_county$Date <- mdy(dph_county$Date)
dph_county <- dph_county %>% filter(Date == "2020-07-01")

sample_set_county <- data.frame(table(sample_set_metadata$county[sample_set_metadata$sample_id %in% sample_set_unique_hq$sample_id]))
names(sample_set_county) <- c("County", "Study Count")

county_counts <- left_join(sample_set_county, dph_county, by = "County")
p3 <- ggplot(county_counts, aes(x = Count, y = `Study Count`, label = County)) + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw()  + geom_label() 
#ggsave("data/figures/DPH_counts_vs_study_counts_by_county.pdf", height =5, width = 5)
#ggsave("data/figures/DPH_counts_vs_study_counts_by_county.jpg", height =5, width = 5)


# make map of DPH counts
dph_county_df <- dph_county[c("County", "Count")]
dph_county_df$County <- tolower(dph_county_df$County)
names(dph_county_df) <- c("county.name", "value1")
dph_county_df <- left_join(dph_county_df, county.regions %>% filter(state.name == "massachusetts"))
dph_county_df$value1 <- as.numeric(dph_county_df$value1)
dph_county_df <- dph_county_df[,c("region", "value1")]

# plot map with country info

dph_pop_county <- left_join(df_pop_county, dph_county_df, by = "region")
dph_pop_county$value1[is.na(dph_pop_county$value1)] <- 0
dph_pop_county <- dph_pop_county[,c("region", "value1")]
names(dph_pop_county)[2] <- "value"

dph_pop_county_log <- data.frame(region = dph_pop_county$region, value = log10(dph_pop_county$value + 1))

# plot DPH case counts
p2 <- county_choropleth(dph_pop_county_log, state_zoom = "massachusetts", num_colors = 1) + 
  scale_fill_gradient(labels = formatC(round(10^c(seq(2, 5, by = 1)), 0), format="e", digits = 1), high = "#0047ab", low = "white") + 
  theme(legend.title = element_blank())


#ggsave("data/figures/DPH_counts_July_1_by_county.pdf", height =5, width = 5)
#ggsave("data/figures/DPH_counts_July_1_by_county.jpg", height =5, width = 5)

# plot sampling proportion
sampling_prop <- data.frame(region = dph_pop_county$region, value = ma_pop_county$value / dph_pop_county$value)
sampling_prop$value[is.nan(sampling_prop$value)] <- 0
p4 <- county_choropleth(sampling_prop, state_zoom = "massachusetts", num_colors = 1) 
#ggsave("data/figures/sampling_proportion_by_county.pdf", height =5, width = 5)
#ggsave("data/figures/sampling_proportion_by_county.jpg", height =5, width = 5)
plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), label_size = 20)
ggsave("data/figures/S3.pdf", height =8, width = 10)
ggsave("data/figures/S3.jpg", height =8, width = 10)


# 
dph_counts <- read_csv("data/metadata/Cases_05_31_2020.csv")
names(dph_counts)[5] <- "Statewide"
p6 <- ggplot(subset(dph_counts, dph_counts$Date < "2020-05-11"), aes(x = Date, y = Statewide)) + geom_point() + geom_smooth(color = cbPalette[1], span = 0.7) + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank()) + 
  scale_color_manual(values=cbPalette[1])

p6  




## Figure 2 - allele frequency over time for particular mutations

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u")]
var_tbl$variant <- rep("C2416T (Conf, Homeless)", nrow(var_tbl))
var_tbl2 <- var_tbl
var_tbl2$variant <- gsub("(Conf, Homeless)", "Conference", var_tbl2$variant)
p7 <- ggplot(var_tbl2, aes(x = date, y = af, col = variant)) + geom_point(size = 2) + geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank())+
  scale_color_manual(values=cbPalette[2])

#ggsave("data/figures/C2416T.pdf", height = 5, width = 5)
#ggsave("data/figures/C2416T.jpg", height = 5, width = 5)

plot_grid(p7,p6, nrow = 2)
#ggsave("data/figures/C2416T_vs_epidemic.pdf", height = 5, width = 5)
#ggsave("data/figures/C2416T_vs_epidemic.jpg", height = 5, width = 5)

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,105], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- rbind(var_tbl,data.frame(date = var_by_date$date, af = var_by_date$af, ci_l = var_by_date$ci_l, ci_u = var_by_date$ci_u, variant = rep("G105T (Homeless)", nrow(var_by_date))))

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,23403], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "A"), count_derived = sum(variant == "G"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- rbind(var_tbl,data.frame(date = var_by_date$date, af = var_by_date$af, ci_l = var_by_date$ci_l, ci_u = var_by_date$ci_u, variant = rep("A23403G (D614G)", nrow(var_by_date))))

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- rbind(var_tbl,data.frame(date = var_by_date$date, af = var_by_date$af, ci_l = var_by_date$ci_l, ci_u = var_by_date$ci_u, variant = rep("G26233T (Conf, Homeless)", nrow(var_by_date))))

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,3892], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- rbind(var_tbl,data.frame(date = var_by_date$date, af = var_by_date$af, ci_l = var_by_date$ci_l, ci_u = var_by_date$ci_u, variant = rep("G3892T (SNF)", nrow(var_by_date))))

p2 <- ggplot(var_tbl, aes(x = date, y = af, col = variant)) + geom_point(size = 2) + geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.72, 0.68), legend.text=element_text(size = 12), legend.title = element_blank(), legend.background = element_blank())+
  scale_color_manual(values=cbPalette[c(4,2,3,1,5)])
p2
ggsave("data/figures/F2B.pdf", height = 5, width = 5)
ggsave("data/figures/F2B.jpg", height = 5, width = 5)

## Figure 3 
#beast <- read.beast("data/sequences_and_trees/MCC.tree") # MCC tree from BEAST
beast <- read.beast("data/sequences_and_trees/MCC_GTR4_2020_08_07.tree")
#mcmc_trace <- readLog("data/sequences_and_trees/trimmed_alignment_with_clades.log", burnin = 0.3) # BEAST log file
mcmc_trace <- readLog("COVID/2020_08_07/BEAST/GTRG4_coal_exp_pop/trimmed_alignment_with_clades.log", burnin = 0.3)
mcmc_df <- data.frame(mcmc_trace)
mcmc_df$iteration <- rownames(mcmc_df)

names(clade_set)[1] <- "strain"
nextstrain_metadata <- left_join(nextstrain_metadata, clade_set, by="strain")
nextstrain_metadata$clade <- gsub("Clade_", "Clade ", nextstrain_metadata$clade)


sample_set_metadata$clade[sample_set_metadata$clade == "None"] <- NA
sample_set_metadata$clade <- gsub("Clade_", "Clade ", sample_set_metadata$clade)
sample_set_metadata$clade <- gsub("Clade 1", "C20099T (BHCHP)", sample_set_metadata$clade)
sample_set_metadata$clade <- gsub("Clade 2", "G3892T (SNF)", sample_set_metadata$clade)
sample_set_metadata$clade <- gsub("Clade 3", "C2416T (Conf, BHCHP)", sample_set_metadata$clade)
sample_set_metadata$clade <- gsub("Clade 4", "G105T (BHCHP)", sample_set_metadata$clade)
sample_set_metadata$clade <- gsub("Clade 5", "G28899T", sample_set_metadata$clade)

sample_set_metadata_beast <- sample_set_metadata[,c("label", "exposure", "clade", "nosocomial_cluster", "zip_code_of_residence", "clinic")]


# clean up data to plot
MRCA_df <- mcmc_df[,c("mrca.date.Full_Set.","mrca.date.Clade_1.","mrca.date.Clade_2.","mrca.date.Clade_3.","mrca.date.Clade_4.","mrca.date.Clade_5.","iteration")]
names(MRCA_df) <- c("Root","C20099T (BHCHP)", "G3892T (SNF)" , "C2416T (Conf, BHCHP)","G105T (BHCHP)", "G28899T", "iteration")
MRCA_reshaped <- melt(MRCA_df, id.vars = "iteration")
MRCA_reshaped$date <- date_decimal(MRCA_reshaped$value)
sample_set_metadata_beast$exposure_clean <- rep(NA, nrow(sample_set_metadata_beast))
sample_set_metadata_beast$exposure_clean[sample_set_metadata$exposure == "ConfA"] <- "Conference" 
sample_set_metadata_beast$exposure_clean[sample_set_metadata$exposure == "SNFA"] <- "SNF" 
sample_set_metadata_beast$exposure_clean[sample_set_metadata$exposure == "Homeless"] <- "BHCHP"

p1 <- ggtree(beast, mrsd="2020-05-09") %<+% sample_set_metadata_beast +
  theme_tree2() +
  geom_tippoint(aes(color=exposure_clean, subset=!is.na(exposure_clean)), size = 2) +
  geom_text2(aes(subset=!isTip & as.numeric(posterior) > 0.8, label=round(posterior, 2)), hjust = 1.5, vjust = -0.3) +
  theme(legend.position = c(0.2, 0.85), legend.title = element_blank(), legend.text=element_text(size = 18)) +
#  geom_cladelabel(node=785, label="", fill = "firebrick", offset = .1) +
#  geom_hilight(node=1277, fill=cbPalette[3], alpha=.3) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size = 20)) +
  scale_color_manual(values=c("firebrick", "purple", "gray")) +
  guides(colour = guide_legend(override.aes = list(size=5)))
p1

#ggsave("data/figures/exposure_groups.pdf", width = 12, height = 8)
#ggsave("data/figures/exposure_groups.jpg", width = 12, height = 8)

p2 <- ggtree(beast, mrsd="2020-05-09") %<+% sample_set_metadata_beast + 
  theme_tree2() + #+ 
  geom_tippoint(aes(color=clade, subset=!is.na(clade)), size = 2) +
  geom_text2(aes(subset=!isTip & as.numeric(posterior) > 0.8, label=round(posterior, 2)), hjust = 1.5, vjust = -0.3) + 
  theme(legend.position = c(0.2, 0.85), legend.title = element_blank(), legend.text=element_text(size = 14)) + 
  #  theme(legend.position="none") + 
  scale_color_manual(values=cbPalette[c(1,3,5,2,4)]) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size = 20))

p2

p3 <- ggplot(MRCA_reshaped, aes(x = date, y= variable, fill = variable)) +
  geom_violin(adjust = 2, draw_quantiles = c(0.5)) +
  theme_classic() +
#  theme(axis.text=element_text(size=15, angle = 90), axis.title=element_text(size = 20)) + xlab("") + ylab("Date") + 
  theme(axis.text = element_text( hjust = 0.5, size = 14), axis.title=element_text(size = 20)) + 
  xlab("Date") + ylab("") +
  theme(legend.position = "none") + 
  scale_fill_manual(values=cbPalette[c(6,1:5)]) + scale_x_datetime(date_labels = "%b %d %Y", date_breaks = "2 months")
p3

plot_grid(p1, p2, labels = c("A", "B"), nrow = 1 , label_size = 20)
ggsave("data/figures/F3A-B.pdf", width = 12, height = 8)
ggsave("data/figures/F3A-B.jpg", width = 12, height = 8)

p3
plot_grid(p3, labels = c("C"), nrow = 1 , label_size = 20)

ggsave("data/figures/F3C.pdf", width = 12, height = 4)
ggsave("data/figures/F3C.jpg", width = 12, height = 4)

# Create summary table of dates: 
tmrca_summary <- apply(MRCA_df[,1:6], 2, function(x) quantile(x, c(0.05, 0.5, 0.95)))
tmrca_dates <- matrix(as.character(round_date(date_decimal(tmrca_summary), unit = "day")), nrow = nrow(tmrca_summary), dimnames = list(rownames(tmrca_summary),colnames(tmrca_summary)))
for(i in 1:ncol(tmrca_dates)){
  tmrca_dates[2,i] <- paste(tmrca_dates[2,i], " (", tmrca_dates[1,i], " - ", tmrca_dates[3,i], ") ", sep="")
}

mutations <- c("C20099T", "G3892T", "C2416T", "G105T", "G28899T")
colnames(tmrca_dates)[2:6] <- mutations
row_label <- c("Number of Genomes","Epidemiology", "Amino Acid substitution" ,"Median tMRCA (95% HPD)")
num_genomes <- c(772, 21, 77, 288, 98, 34)
epi_events <- c("", "BHCHP", "SNF", "Conference, BHCHP", "BHCHP", "")
aa_change <- c("", "ORF1b: A2211V; NSP15: A160V",  "ORF1a: E1209D; NSP3: E391D", "", "", "N: R56I, ORF14: E56*")
tmrca_table <- rbind(num_genomes, epi_events, aa_change, tmrca_dates[2,])
tmrca_table <- cbind(row_label, tmrca_table)
colnames(tmrca_table)[1] <- "Lineage"
t1 <- gt(data.frame(tmrca_table)) %>% tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels(1:7))
t1 
gtsave(t1, "data/tables/table_2.pdf")
# Outplot figure with nodes labeled
ggtree(beast, mrsd="2020-05-09") %<+% sample_set_metadata_beast +
  theme_tree2() +
  geom_tippoint(aes(color=exposure, subset=!is.na(exposure)), size = 2) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  scale_fill_manual(cbPalette)
ggsave("data/convenience_plots/beast_tree_with_nodes_labeled.pdf", height = 100, width = 50, limitsize=FALSE)

pfoo <- ggtree(beast, mrsd="2020-05-09") %<+% left_join(sample_set_metadata, sample_set[,c("sample_id", "Batch")], by = "sample_id")[,c("label", "Batch")] +
  theme_tree2() + #+ 
  geom_tippoint(aes(color=ifelse(Batch == "Batch11", "Batch11", "Other_Batches"), subset=!is.na(Batch)), size = 2) +
#  geom_text2(aes(subset=!isTip & as.numeric(posterior) > 0.8, label=round(posterior, 2)), hjust = 1.5, vjust = -0.3) + 
  theme(legend.position = c(0.2, 0.85), legend.title = element_blank(), legend.text=element_text(size = 14)) + 
  #  theme(legend.position="none") + 
#  scale_color_manual(values=cbPalette[c(1,3,5,2,4)]) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size = 20))

pfoo

ggsave("data/figures/Tree_with_batch11.pdf")
###################################
# Generate Summary Stats
###################################

# Number of unique sequences:
print(paste("Number of sequences:", nrow(sample_set)))
# Number of unique patients:
patient_serial_counts <- sample_set %>% group_by(record_id) %>% summarise(n = n())
print(paste("Number of unique patients:", nrow(patient_serial_counts)-1))
# Number of unique patients in high quality genomes
print(paste("Number of unique patients with high quality genomes:", nrow(sample_set_unique_hq)))
# print out summary stats for sequencing success
threshold <- 0.98
threshold2 <- 0.8
reference_length = 29903
num_hq <- sum(sample_set$assembly_length_unambiguous > (threshold*reference_length), na.rm=T)
fraq_hq <- round(sum(sample_set$assembly_length_unambiguous > (threshold*reference_length), na.rm=T)/nrow(sample_set),2)
num_hq2 <- sum(sample_set$assembly_length_unambiguous > (threshold2*reference_length), na.rm=T)
fraq_hq2 <- round(sum(sample_set$assembly_length_unambiguous > (threshold2*reference_length), na.rm=T)/nrow(sample_set),2)
print(paste("High quality assemblies: ", num_hq, sep=""))
print(paste("Fraction high quality assemblies: ", fraq_hq, sep=""))
print(paste("Good quality assemblies: ", num_hq2, sep=""))
print(paste("Fraction good quality assemblies: ", fraq_hq2, sep=""))
seq_stats <- data.frame(Min_frac_complete = c(threshold, threshold2), Number_genomes = c(num_hq, num_hq2), Proportion_complete = c(fraq_hq, fraq_hq2))
t1 <- seq_stats %>% gt()
t1
gtsave(t1, "data/figures/T0.pdf")

########################
# Supplemental Figures
########################

# Count the number of high quality assemblies by genome threshold
reference_length = 29903
percent_complete <- seq(1, 99, by = 1)
yield_tbl <- as_tibble(percent_complete)
yield_tbl$assembly_number <- rep(NA, nrow(yield_tbl))
for(i in 1:length(percent_complete)){
  yield_tbl$assembly_number[i] <- sum(sample_set$assembly_length_unambiguous >= ((percent_complete[i]/100) * reference_length))
}
names(yield_tbl) <- c("Percent Complete", "Number of Assemblies")
t1 <- yield_tbl %>% gt()
gtsave(t1, "data/figures/ST1.pdf")

#cor.test(sample_set$ViralCT, log10(sample_set$assembly_mean_coverage), use = "complete.obs", "spearman", alternative = "two.sided")
cor.test(sample_set$ViralCT, log10(sample_set$assembly_mean_coverage), use = "complete.obs", "pearson", alternative = "two.sided")

## Supplemental Figure 1

p1 <- ggplot(sample_set, aes(x = ViralCT, y = assembly_mean_coverage)) + 
  scale_y_log10() +  scale_x_log10() +
  geom_point() + 
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size = 16))+ 
  xlab("Viral Ct")+ ylab("Mean Coverage") + 
  geom_smooth(method=MASS::rlm)
#p1
p2 <- ggplot(sample_set, aes(x = ViralCT, y = assembly_length_unambiguous/29903))+ geom_point() + 
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size = 16))+ 
  xlab("Viral Ct")+ ylab("Fraction Unambiguous") + 
  stat_smooth(method="glm", method.args = list(family = "binomial"), se=TRUE) 
p3 <- ggplot(yield_tbl, aes(x = `Percent Complete`, y = `Number of Assemblies`)) + geom_point() + 
  theme_bw() + 
  theme(axis.text=element_text(size=16), axis.title=element_text(size = 16))
p4 <- ggplot(sample_set, aes(x=assembly_length_unambiguous/29903)) + geom_histogram() + 
  theme_bw()  + 
  xlab("Percent Complete") + ylab("Number of Assemblies") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size = 16))
#pdf("data/figures/Supplemental_figure_1.pdf")
plot_grid(p1, p2,p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 20)
ggsave("data/figures/S1.pdf")
ggsave("data/figures/S1.jpg")

# Supplemental Figure 2

RC_clean <- read_csv("COVID/SARS-CoV-2/data/metadata/cobas_ct_by_sample.csv")
RC_clean$ViralCT[RC_clean$ViralCT == 45] <- NA
RC_clean <- left_join(RC_clean, sample_set[,c("sample_id", "assembly_mean_coverage")])

p1 <- ggplot(RC_clean, aes(x = ViralCT, y = `Target 1 SARS-CoV-2`)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p2 <- ggplot(RC_clean, aes(x = ViralCT, y = `Target 2  Pan-Sabecovirus`))+ 
  geom_point() + 
  theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p3 <- ggplot(RC_clean, aes(x = `Target 1 SARS-CoV-2`, y = `Target 2  Pan-Sabecovirus`)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p4 <- ggplot(data = subset(sample_set, !is.na(DPH_N1)), aes(x = ViralCT, y = DPH_N1)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p5 <- ggplot(data = subset(sample_set, !is.na(DPH_N1)), aes(x = ViralCT, y = DPH_N2)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p6 <- ggplot(data = subset(sample_set, !is.na(DPH_N1)), aes(x = DPH_N1, y = DPH_N2)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p7 <- ggplot(RC_clean, aes(x = assembly_mean_coverage, y = `Target 1 SARS-CoV-2`)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_log10() + 
  xlab("Mean Coverage") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p8 <- ggplot(RC_clean, aes(x = assembly_mean_coverage, y = `Target 2  Pan-Sabecovirus`))+ 
  geom_point() + 
  theme_classic() + 
  scale_x_log10() + 
  xlab("Mean Coverage") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p9 <- ggplot(RC_clean, aes(x = assembly_mean_coverage, y = ViralCT)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_log10() + 
  xlab("Mean Coverage") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p10 <- ggplot(data = subset(sample_set, !is.na(DPH_N1)), aes(x = assembly_mean_coverage, y = DPH_N1)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_log10() + 
  xlab("Mean Coverage") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p11 <- ggplot(data = subset(sample_set, !is.na(DPH_N1)), aes(x = assembly_mean_coverage, y = DPH_N2)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_log10() + 
  xlab("Mean Coverage") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

p12 <- ggplot(data = subset(sample_set, !is.na(DPH_N1)), aes(x = assembly_mean_coverage, y = ViralCT)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_log10() + 
  xlab("Mean Coverage") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))

plot_grid(p1, p2, p3, p4, p5, p6,p7, p8, p9, p10, p11, p12, labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), nrow=2)
ggsave("data/figures/S2.pdf", width = 18)
ggsave("data/figures/S2.jpg", width = 18)

## Supplemental Figure 3

# x.dist <- DistanceMatrix(DNAss)*width(DNAss)
# pdf("data/figures/S3A.pdf", height = 8, width = 8)
# heatmap.2(log10(x.dist+1), trace="none", col="bluered", labRow=NA, labCol=NA)
# dev.off()
# 
# 
# pairwise_differences <- rep(0, 0)
# for(i in 1:nrow(x.dist)){
#   for(j in 1:nrow(x.dist)){
#     if(i < j){
#       pairwise_differences <- c(pairwise_differences, x.dist[i,j])
#     }
#   }
# }
# pairdiff <- data.frame(pairwise_differences = pairwise_differences)
# print(median(pairdiff$pairwise_differences))
# print(quantile(pairdiff$pairwise_differences, probs = c(0.25, 0.75)))
# print(range(pairdiff$pairwise_differences))
# 
# ggplot(pairdiff, aes(x = pairwise_differences)) +
#   geom_histogram() +
#   theme_bw() +
#   theme(axis.text=element_text(size=14), axis.title=element_text(size = 14))+
#   xlab("Pairwise Differences")
# ggsave("data/figures/S3B.pdf", height = 5, width = 5)
# ggsave("data/figures/S3B.jpg", height = 5, width = 5)
# 
# 
# ## Supplemental Figure 3C
# 
# GENOME.class@n.sites
# GENOME.class <- neutrality.stats(GENOME.class)
# get.neutrality(GENOME.class)[[1]]
# GENOME.class <- F_ST.stats(GENOME.class)
# GENOME.class <- diversity.stats(GENOME.class)
# get.diversity(GENOME.class)[[1]]
# get.diversity(GENOME.class)[[1]] / GENOME.class@n.sites
# 
# windowSize = 500
# slide <- sliding.window.transform(GENOME.class,windowSize, windowSize, type=2)
# slide <- diversity.stats(slide)
# nucdiv <- slide@nuc.diversity.within
# head(nucdiv)
# slide <- neutrality.stats(slide)
# tajD <- slide@Tajima.D
# popgen_stats <- data.frame(base_pair = 1:nrow(nucdiv)*windowSize/1000, nucleotide_diversity =  nucdiv[,1]/windowSize, tajima_D = tajD[,1])
# 
# p1 <- ggplot(popgen_stats, aes(x = base_pair, y = tajima_D)) + geom_point() +
#   theme_bw()  +
#   theme(axis.text=element_text(size=14), axis.title=element_text(size = 20))+
#   xlab("Position (kilobases)") + ylab("Tajima's D") +
#   ylim(c(-3, 3))
# p1
# ggsave("data/figures/S3C.pdf", width = 10)
# ggsave("data/figures/S3C.jpg", width = 10)
  
## Supplemental Figure 5 - Root-to-tip regressions
  
m1 <- lm(root_to_tip$distance ~ root_to_tip$date)
p1 <- ggplot(root_to_tip, aes(x = date, y = distance)) + 
  geom_point() + 
  theme_bw() +
  theme(axis.text = element_text(size = 16), axis.title=element_text(size = 16)) +
  geom_smooth(method="lm", fullrange=TRUE) + 
  xlim(2019.6, 2020.4)
p1
ggsave("data/figures/S5.pdf", height = 3, width = 8)
ggsave("data/figures/S5.jpg", height = 3, width = 8)

## Supplemental Figure 9 - Western MA Cases

## Western MA Cases

## Study Western MA cases
W_MA_ss <- sample_set[5:9,c("sample_id", "assembly_length_unambiguous", "assembly_fasta")]
W_MA_ss$exposure <- rep("Berkshire County Cluster", nrow(W_MA_ss))
W_MA_ss_partial <- W_MA_ss[!(W_MA_ss$sample_id %in% sample_set_unique_hq$sample_id),]
# W_MA_seqs <- pull_fastas_and_concatenate(W_MA_ss_partial, "data/wma_fasta", "data/sequences_and_trees/wma_concatenated_unaligned.fasta", do_not_download = FALSE, clean_dir = TRUE)
# merged_DNAss <- c(fasta_seqs, W_MA_seqs)
#merged_mafft_aln <- mafft(as.DNAbin(merged_DNAss), thread = 10)
#write.fas(merged_mafft_aln, "data/sequences_and_trees/merged_mafft_aligned.fasta") # write alignment
#merged_mafft_alnDNAss <-readDNAMultipleAlignment("data/sequences_and_trees/merged_mafft_aligned.fasta", format="fasta")
#DNAss <- as(substr(merged_mafft_alnDNAss, 268, ncol(merged_mafft_alnDNAss) - 230), "DNAStringSet")
system("iqtree -s data/sequences_and_trees/merged_mafft_aligned.fasta -bb 10000 -nt 10") # iqtree with 10000 UF bootstraps


#DPH_combined_tree <- read.newick("data/sequences_and_trees/DPH_combined.tree")
DPH_combined_tree <- read.iqtree("COVID/2020_07_06/alignments/combined.fasta.contree")
#DPH_combined_tree <- read.iqtree("data/sequences_and_trees/merged_mafft_aligned.fasta.contree")
DPH_combined_tree@phylo <- midpoint.root(DPH_combined_tree@phylo)

p1 <- ggtree(DPH_combined_tree) %<+% W_MA_ss + geom_tippoint(aes(color = exposure, subset=!is.na(exposure)), size = 3) + 
  #theme(legend.position = "none") + 
  geom_text2(aes(subset=(as.numeric(UFboot) > 80), label = UFboot), hjust = 1.5, vjust = -0.3) + 
  theme(legend.position = c(0.75, 0.2), legend.title = element_blank(), legend.text=element_text(size = 20)) + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size = 20)) + 
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values= cbPalette) #+ 
  #geom_tiplab()
p1
#ggsave("data/convenience_plots/S9_with_labels.pdf", height = 200, width = 20, limitsize = FALSE)
ggsave("data/figures/S9.pdf", height = 10, width = 8)
ggsave("data/figures/S9.jpg", height = 10, width = 8)

## Supplemental Figure 11 - iqtree with annotations by exposure

# revise exposure data
metadata_iqtree$new_exposure <- metadata_iqtree$exposure
metadata_iqtree$new_exposure[metadata_iqtree$exposure == "Western_MA"] <- NA
metadata_iqtree$new_exposure[metadata_iqtree$exposure == "Travel_associated"] <- NA
metadata_iqtree$new_exposure[metadata_iqtree$exposure == "Homeless"] <- "BHCHP"
metadata_iqtree$new_exposure[metadata_iqtree$exposure == "ConfA"] <- "Conference"
metadata_iqtree$new_exposure[metadata_iqtree$exposure == "SNFA"] <- "SNF"



p1 <- ggtree(iqtree) %<+% metadata_iqtree[,c("label", "new_exposure")] +
  theme_tree2() +
  geom_text2(aes(subset=(as.numeric(UFboot) > 80), label = UFboot), hjust = 1.5, vjust = -0.3) + 
  geom_tippoint( aes(color=new_exposure, subset = !is.na(new_exposure)), size = 2) +
  theme(legend.position = c(0.15, 0.75), legend.title = element_blank(), legend.text=element_text(size = 20)) +
  scale_color_manual(values=cbPalette) + 
  ylim(c(0, 800))
p1
ggsave("data_revision/figures2/S12.pdf", height = 10, width = 8)
ggsave("data_revision/figures2/S12.jpg", height = 10, width = 8)

## Supplemental Figure 13 - cluster investigation, requires addition of amp_seq sample and realignment

mp_seq_ss <- read_tsv("data/metadata/AmpSeq_Data.tsv") 
#fasta_seqs_w_as <- pull_fastas_and_concatenate(amp_seq_ss[amp_seq_ss[,1] == "MA_MGH_00905",], fasta_path = "data/amp_seq_fasta", mfa_path = "data/sequences_and_trees/concatenated_unaligned_w_amp_seq.fasta",clean_dir=TRUE, do_not_download = FALSE)
#fasta_seqs_w_as <- c(fasta_seqs, fasta_seqs_w_as)
#mafft_aln2 <- mafft(as.DNAbin(fasta_seqs_w_as), thread = 8, exec = "/opt/anaconda3/bin/mafft") # align and trim fastas
#write.fas(mafft_aln2, "data/sequences_and_trees/mafft_aligned_w_as.fasta") # write alignment
mafft_alnDNAss_w_as <-readDNAMultipleAlignment("data/sequences_and_trees/mafft_aligned_w_as.fasta", format="fasta")
DNAss_w_as <- as(substr(mafft_alnDNAss_w_as, 268, ncol(mafft_alnDNAss_w_as) - 230), "DNAStringSet")
DNAss_w_as <- relabel_DNAss(DNAss_w_as, sample_set_metadata) # need to fix this
#writeXStringSet(DNAss_w_as, "data/sequences_and_trees/trimmed_alignment_w_as.fasta", format = "fasta")

#system("fasttree -nt -gtr -quote data/sequences_and_trees/trimmed_alignment_w_as.fasta > data/sequences_and_trees/trimmed_alignment_w_as.tree") # fasttree
#system("trimal -in data/sequences_and_trees/trimmed_alignment_w_as.fasta -out data/sequences_and_trees/trimmed_alignment_w_as.phy -phylip")
#system("phyml --no_memory_check -i data/sequences_and_trees/trimmed_alignment_w_as.phy") # phyml tree
#system("rm data/sequences_and_trees/trimmed_alignment_w_as.phy")
#system("iqtree -s data/sequences_and_trees/trimmed_alignment_w_as.fasta -bb 10000 -nt 10") # iqtree with 10000 UF bootstraps
iqtree2 <- read.iqtree("data/sequences_and_trees/trimmed_alignment_w_as.fasta.contree")
iqtree2@phylo <- midpoint.root(iqtree2@phylo) # midpoint root ML tree

p1 <- ggtree(iqtree2) %<+% metadata_iqtree_clade_labels + 
  theme_tree2() + 
  geom_tippoint(aes(subset = !is.na(nosocomial_cluster), color = nosocomial_cluster), size = 6) + 
  geom_text2(aes(subset=!isTip & as.numeric(UFboot) > 80, label=round(UFboot, 2)), hjust = 1.5, vjust = -0.3)+
  theme(legend.position = c(0.15, 0.75), legend.title = element_blank(), legend.text=element_text(size = 20))+
  scale_color_manual(values=cbPalette) + 
  ylim(0,800)
p1
ggsave("data_revision/figures/S16.pdf", height = 12)
ggsave("data_revision/figures/S16.jpg", height = 12)

# Test for associations with homeless cohort
coinfection_counts <- c(12, 8, 302,1109 )
#coinfection_counts <- c(12, 8, 313,760 )
coinfections <- matrix(coinfection_counts,nrow = 2,  dimnames = list(Site = c("DPH", "MGH"), Coinfections = c("Yes", "No")))
fisher.test(coinfections)

# Supplemental Figure S7
p1 <- ggtree(nextstrain_tree, size = 0.3) %<+% nextstrain_metadata +
  theme(legend.position = c(0.1, 0.87), legend.title = element_blank(), legend.text=element_text(size = 12)) + 
  geom_tippoint(aes(color=MA, subset = !is.na(MA)), size = 1) +
  geom_tiplab()+
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  scale_color_manual(values=cbPalette)
p1

p2 <- ggtree(nextstrain_tree) %<+% nextstrain_metadata +
  geom_tippoint(aes(color = region, subset=!is.na(region)), size = 4) +
  geom_tiplab() + 
  scale_color_manual(values=cbPalette)
p2

ggsave("data/convenience_plots/S7.pdf", height = 1000, width =20, limitsize = FALSE)
#ggsave("data/convenience_plots/S7.jpg", height = 100, width = 8)


# Supplemental Figure 8
p2 <- ggtree(nextstrain_tree) %<+% nextstrain_metadata +
  geom_tippoint(aes(color = region, subset=!is.na(region)), size = 4) +
  geom_tiplab() + 
  scale_color_manual(values=cbPalette)
p2
p4 <- viewClade(node = 4825)
p4

ggsave("data/figures/S8.pdf", height = 10, width = 8)
ggsave("data/figures/S8.jpg", height = 10, width = 8)




# plot DPH cases by day

dph_march <- read_csv("data/metadata/DPH_march_reports.csv")
dph_march_longer <- gather(dph_march, Total, Conference_cluster, Travel, Under_Investigation, Berkshire_County_cluster, key = "Type", value = "Cases")
p1 <- ggplot(dph_march_longer, aes(x = date, y = Cases, color= Type)) + geom_point(size = 2.5) + geom_line(size = 1.5) + 
  theme_bw() + 
  theme(legend.position = c(0.4, 0.6), legend.text = element_text(size = 16), legend.title = element_blank(), legend.background=element_blank()) + 
  theme(axis.text=element_text(size = 14), axis.title=element_text(size = 14)) + xlab("") + 
  scale_color_manual(values=cbPalette)
p1
ggsave("data/figures/DPH_reports.pdf", height = 5, width = 5)
ggsave("data/figures/DPH_reports.jpg", height = 5, width = 5)

# plot beast tree with zip codes, Supplemental Figure 14

sample_set_metadata_beast$chelsea <- ifelse(sample_set_metadata_beast$zip_code_of_residence == "02150", "Chelsea", NA)
sample_set_metadata_beast$revere <- ifelse(sample_set_metadata_beast$zip_code_of_residence == "02151", "Revere", NA)
sample_set_metadata_beast$everett <- ifelse(sample_set_metadata_beast$zip_code_of_residence == "02149", "Everett", NA)
sample_set_metadata_beast$zip_bin <- rep(NA, nrow(sample_set_metadata_beast))
sample_set_metadata_beast$zip_bin[sample_set_metadata_beast$zip_code_of_residence == "02150"] <- "Chelsea"
sample_set_metadata_beast$zip_bin[sample_set_metadata_beast$zip_code_of_residence == "02151"] <- "Revere"
sample_set_metadata_beast$zip_bin[sample_set_metadata_beast$zip_code_of_residence == "02149"] <- "Everett"
p1 <- ggtree(beast, mrsd="2020-05-09") %<+% sample_set_metadata_beast + 
  theme_tree2() + #+ 
  geom_tippoint(aes(color=zip_bin, subset=!is.na(zip_bin)), size = 2) +
  geom_text2(aes(subset=!isTip & as.numeric(posterior) > 0.8, label=round(posterior, 2)), hjust = 1.5, vjust = -0.3) + 
  theme(legend.position = c(0.2, 0.85), legend.title = element_blank(), legend.text=element_text(size = 18)) + 
  #  theme(legend.position="none") + 
  scale_color_manual(values=cbPalette[c(1,3,5,2,4)]) + 
  theme(axis.text=element_text(size=15), axis.title=element_text(size = 20))

p1
ggsave("data/figures/S14.pdf", height = 10, width = 8)
ggsave("data/figures/S14.jpg", height = 10, width = 8)


# Supplemental Figure 15

# plot diffusion of C2416T

zip_zoom <- c(25017,25025, 25021, 25009, 25023, 25019, 25001, 25005)

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata[sample_set_metadata$exposure != "SNFA",], by = "sample_id")


zip_df <- data.frame(table(variant_association$zip_code_of_residence[variant_association$variant == "T"]))
names(zip_df) <- c("region", "value1")
zip_df$region <- as.character(zip_df$region)
zip_df$value1 <- as.numeric(zip_df$value1)

# plot map with zip code info
data(df_pop_zip)
ma_pop_zip <- df_pop_zip
ma_pop_zip <- left_join(ma_pop_zip, zip_df, by = "region")
ma_pop_zip2 <- df_pop_zip
ma_pop_zip2$value <- 0
ma_pop_zip2 <- left_join(ma_pop_zip2, zip_df, by = "region")
ma_pop_zip2$value1[is.na(ma_pop_zip2$value1)] <- 0
ma_pop_zip2 <- ma_pop_zip2[,c("region", "value1")]
ma_pop_ziplog <- data.frame(region = ma_pop_zip2$region, value = log10(ma_pop_zip2$value + 1))
names(ma_pop_zip2)[2] <- "value"

p1 <- zip_choropleth(ma_pop_ziplog, county_zoom = zip_zoom, state_zoom = "massachusetts", num_colors = 1, title = "       C2416T")
p1

# plot diffusion of G26233T

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata[sample_set_metadata$exposure != "SNFA",], by = "sample_id")


zip_df <- data.frame(table(variant_association$zip_code_of_residence[variant_association$variant == "T"]))
names(zip_df) <- c("region", "value1")
zip_df$region <- as.character(zip_df$region)
zip_df$value1 <- as.numeric(zip_df$value1)

# plot map with zip code info
data(df_pop_zip)
ma_pop_zip <- df_pop_zip
ma_pop_zip <- left_join(ma_pop_zip, zip_df, by = "region")
ma_pop_zip2 <- df_pop_zip
ma_pop_zip2$value <- 0
ma_pop_zip2 <- left_join(ma_pop_zip2, zip_df, by = "region")
ma_pop_zip2$value1[is.na(ma_pop_zip2$value1)] <- 0
ma_pop_zip2 <- ma_pop_zip2[,c("region", "value1")]
ma_pop_ziplog <- data.frame(region = ma_pop_zip2$region, value = log10(ma_pop_zip2$value + 1))
names(ma_pop_zip2)[2] <- "value"

p2 <- zip_choropleth(ma_pop_ziplog, county_zoom = zip_zoom, state_zoom = "massachusetts", num_colors = 1, title = "       G26233T")
p2

# plot diffusion of 1059T

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1059], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata[sample_set_metadata$exposure != "SNFA",], by = "sample_id")


zip_df <- data.frame(table(variant_association$zip_code_of_residence[variant_association$variant == "T"]))
names(zip_df) <- c("region", "value1")
zip_df$region <- as.character(zip_df$region)
zip_df$value1 <- as.numeric(zip_df$value1)

# plot map with zip code info
data(df_pop_zip)
ma_pop_zip <- df_pop_zip
ma_pop_zip <- left_join(ma_pop_zip, zip_df, by = "region")
ma_pop_zip2 <- df_pop_zip
ma_pop_zip2$value <- 0
ma_pop_zip2 <- left_join(ma_pop_zip2, zip_df, by = "region")
ma_pop_zip2$value1[is.na(ma_pop_zip2$value1)] <- 0
ma_pop_zip2 <- ma_pop_zip2[,c("region", "value1")]
ma_pop_ziplog <- data.frame(region = ma_pop_zip2$region, value = log10(ma_pop_zip2$value + 1))
names(ma_pop_zip2)[2] <- "value"

p3 <- zip_choropleth(ma_pop_ziplog, county_zoom = zip_zoom, state_zoom = "massachusetts", num_colors = 1, title = "       C1059T")
p3

# plot diffusion of 28603

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,28899], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata[sample_set_metadata$exposure != "SNFA",], by = "sample_id")


zip_df <- data.frame(table(variant_association$zip_code_of_residence[variant_association$variant == "T"]))
names(zip_df) <- c("region", "value1")
zip_df$region <- as.character(zip_df$region)
zip_df$value1 <- as.numeric(zip_df$value1)

# plot map with zip code info
data(df_pop_zip)
ma_pop_zip <- df_pop_zip
ma_pop_zip <- left_join(ma_pop_zip, zip_df, by = "region")
ma_pop_zip2 <- df_pop_zip
ma_pop_zip2$value <- 0
ma_pop_zip2 <- left_join(ma_pop_zip2, zip_df, by = "region")
ma_pop_zip2$value1[is.na(ma_pop_zip2$value1)] <- 0
ma_pop_zip2 <- ma_pop_zip2[,c("region", "value1")]
ma_pop_ziplog <- data.frame(region = ma_pop_zip2$region, value = log10(ma_pop_zip2$value + 1))
names(ma_pop_zip2)[2] <- "value"

p4 <- zip_choropleth(ma_pop_ziplog, county_zoom = zip_zoom, state_zoom = "massachusetts", num_colors = 1, title = "       G28899T")
p4


plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), label_size = 16)
ggsave("data/figures/S15.pdf", height = 5, width = 5)
ggsave("data/figures/S15.jpg", height = 5, width = 5)


### Look at allele frequency by county

## 2416T
# below is among all samples
county_props <- rep(0, 4)
# Suffolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Suffolk" & !is.na(sample_set_metadata$county),], by = "sample_id")

var_by_suffolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[1] <- sum(var_by_suffolk$count_derived)/(sum(var_by_suffolk$count_ancestral) + sum(var_by_suffolk$count_derived))

# Middlesex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Middlesex" & !is.na(sample_set_metadata$county),], by = "sample_id")

var_by_middlesex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[2] <- sum(var_by_middlesex$count_derived)/(sum(var_by_middlesex$count_ancestral) + sum(var_by_middlesex$count_derived))

# Essex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Essex" & !is.na(sample_set_metadata$county),], by = "sample_id")

var_by_essex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[3] <-sum(var_by_essex$count_derived)/(sum(var_by_essex$count_ancestral) + sum(var_by_essex$count_derived))

# Norfolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Norfolk" & !is.na(sample_set_metadata$county),], by = "sample_id")

var_by_norfolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[4] <-sum(var_by_norfolk$count_derived)/(sum(var_by_norfolk$count_ancestral) + sum(var_by_norfolk$count_derived))

county_counts_may_9 <- c(15423, 17694, 11429, 7172)
county_counts_aug_1 <- c(17305, 25801, 17305, 10305)

sum(county_props*county_counts_may_9)
sum(county_props*county_counts_aug_1)

## 2416T
# below is among samples excluding conference and SNF
county_props <- rep(0, 4)
# Suffolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Suffolk" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")
table(variant_association$county, variant_association$variant)
var_by_suffolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
sum(var_by_suffolk$count_derived)
sum(var_by_suffolk$count_ancestral)
county_props[1] <- sum(var_by_suffolk$count_derived)/(sum(var_by_suffolk$count_ancestral) + sum(var_by_suffolk$count_derived))

# Middlesex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Middlesex" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")
table(variant_association$county, variant_association$variant)
var_by_middlesex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[2] <- sum(var_by_middlesex$count_derived)/(sum(var_by_middlesex$count_ancestral) + sum(var_by_middlesex$count_derived))

sum(var_by_middlesex$count_derived)
sum(var_by_middlesex$count_ancestral)

# Essex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Essex" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")
table(variant_association$county, variant_association$variant)
var_by_essex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[3] <-sum(var_by_essex$count_derived)/(sum(var_by_essex$count_ancestral) + sum(var_by_essex$count_derived))

sum(var_by_essex$count_derived)
sum(var_by_essex$count_ancestral)

# Norfolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Norfolk" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")
table(variant_association$county, variant_association$variant)
var_by_norfolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[4] <-sum(var_by_norfolk$count_derived)/(sum(var_by_norfolk$count_ancestral) + sum(var_by_norfolk$count_derived))

sum(var_by_norfolk$count_derived)
sum(var_by_norfolk$count_ancestral)

county_counts_may_9 <- c(15423, 17694, 11429, 7172)
county_counts_aug_1 <- c(17305, 25801, 17305, 10305)

sum(county_props*county_counts_may_9)
sum(county_props*county_counts_aug_1)

# below is with "NoKnown" exposure

county_props <- rep(0, 4)
# Suffolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Suffolk" & sample_set_metadata$exposure == "NoKnown",], by = "sample_id")

var_by_suffolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
sum(var_by_suffolk$count_derived)
sum(var_by_suffolk$count_ancestral)
county_props[1] <- sum(var_by_suffolk$count_derived)/(sum(var_by_suffolk$count_ancestral) + sum(var_by_suffolk$count_derived))

# Middlesex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Middlesex" & sample_set_metadata$exposure == "NoKnown",], by = "sample_id")

var_by_middlesex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[2] <- sum(var_by_middlesex$count_derived)/(sum(var_by_middlesex$count_ancestral) + sum(var_by_middlesex$count_derived))

sum(var_by_middlesex$count_derived)
sum(var_by_middlesex$count_ancestral)

# Essex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Essex" & sample_set_metadata$exposure == "NoKnown",], by = "sample_id")

var_by_essex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[3] <-sum(var_by_essex$count_derived)/(sum(var_by_essex$count_ancestral) + sum(var_by_essex$count_derived))

sum(var_by_essex$count_derived)
sum(var_by_essex$count_ancestral)

# Norfolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Norfolk" & sample_set_metadata$exposure == "NoKnown",], by = "sample_id")

var_by_norfolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[4] <-sum(var_by_norfolk$count_derived)/(sum(var_by_norfolk$count_ancestral) + sum(var_by_norfolk$count_derived))

sum(var_by_norfolk$count_derived)
sum(var_by_norfolk$count_ancestral)

county_counts_may_9 <- c(15423, 17694, 11429, 7172)
county_counts_aug_1 <- c(17305, 25801, 17305, 10305)

sum(county_props*county_counts_may_9)
sum(county_props*county_counts_aug_1)


## 26233

county_props <- rep(0, 4)

# Suffolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Suffolk" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")

var_by_suffolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[1] <- sum(var_by_suffolk$count_derived)/(sum(var_by_suffolk$count_ancestral) + sum(var_by_suffolk$count_derived))

# Middlesex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Middlesex" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")

var_by_middlesex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[2] <- sum(var_by_middlesex$count_derived)/(sum(var_by_middlesex$count_ancestral) + sum(var_by_middlesex$count_derived))

# Essex
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Essex" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")

var_by_essex <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[3] <-sum(var_by_essex$count_derived)/(sum(var_by_essex$count_ancestral) + sum(var_by_essex$count_derived))

# Norfolk
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- inner_join(variant_association, sample_set_metadata[sample_set_metadata$county == "Norfolk" & !is.na(sample_set_metadata$county) & sample_set_metadata$exposure != "ConfA" & sample_set_metadata$exposure != "SNFA",], by = "sample_id")

var_by_norfolk <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"), prop_derived = sum(variant == "T")/ (sum(variant == "T") + sum(variant == "C")))
county_props[4] <-sum(var_by_norfolk$count_derived)/(sum(var_by_norfolk$count_ancestral) + sum(var_by_norfolk$count_derived))

county_counts_may_9 <- c(15423, 17694, 11429, 7172)
county_counts_aug_1 <- c(17305, 25801, 17305, 10305)

sum(county_props*county_counts_may_9)
sum(county_props*county_counts_aug_1)


# Supplemental Figure S16: Allele frequency by site, including RIC samples

ggtree(beast, mrsd="2020-05-09") %<+% sample_set_metadata_beast +
  theme_tree2() +
  geom_tippoint(aes(color=clinic, subset=!is.na(clinic)), size = 2) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  scale_fill_manual(cbPalette)


variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
variant_association <- variant_association[variant_association$clinic == "Chelsea RIC" & !is.na(variant_association$clinic),]

var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"))
#var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
#var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$count_derived/(var_by_date$count_ancestral + var_by_date$count_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$count_derived[i], var_by_date$count_ancestral[i] + var_by_date$count_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u", "n", "count_derived", "count_ancestral")]
var_tbl$variant <- rep("C2416T", nrow(var_tbl))
var_tbl2 <- var_tbl
p7 <- ggplot(data = subset(var_tbl2, var_tbl2$n > 2), aes(x = date, y = af, col = variant)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  geom_hline(yintercept = sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)), color = cbPalette[2])+
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank(), plot.title=element_text(size = 20, hjust = 0.5))+
  scale_color_manual(values=cbPalette[2]) + 
  ggtitle("Chelsea RIC")
p7
print(paste("Chelsea RIC C2416T:", sum(var_tbl2$count_derived)))
print(paste("Chelsea RIC total:", sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)))
print(paste("Chelsea RIC proportion C2416T:", round(sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)) ,3)))

# plot for 26233
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
variant_association <- variant_association[variant_association$clinic == "Chelsea RIC" & !is.na(variant_association$clinic),]

var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
#var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
#var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$count_derived/(var_by_date$count_ancestral + var_by_date$count_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$count_derived[i], var_by_date$count_ancestral[i] + var_by_date$count_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u", "n", "count_derived", "count_ancestral")]
var_tbl$variant <- rep("G26233T", nrow(var_tbl))
var_tbl2 <- var_tbl
p8 <- ggplot(data = subset(var_tbl2, var_tbl2$n > 2), aes(x = date, y = af, col = variant)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  geom_hline(yintercept = sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)), color = cbPalette[3])+
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank(), plot.title=element_text(size = 20, hjust = 0.5))+
  scale_color_manual(values=cbPalette[3])+
  ggtitle("Chelsea RIC")
p8
print(paste("Chelsea RIC G26233T:", sum(var_tbl2$count_derived)))
print(paste("Chelsea RIC total:", sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)))
print(paste("Chelsea RIC proportion G26233T:", round(sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)) ,3)))
 
# look in the homeless population

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
variant_association <- variant_association[variant_association$exposure == "Homeless",]

var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"))
#var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
#var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$count_derived/(var_by_date$count_ancestral + var_by_date$count_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$count_derived[i], var_by_date$count_ancestral[i] + var_by_date$count_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u", "n", "count_derived", "count_ancestral")]
var_tbl$variant <- rep("C2416T", nrow(var_tbl))
var_tbl2 <- var_tbl
p9 <- ggplot(data = subset(var_tbl2, var_tbl2$n > 2), aes(x = date, y = af, col = variant)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  geom_hline(yintercept = sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)), color = cbPalette[2])+
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank(), plot.title=element_text(size = 20, hjust = 0.5))+
  scale_color_manual(values=cbPalette[2]) +
  ggtitle("Homeless Cohort")
p9
print(paste("Homeless C2416T:", sum(var_tbl2$count_derived)))
print(paste("Homeless total:", sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)))
print(paste("Homeless proportion C2416T:", round(sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)) ,3)))

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
variant_association <- variant_association[variant_association$exposure == "Homeless",]

var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
#var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
#var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$count_derived/(var_by_date$count_ancestral + var_by_date$count_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$count_derived[i], var_by_date$count_ancestral[i] + var_by_date$count_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u", "n", "count_derived", "count_ancestral")]
var_tbl$variant <- rep("G26233T", nrow(var_tbl))
var_tbl2 <- var_tbl
p10 <- ggplot(data = subset(var_tbl2, var_tbl2$n > 2), aes(x = date, y = af, col = variant)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  geom_hline(yintercept = sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)), color = cbPalette[3])+
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank(), plot.title=element_text(size = 20, hjust = 0.5))+
  scale_color_manual(values=cbPalette[3]) + 
  ggtitle("Homeless Cohort")
p10
print(paste("Homeless G26233T:", sum(var_tbl2$count_derived)))
print(paste("Homeless total:", sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)))
print(paste("Homeless proportion G26233T:", round(sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)) ,3)))

# in the non-SNF, non-homeless, non-Chelsea RIC, non-Conf population

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
variant_association <- variant_association[variant_association$exposure == "NoKnown" & is.na(variant_association$clinic),]

var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"))
#var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
#var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$count_derived/(var_by_date$count_ancestral + var_by_date$count_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$count_derived[i], var_by_date$count_ancestral[i] + var_by_date$count_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u", "n", "count_derived", "count_ancestral")]
var_tbl$variant <- rep("C2416T", nrow(var_tbl))
var_tbl2 <- var_tbl
p11 <- ggplot(data = subset(var_tbl2, var_tbl2$n > 2), aes(x = date, y = af, col = variant)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  geom_hline(yintercept = sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)), color = cbPalette[2])+
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank(), plot.title=element_text(size = 20, hjust = 0.5))+
  scale_color_manual(values=cbPalette[2]) + 
  ggtitle("MGH (not in a cluster)")
p11
print(paste("MGH Background C2416T:", sum(var_tbl2$count_derived)))
print(paste("MGH Background total:", sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)))
print(paste("MGH Background proportion C2416T:", round(sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)) ,3)))

variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
variant_association <- variant_association[variant_association$exposure == "NoKnown" & is.na(variant_association$clinic),]

var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
#var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
#var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$count_derived/(var_by_date$count_ancestral + var_by_date$count_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$count_derived[i], var_by_date$count_ancestral[i] + var_by_date$count_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u", "n", "count_derived", "count_ancestral")]
var_tbl$variant <- rep("G26233T", nrow(var_tbl))
var_tbl2 <- var_tbl
p12 <- ggplot(data = subset(var_tbl2, var_tbl2$n > 2), aes(x = date, y = af, col = variant)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  geom_hline(yintercept = sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)), color = cbPalette[3])+
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 12), legend.title = element_blank(), plot.title=element_text(size = 20, hjust = 0.5))+
  scale_color_manual(values=cbPalette[3]) +
  ggtitle("MGH (not in a cluster)")
p12
print(paste("MGH Background G26233T:", sum(var_tbl2$count_derived)))
print(paste("MGH Background total:", sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)))
print(paste("MGH Background proportion G26233T:", round(sum(var_tbl2$count_derived) / (sum(var_tbl2$count_ancestral) + sum(var_tbl2$count_derived)) ,3)))

plot_grid(p7, p8, p9, p10, p11, p12, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"), label_size = 12)
ggsave("data/figures/S16.pdf", height = 12, width = 12)
ggsave("data/figures/S16.jpg", height = 12, width = 12)


# Figure S17
p6 <- ggplot(subset(dph_counts, dph_counts$Date < "2020-05-11"), aes(x = Date, y = Statewide)) + geom_point() + geom_smooth(color = "gray", span = 0.7) + theme_bw() + theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.65), legend.text=element_text(size = 14), legend.title = element_blank()) + 
  ylab("New Cases (Statewide)") + 
  theme(axis.title = element_text(size = 12)) +
  scale_color_manual(values=cbPalette[1])

p6  


variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,2416], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "C"), count_derived = sum(variant == "T"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u")]
var_tbl$variant <- rep("C2416T", nrow(var_tbl))
var_tbl2 <- var_tbl


variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,26233], sample_id = rownames(mafft_alnDNAss)))
variant_association <- left_join(variant_association, sample_set_metadata, by = "sample_id")
var_by_date <- variant_association %>% group_by(date) %>% summarise(n = n(), count_ancestral = sum(variant == "G"), count_derived = sum(variant == "T"))
var_by_date$cumsum_ancestral <- cumsum(var_by_date$count_ancestral)
var_by_date$cumsum_derived <- cumsum(var_by_date$count_derived)
var_by_date$af <- var_by_date$cumsum_derived/(var_by_date$cumsum_ancestral + var_by_date$cumsum_derived)
var_by_date$ci_l <- rep(NA, nrow(var_by_date))
var_by_date$ci_u <- rep(NA, nrow(var_by_date))
for(i in 1:nrow(var_by_date)){
  bin_test <- binom.test(var_by_date$cumsum_derived[i], var_by_date$cumsum_ancestral[i] + var_by_date$cumsum_derived[i])
  var_by_date$ci_l[i] <- bin_test$conf.int[1]
  var_by_date$ci_u[i] <- bin_test$conf.int[2]
}

var_tbl <- var_by_date[,c("date", "af", "ci_l", "ci_u")]
var_tbl$variant <- rep("G26233T", nrow(var_tbl))
var_tbl2 <- rbind(var_tbl, var_tbl2)

p7 <- ggplot(var_tbl2, aes(x = date, y = af, col = variant)) + geom_point(size = 2) + geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, col = variant), alpha = 0.3) + 
  geom_vline(xintercept = as_date("2020-02-28")) + theme_bw() + 
  xlab("Date") + ylab("Allele Frequency") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.75, 0.75), legend.text=element_text(size = 14), legend.title = element_blank())+
  scale_color_manual(values=cbPalette[c(2:1)])

dph_march <- read_csv("data/metadata/DPH_march_reports.csv")
dph_march_longer <- gather(dph_march, Total, Conference_cluster, Travel, Under_Investigation, key = "Type", value = "Cases")
p8 <- ggplot(dph_march_longer, aes(x = date, y = Cases, color= Type)) + geom_point(size = 2.5) + geom_line(size = 1.5) + 
  theme_bw() + 
  theme(legend.position = c(0.8, 0.5), legend.text = element_text(size = 14), legend.title = element_blank(), legend.background=element_blank()) + 
  theme(axis.text=element_text(size = 12), axis.title=element_text(size = 12)) + xlab("") + 
  scale_color_manual(values=cbPalette[c(5,2,8,6,3)]) +
  coord_cartesian(xlim =c(as_date("2020-02-28"), as_date("2020-05-13"))) +
  ylab("Cumulative Cases (Statewide)")

plot_grid(p8, p7, p6, nrow = 3, labels = c("A", "B", "C"), label_size = 16)
ggsave("data_revision/figures2/S13.pdf", height = 10, width = 8)
ggsave("data_revision/figures2/S13.jpg", height = 10, width = 8)

# Compare different substition models in BEAST

#beast_model2 <- read.beast("data/sequences_and_trees/MCC.tree") # MCC tree from BEAST
#mcmc_trace2 <- readLog("COVID/2020_07_13/coal_exp_prior_HKY_G4/trimmed_alignment_with_clades.log", burnin = 0.3) # BEAST log file
mcmc_trace2 <- readLog("COVID/2020_08_11/BEAST/HKY_G4_coal_exp_prior_2/trimmed_alignment_with_clades.log", burnin = 0.3)
mcmc_df2 <- data.frame(mcmc_trace2)
mcmc_df2$iteration <- rownames(mcmc_df2)

# clean up data to plot
MRCA_df2 <- mcmc_df2[,c("mrca.date.Full_Set.","mrca.date.Clade_1.","mrca.date.Clade_2.","mrca.date.Clade_3.","mrca.date.Clade_4.","mrca.date.Clade_5.","iteration")]
names(MRCA_df2) <- c("Root","C20099T (BHCHP)", "G3892T (SNF)" , "C2416T (Conf, BHCHP)","G105T (BHCHP)", "G28899T", "iteration")
MRCA_reshaped2 <- melt(MRCA_df2, id.vars = "iteration")

p2_2 <- ggplot(MRCA_reshaped, aes(y = value, x= variable, fill = variable)) +
  geom_violin(adjust = 2, draw_quantiles = c(0.5)) +
  theme_classic() +
  #  theme(axis.text=element_text(size=15, angle = 90), axis.title=element_text(size = 20)) + xlab("") + ylab("Date") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95, angle = 90, size = 14), axis.title=element_text(size = 20), plot.title=element_text(hjust = 0.5)) + 
  xlab("") + ylab("Date") +
  theme(legend.position = "none") + 
  scale_fill_manual(values=cbPalette[c(6,1:5)]) + 
  ggtitle("GTR4G")
p2_2

p3_2 <- ggplot(MRCA_reshaped2, aes(y = value, x= variable, fill = variable)) +
  geom_violin(adjust = 2, draw_quantiles = c(0.5)) +
  theme_classic() +
  #  theme(axis.text=element_text(size=15, angle = 90), axis.title=element_text(size = 20)) + xlab("") + ylab("Date") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95, angle = 90, size = 14), axis.title=element_text(size = 20), plot.title=element_text(hjust = 0.5)) + 
  xlab("") + ylab("Date") +
  theme(legend.position = "none") + 
  scale_fill_manual(values=cbPalette[c(6,1:5)]) + 
  ggtitle("HKY4G") 
p3_2

p4_1 <- ggplot(data.frame(mcmc_trace), aes( x = clockRate)) + geom_density(fill = "grey") + theme_bw() + 
  ggtitle("GTR4G") + theme(plot.title=element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(0.0008, 0.001, 0.0012))
p4_1
p4_2 <- ggplot(data.frame(mcmc_trace2), aes( x = clockRate)) + geom_density(fill = "grey") + theme_bw() + 
  ggtitle("HKY4G") + 
  theme(plot.title=element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(0.0008, 0.001, 0.0012))
p4_2


p5_1 <- ggplot(data.frame(mcmc_trace), aes( x = ePopSize)) + geom_density(fill = "grey") + theme_bw() + 
  ggtitle("GTR4G") + theme(plot.title=element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(5, 7.5, 10, 12.5))
p5_2 <- ggplot(data.frame(mcmc_trace2), aes( x = ePopSize)) + geom_density(fill = "grey") + theme_bw() + 
  ggtitle("HKY4G") + 
  theme(plot.title=element_text(hjust = 0.5))  + 
  scale_x_continuous(breaks = c(5, 7.5, 10, 12.5))
p6_1 <- ggplot(data.frame(mcmc_trace), aes( x = growthRate)) + geom_density(fill = "grey") + theme_bw() + 
  ggtitle("GTR4G") + theme(plot.title=element_text(hjust = 0.5))  + 
  scale_x_continuous(breaks = c(10, 15))
p6_2 <- ggplot(data.frame(mcmc_trace2), aes( x = growthRate)) + geom_density(fill = "grey") + theme_bw() + 
  ggtitle("HKY4G") + 
  theme(plot.title=element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = c(10, 15))

plot_grid(p4_1, p4_2, p5_1, p5_2, p6_1, p6_2,p2_2, p3_2, labels = c("B", "C", "D", "E", "F", "G", "H", "I"), nrow = 4 , 
          label_size = 20, rel_heights = c(1,1, 1, 3))


ggsave("data_revision/figures/S5_B-I.pdf", height = 8, width = 6)
ggsave("data_revision/figures/S5_B-I.jpg", height = 8, width = 6)


