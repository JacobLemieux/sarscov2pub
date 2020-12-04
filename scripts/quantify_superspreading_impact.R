### Quantify impact of superspreading state-by-state
# Jacob E. Lemieux
# November 09 2020

library(Biostrings)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(cowplot)
library(dplyr)
library(MCMCglmm)
library(choroplethr)
library(choroplethrMaps)
data(state.map)

setwd("COVID/SARS-CoV-2/")
source("scripts/SARS-CoV-2_functions.R")

# colorblind palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


# initialize list of simulations for total
state_totals <- list(length(10))

## Read in alignment and trim UTRs

#### Look at USA Data
mafft_alnDNAss <-readDNAMultipleAlignment("COVID/2020_11_02/msa_1102/subset2416_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

# filter
variant_association <- variant_association %>% filter(variant_association$country == "USA")
variant_association$state <- substr(variant_association$seq_name, 1, 2)
variants_by_state <- as_tibble(as.data.frame.matrix(table(variant_association$state, variant_association$variant)[,c(2,4)]), rownames="state")
variants_by_state$AF <- variants_by_state$T/(variants_by_state$C + variants_by_state$T)
variants_by_state$ci_l <- rep(NA, nrow(variants_by_state))
variants_by_state$ci_u <- rep(NA, nrow(variants_by_state))
for(i in 1:nrow(variants_by_state)){
  bin_test <- binom.test(variants_by_state$T[i], variants_by_state$C[i] + variants_by_state$T[i])
  variants_by_state$ci_l[i] <- bin_test$conf.int[1]
  variants_by_state$ci_u[i] <- bin_test$conf.int[2]
}

# case counts without division by epoch
doi <- "2020-11-01"
case_counts <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == doi)
case_counts$state <- sapply(case_counts$state, state_to_abbrev)
variants_by_state <- left_join(variants_by_state, case_counts, by = "state")
variants_by_state$linked_cases <- variants_by_state$AF * variants_by_state$cases
variants_by_state <- variants_by_state %>% filter(!is.na(cases))
variants_by_state$alleles <- variants_by_state$C + variants_by_state$T


# create choropleth map
names(variants_by_state)[8] <- "STATE"
state_key <- state.map[,c("STATE", "region")]
state_key <- state_key[!duplicated(state_key$STATE),]
s_df <- left_join(variants_by_state, state_key, by = "STATE")
s_df <- s_df[,c("region", "AF")]
names(s_df)[2] <- "value"
s_df <- rbind(s_df, c("north dakota", "0"))
s_df$value <- as.numeric(s_df$value)
s0 <- state_choropleth(s_df, num_colors = 1, legend = "Allele Frequency")
s0t <- state_choropleth(s_df, num_colors = 1, title = "      C2416T", legend = "Allele Frequency")
s0t
ggsave("data_revision/figures2/C2416T_with_abbreviations.pdf")
ggsave("data_revision/figures2/C2416T_with_abbreviations.jpg")

s1 <- StateChoropleth$new(s_df)
s1$show_labels = FALSE
s1$set_num_colors(1)
s1$render()
ggsave("data_revision/figures2/C2416T_no_abbreviations.pdf")
ggsave("data_revision/figures2/C2416T_no_abbreviations.jpg")

# filter for only in those places that have alleles
variants_by_state <- variants_by_state %>% filter(T > 0 & alleles > 10)

# simulate
C2416T_sims <- sim2416(10000, variants_by_state)

# plot
p1 <- ggplot(variants_by_state, aes(x = state, y = AF)) + geom_point() + 
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3) + 
  theme_bw() + 
  ylab("Allele Frequency") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) + 
  ggtitle("C2416T")

p2 <- ggplot(variants_by_state, aes(x = state, y = linked_cases)) + geom_point() + 
  geom_errorbar(aes(ymin = cases*ci_l, ymax = cases*ci_u), alpha = 0.3) + 
  theme_bw() + 
  ylab("Cases") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  ggtitle("C2416T")
#plot_grid(p1, p2, labels = c("A", "B"), nrow = 2)


# Look at 26233

mafft_alnDNAss <-readDNAMultipleAlignment("COVID/2020_11_02/msa_1102/subset26233_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

variant_association <- variant_association %>% filter(country == "USA")
variant_association$state <- substr(variant_association$seq_name, 1, 2)
variants_by_state <- as_tibble(as.data.frame.matrix(table(variant_association$state, variant_association$variant)[,c(2,4)]), rownames="state")
variants_by_state$AF <- variants_by_state$T/(variants_by_state$G + variants_by_state$T)
variants_by_state$ci_l <- rep(NA, nrow(variants_by_state))
variants_by_state$ci_u <- rep(NA, nrow(variants_by_state))
for(i in 1:nrow(variants_by_state)){
  bin_test <- binom.test(variants_by_state$T[i], variants_by_state$G[i] + variants_by_state$T[i])
  variants_by_state$ci_l[i] <- bin_test$conf.int[1]
  variants_by_state$ci_u[i] <- bin_test$conf.int[2]
}


# add in total counts of population

doi <- "2020-11-01"
case_counts <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == doi)
case_counts$state <- sapply(case_counts$state, state_to_abbrev)
variants_by_state <- left_join(variants_by_state, case_counts, by = "state")
variants_by_state$linked_cases <- variants_by_state$AF * variants_by_state$cases
variants_by_state <- variants_by_state %>% filter(!is.na(cases))
variants_by_state$alleles <- variants_by_state$G + variants_by_state$T

# create choropleth map

names(variants_by_state)[8] <- "STATE"
state_key <- state.map[,c("STATE", "region")]
state_key <- state_key[!duplicated(state_key$STATE),]
s_df <- left_join(variants_by_state, state_key, by = "STATE")
s_df <- s_df[,c("region", "AF")]
names(s_df)[2] <- "value"
s_df <- rbind(s_df, c("north dakota", "0"))
s_df$value <- as.numeric(s_df$value)
s2 <- state_choropleth(s_df, num_colors = 1, legend = "Allele Frequency")
s2t <- state_choropleth(s_df, num_colors = 1, legend = "Allele Frequency", title = "      G26233T")
ggsave("data_revision/figures2/G26233T_with_abbreviations.pdf")
ggsave("data_revision/figures2/G26233T_with_abbreviations.jpg")

s3 <- StateChoropleth$new(s_df)
s3$show_labels = FALSE
s3$set_num_colors(1)
s3$render()
ggsave("data_revision/figures2/G26233T_no_abbreviations.pdf")
ggsave("data_revision/figures2/G26233T_no_abbreviations.jpg")

# simulate only in those plates that have alleles
variants_by_state <- variants_by_state %>% filter(variants_by_state$T > 0 & variants_by_state$alleles > 10)

# simulate
C26233T_sims <- sim26233(10000, variants_by_state)

# plot

p4 <- ggplot(variants_by_state, aes(x = state, y = AF)) + geom_point() + 
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3) + 
  theme_bw() + 
  ylab("Allele Frequency") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) + 
  ggtitle("G26233T")

p5 <- ggplot(variants_by_state, aes(x = state, y = linked_cases)) + geom_point() + 
  geom_errorbar(aes(ymin = cases*ci_l, ymax = cases*ci_u), alpha = 0.3) + 
  theme_bw() + 
  ylab("Cases") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  ggtitle("G26233T")

# Look at global 26233

mafft_alnDNAss <-readDNAMultipleAlignment("COVID/2020_11_02/msa_1102/subset26233_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

variants_by_country <- as_tibble(as.data.frame.matrix(table(variant_association$country, variant_association$variant)[,c(2,5)]), rownames="country")
variants_by_country <- variants_by_country %>% filter(variants_by_country$T > 0)
variants_by_country$AF <- variants_by_country$T/(variants_by_country$G + variants_by_country$T)
variants_by_country$ci_l <- rep(NA, nrow(variants_by_state))
variants_by_country$ci_u <- rep(NA, nrow(variants_by_state))
for(i in 1:nrow(variants_by_country)){
  bin_test <- binom.test(variants_by_country$T[i], variants_by_country$G[i] + variants_by_country$T[i])
  variants_by_country$ci_l[i] <- bin_test$conf.int[1]
  variants_by_country$ci_u[i] <- bin_test$conf.int[2]
}

# add in total counts of population
doi <- "2020-11-01"
country_counts <- read_csv("data_revision/metadata/JHU_covid_tracker_global_data_11_09_2020.csv")
country_counts <- country_counts %>% group_by(`Country/Region`) %>% summarise(nov_1_total = sum(`11/1/20`))
names(country_counts) <- c("country", "cases")
country_counts$country <- gsub("US", "USA",country_counts$country)
country_counts$country <- gsub("Czechia", "Czech Republic",country_counts$country)
country_counts$country <- gsub("United Kingdom", "England",country_counts$country)

variants_by_country <- left_join(variants_by_country, country_counts, by = "country")
variants_by_country$linked_cases <- variants_by_country$AF * variants_by_country$cases
variants_by_country <- variants_by_country %>% filter(!is.na(cases))
variants_by_country <- variants_by_country %>% filter(country != "USA")
sum(variants_by_country$linked_cases, na.rm=T)

# simulate number of cases based on estimated allele frequency and potential bias
variants_by_country$alleles <- variants_by_country$G + variants_by_country$T

C26233T_global <- sim26233(10000, variants_by_country)

# plot
p7 <- ggplot(variants_by_country, aes(x = country, y = AF)) + geom_point() + 
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3) + 
  theme_bw() + 
  ylab("Allele Frequency of G26233T") +
  xlab("")+
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) 

p8 <- ggplot(variants_by_country, aes(x = country, y = linked_cases)) + geom_point() + 
  geom_errorbar(aes(ymin = cases*ci_l, ymax = cases*ci_u), alpha = 0.3) + 
  theme_bw() + 
  ylab("Case Counts of G26233T") +
  xlab("") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) 

# Time-dependent model 2416
mafft_alnDNAss <-readDNAMultipleAlignment("COVID/2020_11_02/msa_1102/subset2416_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
date_zero <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
for(dz in date_zero){
  variant_association$sample_id <- gsub(paste("2020-",dz,"-00", sep=""), paste("2020-",dz,"-15", sep=""), variant_association$sample_id)
}
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

# filter
epoch1 <- "2020-05-13"
variant_association$epoch <- ifelse(variant_association$date <= epoch1, 0,1)
variant_association <- variant_association %>% 
  filter(variant_association$country == "USA")

variant_association$state <- substr(variant_association$seq_name, 1, 2)
var_by_epoch <- table(variant_association$state, variant_association$variant, variant_association$epoch)
variants_by_state <- as_tibble(as.data.frame.matrix(var_by_epoch[,c(2,4),1]), rownames="state")
names(variants_by_state) <- c("state", "C0", "T0")
variants_by_state_epoch_2 <- as_tibble(as.data.frame.matrix(var_by_epoch[,c(2,4),2]), rownames="state")
names(variants_by_state_epoch_2) <- c("state", "C1", "T1")
variants_by_state <- left_join(variants_by_state, variants_by_state_epoch_2, by = "state")
variants_by_state$AF0 <- variants_by_state$T0/(variants_by_state$C0 + variants_by_state$T0)
variants_by_state$ci_l0 <- rep(NA, nrow(variants_by_state))
variants_by_state$ci_u0 <- rep(NA, nrow(variants_by_state))
variants_by_state$AF1 <- variants_by_state$T1/(variants_by_state$C1 + variants_by_state$T1)
variants_by_state$ci_l1 <- rep(NA, nrow(variants_by_state))
variants_by_state$ci_u1 <- rep(NA, nrow(variants_by_state))

for(i in 1:nrow(variants_by_state)){
  if((variants_by_state$C0[i] + variants_by_state$T0[i]) > 0){
    bin_test0 <- binom.test(variants_by_state$T0[i], variants_by_state$C0[i] + variants_by_state$T0[i])
    variants_by_state$ci_l0[i] <- bin_test0$conf.int[1]
    variants_by_state$ci_u0[i] <- bin_test0$conf.int[2]
  }
    if((variants_by_state$C1[i] + variants_by_state$T1[i]) > 0){
    bin_test1 <- binom.test(variants_by_state$T1[i], variants_by_state$C1[i] + variants_by_state$T1[i])
    variants_by_state$ci_l1[i] <- bin_test1$conf.int[1]
    variants_by_state$ci_u1[i] <- bin_test1$conf.int[2]
  }
}

doi <- "2020-11-01"
case_counts <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == doi)
case_counts$state <- sapply(case_counts$state, state_to_abbrev)
case_counts_epoch1 <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == epoch1)
case_counts_epoch1$state <- sapply(case_counts_epoch1$state, state_to_abbrev)
case_counts$cases_epoch1 <- case_counts_epoch1$cases
case_counts$cases_epoch2 <- case_counts$cases - case_counts_epoch1$cases

variants_by_state <- left_join(variants_by_state, case_counts, by = "state")
variants_by_state$linked_cases_epoch1 <- variants_by_state$AF0 * variants_by_state$cases_epoch1
variants_by_state$linked_cases_epoch2 <- variants_by_state$AF1 * variants_by_state$cases_epoch2
variants_by_state$linked_cases <- apply(cbind(variants_by_state$linked_cases_epoch1, variants_by_state$linked_cases_epoch2), 1, function(x) sum(x, na.rm=T))
variants_by_state <- variants_by_state %>% filter(!is.na(cases))
sum(variants_by_state$linked_cases, na.rm=T)

# count up alleles per epoch

variants_by_state$alleles_epoch1 <- variants_by_state$C0 + variants_by_state$T0
variants_by_state$alleles_epoch2 <- variants_by_state$C1 + variants_by_state$T1

# filter
variants_by_state$alleles_total <- variants_by_state$alleles_epoch1 + variants_by_state$alleles_epoch2
variants_by_state$alleles_T_total <- variants_by_state$T0 + variants_by_state$T1
variants_by_state <- variants_by_state %>% filter(alleles_total > 9 & alleles_T_total > 0)

# simulate cases by state using a time-dependent model
td1 <- simTimeDep2416(10000, variants_by_state)

# plot the raw data
vbs_e1 <- variants_by_state[,c(1,6:8, 16, 18)]
vbs_e1 <- cbind(vbs_e1, rep("Epoch 1", nrow(vbs_e1)))
names(vbs_e1) <- c("state","AF", "ci_l", "ci_u", "cases", "linked_cases", "Epoch")
vbs_e2 <- variants_by_state[,c(1,9:11,17, 19)]
vbs_e2 <- cbind(vbs_e2, rep("Epoch 2", nrow(vbs_e2)))
names(vbs_e2) <- c("state","AF", "ci_l", "ci_u", "cases", "linked_cases", "Epoch")
vbs_df <- rbind(vbs_e1, vbs_e2)

p11a <- ggplot(vbs_df, aes(x = state, y = AF, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Allele Frequency") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.8, 0.815), legend.title=element_blank(), legend.background = element_blank(), legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) + 
  scale_color_manual(labels = c("Through May 13 2020", "After May 13 2020"), values = cbPalette) +
  ggtitle("C2416T")
p11a

p12a <- ggplot(vbs_df, aes(x = state, y = linked_cases, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = cases*ci_l, ymax = cases*ci_u), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Cases") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  scale_color_manual(labels = c("Through May 13 2020", "After May 13 2020"), values = cbPalette) +
  ggtitle("C2416T") + 
  coord_cartesian(ylim = c(0, 280000)) 
p12a
plot_grid(p11a, p12a, nrow = 2)

names(variants_by_state)[1] <- "STATE"
state_key <- state.map[,c("STATE", "region")]
state_key <- state_key[!duplicated(state_key$STATE),]
s_df <- left_join(variants_by_state, state_key, by = "STATE")
s_df <- s_df[,c("region", "AF1")]
names(s_df)[2] <- "value"
s_df <- rbind(s_df, c("north dakota", "0"))
s_df$value <- as.numeric(s_df$value)
#s0 <- state_choropleth(s_df, num_colors = 1, legend = "Allele Frequency")
#s0t <- state_choropleth(s_df, num_colors = 1, title = "      C2416T", legend = "Allele Frequency")
#s0t


sim_by_state_epoch1 <- as.data.frame(t(apply(td1[[4]],2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)) )))
sim_by_state_epoch2 <- as.data.frame(t(apply(td1[[5]],2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)) )))
sim_by_state_epoch1$state <- variants_by_state$STATE
sim_by_state_epoch2$state <- variants_by_state$STATE
sim_by_state_epoch1$Epoch <- rep("Epoch1", nrow(sim_by_state_epoch1))
sim_by_state_epoch2$Epoch <- rep("Epoch2", nrow(sim_by_state_epoch2))
sim_by_state <- rbind(sim_by_state_epoch1, sim_by_state_epoch2)
colnames(sim_by_state) <- c("sim_l_ci", "sim_median", "sim_u_ci", "state", "Epoch")

p13a <- ggplot(sim_by_state, aes(x = state, y = sim_median, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = sim_l_ci, ymax = sim_u_ci), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Cases (model)") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  scale_color_manual(labels = c("Through May 13 2020", "After May 13 2020"), values = cbPalette) +
  ggtitle("C2416T") + 
  coord_cartesian(ylim = c(0, 280000))
p13a
plot_grid(p11a, p12a, p13a, nrow = 3)


# are allele frequencies different in Florida between the two periods?
row1 <- variants_by_state[8,c(2,3)]
names(row1) <- c("Success", "Total")

row2 <- variants_by_state[8,c(4,5)]
names(row2) <- c("Success", "Total")

fl_af <- rbind(row1, row2)
fisher.test(t(fl_af))


# Time dependent model for 26233

mafft_alnDNAss <-readDNAMultipleAlignment("COVID/2020_11_02/msa_1102/subset26233_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
date_zero <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
for(dz in date_zero){
  variant_association$sample_id <- gsub(paste("2020-",dz,"-00", sep=""), paste("2020-",dz,"-15", sep=""), variant_association$sample_id)
}
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

# filter

epoch1 <- "2020-05-13"
variant_association$epoch <- ifelse(variant_association$date <= epoch1, 0,1)
variant_association <- variant_association %>% 
  filter(variant_association$country == "USA")

variant_association$state <- substr(variant_association$seq_name, 1, 2)
var_by_epoch <- table(variant_association$state, variant_association$variant, variant_association$epoch)
variants_by_state <- as_tibble(as.data.frame.matrix(var_by_epoch[,c(2,4),1]), rownames="state")
names(variants_by_state) <- c("state", "G0", "T0")
variants_by_state_epoch_2 <- as_tibble(as.data.frame.matrix(var_by_epoch[,c(2,4),2]), rownames="state")
names(variants_by_state_epoch_2) <- c("state", "G1", "T1")
variants_by_state <- left_join(variants_by_state, variants_by_state_epoch_2, by = "state")


variants_by_state$AF0 <- variants_by_state$T0/(variants_by_state$G0 + variants_by_state$T0)
variants_by_state$ci_l0 <- rep(NA, nrow(variants_by_state))
variants_by_state$ci_u0 <- rep(NA, nrow(variants_by_state))
variants_by_state$AF1 <- variants_by_state$T1/(variants_by_state$G1 + variants_by_state$T1)
variants_by_state$ci_l1 <- rep(NA, nrow(variants_by_state))
variants_by_state$ci_u1 <- rep(NA, nrow(variants_by_state))

for(i in 1:nrow(variants_by_state)){
  if((variants_by_state$G0[i] + variants_by_state$T0[i]) > 0){
    bin_test0 <- binom.test(variants_by_state$T0[i], variants_by_state$G0[i] + variants_by_state$T0[i])
    variants_by_state$ci_l0[i] <- bin_test0$conf.int[1]
    variants_by_state$ci_u0[i] <- bin_test0$conf.int[2]
  }
  if((variants_by_state$G1[i] + variants_by_state$T1[i]) > 0){
    bin_test1 <- binom.test(variants_by_state$T1[i], variants_by_state$G1[i] + variants_by_state$T1[i])
    variants_by_state$ci_l1[i] <- bin_test1$conf.int[1]
    variants_by_state$ci_u1[i] <- bin_test1$conf.int[2]
  }
}

# merge in case counts

doi <- "2020-11-01"
case_counts <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == doi)
case_counts$state <- sapply(case_counts$state, state_to_abbrev)
case_counts_epoch1 <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == epoch1)
case_counts_epoch1$state <- sapply(case_counts_epoch1$state, state_to_abbrev)
case_counts$cases_epoch1 <- case_counts_epoch1$cases
case_counts$cases_epoch2 <- case_counts$cases - case_counts_epoch1$cases

variants_by_state <- left_join(variants_by_state, case_counts, by = "state")
variants_by_state$linked_cases_epoch1 <- variants_by_state$AF0 * variants_by_state$cases_epoch1
variants_by_state$linked_cases_epoch2 <- variants_by_state$AF1 * variants_by_state$cases_epoch2
variants_by_state$linked_cases <- apply(cbind(variants_by_state$linked_cases_epoch1, variants_by_state$linked_cases_epoch2), 1, function(x) sum(x, na.rm=T))
variants_by_state <- variants_by_state %>% filter(!is.na(cases))
sum(variants_by_state$linked_cases, na.rm=T)

# count up alleles per epoch

variants_by_state$alleles_epoch1 <- variants_by_state$G0 + variants_by_state$T0
variants_by_state$alleles_epoch2 <- variants_by_state$G1 + variants_by_state$T1

# filter
variants_by_state$alleles_total <- variants_by_state$alleles_epoch1 + variants_by_state$alleles_epoch2
variants_by_state$alleles_T_total <- variants_by_state$T0 + variants_by_state$T1
variants_by_state <- variants_by_state %>% filter(alleles_total > 9 & alleles_T_total > 0)

# simulate cases by state using a time-dependent model
td2 <- simTimeDep26233(10000, variants_by_state)

# plot the raw data
vbs_e1 <- variants_by_state[,c(1,6:8, 16, 18)]
vbs_e1 <- cbind(vbs_e1, rep("Epoch 1", nrow(vbs_e1)))
names(vbs_e1) <- c("state","AF", "ci_l", "ci_u", "cases", "linked_cases", "Epoch")
vbs_e2 <- variants_by_state[,c(1,9:11,17, 19)]
vbs_e2 <- cbind(vbs_e2, rep("Epoch 2", nrow(vbs_e2)))
names(vbs_e2) <- c("state","AF", "ci_l", "ci_u", "cases", "linked_cases", "Epoch")
vbs_df <- rbind(vbs_e1, vbs_e2)



p11b <- ggplot(vbs_df, aes(x = state, y = AF, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Allele Frequency") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) + 
  scale_color_manual(labels = c("Through May 13 2020", "After May 13 2020"), values = cbPalette) +
  ggtitle("C26233T")
p11b

p12b <- ggplot(vbs_df, aes(x = state, y = linked_cases, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = cases*ci_l, ymax = cases*ci_u), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Cases") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  scale_color_manual(labels = c("Through May 13 2020", "After May 13 2020"), values = cbPalette) +
  ggtitle("26233T") + 
  coord_cartesian(ylim = c(0, 160000))
p12b

names(variants_by_state)[1] <- "STATE"
state_key <- state.map[,c("STATE", "region")]
state_key <- state_key[!duplicated(state_key$STATE),]
s_df <- left_join(variants_by_state, state_key, by = "STATE")
s_df <- s_df[,c("region", "AF1")]
names(s_df)[2] <- "value"
s_df <- rbind(s_df, c("north dakota", "0"))
s_df$value <- as.numeric(s_df$value)
#s0 <- state_choropleth(s_df, num_colors = 1, legend = "Allele Frequency")
#s0t <- state_choropleth(s_df, num_colors = 1, title = "      C2416T", legend = "Allele Frequency")
#s0t


sim_by_state_epoch1 <- as.data.frame(t(apply(td2[[4]],2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)) )))
sim_by_state_epoch2 <- as.data.frame(t(apply(td2[[5]],2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)) )))
sim_by_state_epoch1$state <- variants_by_state$STATE
sim_by_state_epoch2$state <- variants_by_state$STATE
sim_by_state_epoch1$Epoch <- rep("Epoch1", nrow(sim_by_state_epoch1))
sim_by_state_epoch2$Epoch <- rep("Epoch2", nrow(sim_by_state_epoch2))
sim_by_state <- rbind(sim_by_state_epoch1, sim_by_state_epoch2)
colnames(sim_by_state) <- c("sim_l_ci", "sim_median", "sim_u_ci", "state", "Epoch")

p13b <- ggplot(sim_by_state, aes(x = state, y = sim_median, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = sim_l_ci, ymax = sim_u_ci), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Cases (model)") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  ggtitle("C26233T") +
  scale_color_manual(labels = c("Through May 13 2020", "After May 13 2020"), values = cbPalette) + 
  coord_cartesian(ylim = c(0, 160000))

p13b
plot_grid(p11b, p12b, p13b, nrow = 3)






# have objects C2416T_sim, C26233T_sim, C26233T_global_sim, td1, td2
# regular sims have 1 - np, 2 - robust1, 3 - robust 2, 4 - bb
# td sims have 1 - nb1, 2 - nb2, 3 - nb total, 4 - bb1, 5 - bb2, 6 - bb total

# make figures

# beta-binomial 2416T
C2416T_df <- data.frame(a1 = apply(C2416T_sims[[4]], 1, sum), a2 = apply(td1[[4]],1, sum), a3 = apply(td1[[5]],1,sum), a4 = apply(td1[[6]], 1, sum)) %>% 
  pivot_longer(c("a1", "a2", "a3", "a4"))
C2416T_df$Variant = "C2416T"
# beta-binomial C26233T
G26233T_df <- data.frame(b1 = apply(C26233T_sims[[4]], 1, sum), b2 = apply(td2[[4]],1, sum), b3 = apply(td2[[5]],1,sum), b4 = apply(td2[[6]], 1, sum), 
                         b5 = apply(C26233T_global[[4]],1,sum), b6 = (apply(C26233T_global[[4]], 1, sum) + apply(C26233T_sims[[4]], 1, sum))) %>%
  pivot_longer(c("b1", "b2", "b3", "b4", "b5", "b6"))
G26233T_df$Variant = "G26233T"

# normal approx 2416T
C2416T_df_nl <- data.frame(c1 = apply(C2416T_sims[[1]], 1, sum), c2 = apply(td1[[1]],1, sum), c3 = apply(td1[[2]],1,sum), c4 = apply(td1[[3]], 1, sum),
                           c5 = apply(C2416T_sims[[2]], 1, sum), c6 = apply(C2416T_sims[[3]], 1, sum)) %>% 
  pivot_longer(c("c1", "c2", "c3", "c4", "c5", "c6"))
C2416T_df_nl$Variant = "C2416T"
# normal approx 2416T
G26233T_df_nl <- data.frame(d1 = apply(C26233T_sims[[1]], 1, sum), d2 = apply(td2[[1]],1, sum), d3 = apply(td2[[2]],1,sum), d4 = apply(td2[[3]], 1, sum),
                            d5 = apply(C26233T_sims[[2]], 1, sum), d6 = apply(C26233T_sims[[3]], 1, sum),
                            d7 = apply(C26233T_global[[1]],1,sum), d8 = apply(C26233T_global[[2]],1,sum), d9 = apply(C26233T_global[[3]],1,sum), 
                            d91 = (apply(C26233T_global[[1]],1,sum) + apply(C26233T_sims[[1]], 1, sum))) %>% 
  pivot_longer(c("d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d91"))
G26233T_df_nl$Variant = "G26233T"

sum_df1 <- rbind(C2416T_df, G26233T_df) 
sum_df2 <- rbind(C2416T_df_nl, G26233T_df_nl)
p10a <- ggplot(sum_df1, aes(x = name, y = value, fill = Variant, color = Variant)) + 
  geom_violin(adjust = 1.5, draw_quantiles = c(0.5)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95, angle = 90, size = 14), axis.title=element_text(size = 20), plot.title=element_text(hjust = 0.5)) + 
  xlab("") + ylab("Cases") +
  theme(legend.position = c(0.8, 0.8)) + 
  scale_x_discrete(labels = c("C2416T", "C2416T US (Before May 13 2020)", "C2416T US (After May 13 2020)", "C2416T (Time adjusted)",
                   "C26233T", "C26233T US (Before May 13 2020)", "C26233T US (After May 13 2020)", "C26233T (Time adjusted)", "C26233T global", "C26233T total"))+
  scale_fill_manual(values = cbPalette) + 
  scale_color_manual(values = cbPalette) +
  ggtitle("Binomial") + 
  coord_cartesian(ylim = c(0, 500000))
p10a

p10b <- ggplot(sum_df2, aes(x = name, y = value, fill = Variant, color = Variant)) + 
  geom_violin(adjust = 1.5, draw_quantiles = c(0.5)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95, angle = 90, size = 14), axis.title=element_text(size = 20), plot.title=element_text(hjust = 0.5)) + 
  xlab("") + ylab("Cases") +
  theme(legend.position = c(0.8, 0.8)) + 
  scale_x_discrete(labels = c("C2416T", "C2416T US (Before May 13 2020)", "C2416T US (After May 13 2020)", "C2416T (Time adjusted)", "C2416T r1", "C2416T r2",
                              "C26233T", "C26233T US (Before May 13 2020)", "C26233T US (After May 13 2020)", "C26233T (Time adjusted)", 
                              "C26233T r1","C26233T r2","G26233T global", "G26233T global r1", "G26233T global r2", "G26233T total"))+
  scale_fill_manual(values = cbPalette) + 
  scale_color_manual(values = cbPalette) +
  ggtitle("Normal Approximation") +
  coord_cartesian(ylim = c(0, 500000))
p10b


# create figures

plot_grid(p1, p2, p4, p5, labels = c("A", "B", "C", "D"), nrow = 4)
ggsave("data_revision/figures2/S15A-D.pdf", height = 16, width = 12)
ggsave("data_revision/figures2/S15A-D.jpg", height = 16, width = 12)

plot_grid(p11a, p12a, p13a, p11b, p12b,p13b, labels = c("E", "F", "G", "H", "I", "J"), nrow = 6)
ggsave("data_revision/figures2/S15E-J.pdf", height = 16, width = 12)
ggsave("data_revision/figures2/S15E-J.jpg", height = 16, width = 12)


plot_grid(p7, p8, labels = c("K", "L"), nrow = 1)
ggsave("data_revision/figures2/S15_K-L.pdf", height = 5, width = 10)
ggsave("data_revision/figures2/S15_K-L.jpg", height = 5, width = 10)

plot_grid(p10a, p10b, labels = c("M", "N"), nrow = 1, rel_widths = c(1, 1.4))
ggsave("data_revision/figures2/S15M-N.pdf", height = 8, width = 10)
ggsave("data_revision/figures2/S15M-N.jpg", height = 8, width = 10)


plot_grid(s0t, s2t, nrow = 2, labels = c("O", "P"), scale = 0.95)
ggsave("data_revision/figures2/S15_O-P.pdf", height = 12, width = 12)
ggsave("data_revision/figures2/S15_O-P.jpg", height = 12, width = 12 )


# binomial, time-dependent model, total cases C2416T
round(quantile(apply(td1[[6]], 1, sum), prob = c(0.05, 0.5, 0.95)), -3)

# binomial, time-dependent model, epoch 1, total cases C2416T
round(quantile(apply(td1[[4]], 1, sum), prob = c(0.05, 0.5, 0.95)), -3)

# binomial, time-dependent model, epoch 2, total cases C2416T
round(quantile(apply(td1[[5]], 1, sum), prob = c(0.05, 0.5, 0.95)), -3)

# proportion of MA cases, epoch 1, C2416T
quantile(td1[[4]][,12] / apply(td1[[4]], 1, sum), prob = c(0.05, 0.5, 0.95))

# binomial, time-dependent model, total cases G26233T
round(quantile(apply(td2[[6]], 1, sum), prob = c(0.05, 0.5, 0.95)), -3)

# binomial, time-dependent model, epoch 1, total cases G26233T
round(quantile(apply(td2[[4]], 1, sum), prob = c(0.05, 0.5, 0.95)), -3)

# binomial, time-dependent model, epoch 2, total cases G26233T
round(quantile(apply(td2[[5]], 1, sum), prob = c(0.05, 0.5, 0.95)), -3)

# proportion of FL cases, epoch 2, C2416T
quantile(td1[[5]][,8] / apply(td1[[5]], 1, sum), prob = c(0.05, 0.5, 0.95))

# proportion of FL cases, both epochs, C2416T
quantile(td1[[6]][,8] / apply(td1[[6]], 1, sum), prob = c(0.05, 0.5, 0.95))


