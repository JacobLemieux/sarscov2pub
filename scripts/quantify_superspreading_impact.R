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

setwd("~/Dropbox/COVID/SARS-CoV-2/")
source("scripts/SARS-CoV-2_functions.R")

# colorblind palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


# initialize list of simulations for total
state_totals <- list(length(10))

## Read in alignment and trim UTRs

#### Look at USA Data
#mafft_alnDNAss <-readDNAMultipleAlignment("~/Dropbox/COVID/2020_07_14/gisaid/subset_msa_0714.fasta", format="fasta")
#mafft_alnDNAss <-readDNAMultipleAlignment("~/Dropbox/COVID/2020_09_29/GISAID/subset2416_msa_0930.fasta", format="fasta")
mafft_alnDNAss <-readDNAMultipleAlignment("~/Dropbox/COVID/2020_11_02/msa_1102/subset2416_msa_1102.fasta", format="fasta")
#variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,148], sample_id = rownames(mafft_alnDNAss)))
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("~/Dropbox/COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

# filter

epoch1 <- "2020-05-30"

variant_association <- variant_association %>% 
  filter(variant_association$country == "USA")
  #filter(variant_association$country == "USA" & variant_association$date > epoch1)



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


cases <- variant_association %>% filter(variant_association$variant == "T") 

write_tsv(cases, "data_revision/metadata/C2416T_from_gisaid.tsv")


# case counts without division by epoch

doi <- "2020-11-01"
case_counts <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == doi)
case_counts$state <- sapply(case_counts$state, state_to_abbrev)
variants_by_state <- left_join(variants_by_state, case_counts, by = "state")
variants_by_state$linked_cases <- variants_by_state$AF * variants_by_state$cases
variants_by_state <- variants_by_state %>% filter(!is.na(cases))
sum(variants_by_state$linked_cases, na.rm=T)

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

# simulate number of cases based on estimated allele frequency and potential bias
variants_by_state$alleles <- variants_by_state$C + variants_by_state$T

n_sims <- 10000 


k <- c(194, 31, 12, 9) # number of C2416T alleles in each of top 4 counties
N <- c(418, 103, 40, 22) # number of samples in each of the top 4 counties
pop_size <- c(15423, 17694, 11429, 7172) # total number of cases reported through May 9 in each of top 4 counties
phat <- k/N
county_sims <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
county_sims_robust1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
county_sims_robust2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
old_mu_hat <- rep(NA, n_sims)
new_mu_hat <- rep(NA, n_sims)

# simulate only in those plates that have alleles

variants_by_state <- variants_by_state %>% filter(T > 0)

for(i in 1:nrow(variants_by_state)){
  mu_hat <- variants_by_state$AF[i]
  mu_hat.h <- (variants_by_state$T[i] + 0.5) / (variants_by_state$alleles[i] + 1)
  sd_hat <- sqrt((mu_hat.h*(1-mu_hat.h))/(variants_by_state$alleles[i] + 1))
  county_sims[,i] <- (rbinom(n = n_sims, size = variants_by_state$alleles[i], prob = rtnorm(n_sims, mu_hat,sd_hat, lower = 0, upper = 1) )/variants_by_state$alleles[i] *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1)) * variants_by_state$cases[i]
  county_sims_robust1[,i] <- (rbinom(n = n_sims, size = variants_by_state$alleles[i], prob = rtnorm(n_sims, mu_hat,2*sd_hat, lower = 0, upper = 1) )/variants_by_state$alleles[i] *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1)) * variants_by_state$cases[i]
  county_sims_robust2[,i] <- (rbinom(n = n_sims, size = variants_by_state$alleles[i], prob = rtnorm(n_sims, mu_hat,3*sd_hat, lower = 0, upper = 1) )/variants_by_state$alleles[i] *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1)) * variants_by_state$cases[i]  
}


E_MA_sims <- apply(county_sims, 1, function(x) sum(x, na.rm=T))
robust1_sums <- apply(county_sims_robust1, 1, function(x) sum(x, na.rm=T))
robust2_sums <- apply(county_sims_robust2, 1, function(x) sum(x, na.rm=T))

quantile(E_MA_sims, c(0.05, 0.5, 0.95))
mean(E_MA_sims)
state_totals[[1]] <- E_MA_sims
state_totals[[4]] <- robust1_sums
state_totals[[5]] <- robust2_sums

E_MA_sims <- data.frame(Total_Cases = E_MA_sims)

sum(variants_by_state$AF*variants_by_state$cases, na.rm=T)

# filter by epoch

#case_counts <- read_csv("data_revision/us-states_nyt_2020_09_28.csv") %>% filter(date == doi)
#case_counts_epoch1 <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == epoch1)
#case_counts_epoch2 <- case_counts
#case_counts_epoch2$cases <- case_counts_epoch2$cases - case_counts_epoch1$cases
#case_counts_epoch2$state <- sapply(case_counts_epoch2$state, state_to_abbrev)
#variants_by_state <- left_join(variants_by_state, case_counts_epoch2, by = "state")
#variants_by_state$linked_cases <- variants_by_state$AF * variants_by_state$cases
#sum(variants_by_state$linked_cases, na.rm=T)



# Look at 26233


mafft_alnDNAss <-readDNAMultipleAlignment("~/Dropbox/COVID/2020_11_02/msa_1102/subset26233_msa_1102.fasta", format="fasta")
#variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,148], sample_id = rownames(mafft_alnDNAss)))
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("~/Dropbox/COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

# output list of cases in the world with variant
cases <- variant_association %>% filter(variant_association$variant == "T") 
write_tsv(cases, "data_revision/metadata/G26233T_from_gisaid.tsv")

variant_association <- variant_association[variant_association$country == "USA",]

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
#case_counts <- read_csv("data_revision/us-states_nyt_2020_09_28.csv") %>% filter(date == doi)
case_counts <- read_csv("data_revision/metadata/nyt_case_counts_2020_11_02.csv") %>% filter(date == doi)
case_counts$state <- sapply(case_counts$state, state_to_abbrev)
variants_by_state <- left_join(variants_by_state, case_counts, by = "state")
variants_by_state$linked_cases <- variants_by_state$AF * variants_by_state$cases
variants_by_state <- variants_by_state %>% filter(!is.na(cases))
sum(variants_by_state$linked_cases, na.rm=T)

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



#plot_grid(p1, p2, labels = c("A", "B"), nrow = 2)
#ggsave("data_revision/figures/impact_superspreading.pdf", height = 5, width = 12)
#ggsave("data_revision/figures/impact_superspreading.jpg", height = 5, width = 12)

# simulate only in those plates that have alleles

variants_by_state <- variants_by_state %>% filter(T > 0)


# simulate number of cases based on estimated allele frequency and potential bias
variants_by_state$alleles <- variants_by_state$G + variants_by_state$T

n_sims <- 10000 


k <- c(194, 31, 12, 9) # number of C2416T alleles in each of top 4 counties
N <- c(418, 103, 40, 22) # number of samples in each of the top 4 counties
pop_size <- c(15423, 17694, 11429, 7172) # total number of cases reported through May 9 in each of top 4 counties
phat <- k/N
county_sims <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)


for(i in 1:nrow(variants_by_state)){
  mu_hat <- variants_by_state$AF[i]
  mu_hat.h <- (variants_by_state$T[i] + 0.5) / (variants_by_state$alleles[i] + 1)
  sd_hat <- sqrt((mu_hat.h*(1-mu_hat.h))/(variants_by_state$alleles[i] + 1))
  county_sims[,i] <- (rbinom(n = n_sims, size = variants_by_state$alleles[i], prob = rtnorm(n_sims, mu_hat,sd_hat, lower = 0, upper = 1)) /variants_by_state$alleles[i]) * variants_by_state$cases[i]
}

E_MA_sims <- apply(county_sims, 1, function(x) sum(x, na.rm=T))
quantile(E_MA_sims, c(0.05, 0.5, 0.95))
mean(E_MA_sims)
state_totals[[2]] <- E_MA_sims
E_MA_sims <- data.frame(Total_Cases = E_MA_sims)
sum(variants_by_state$AF*variants_by_state$cases, na.rm=T)


# Look at global 26233


mafft_alnDNAss <-readDNAMultipleAlignment("~/Dropbox/COVID/2020_11_02/msa_1102/subset26233_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("~/Dropbox/COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
to_remove <- to_remove$sample_id
to_remove <- gsub("_", "-", to_remove)
variant_association <- variant_association[!(variant_association$seq_name %in% to_remove),]

# output list of cases in the world with variant
cases <- variant_association %>% filter(variant_association$variant == "T") 

#variant_association$state <- substr(variant_association$seq_name, 1, 2)
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
#case_counts <- read_csv("data_revision/us-states_nyt_2020_09_28.csv") %>% filter(date == doi)
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

# simulate number of cases based on estimated allele frequency and potential bias
variants_by_country$alleles <- variants_by_country$G + variants_by_country$T


n_sims <- 10000 


county_sims <- matrix(rep(NA, n_sims*nrow(variants_by_country)), nrow = n_sims)

for(i in 1:nrow(variants_by_country)){
  mu_hat <- variants_by_country$AF[i]
  mu_hat.h <- (variants_by_country$T[i] + 0.5) / (variants_by_country$alleles[i] + 1)
  sd_hat <- sqrt((mu_hat.h*(1-mu_hat.h))/(variants_by_country$alleles[i] + 1))
  county_sims[,i] <- (rbinom(n = n_sims, size = variants_by_country$alleles[i], prob = rtnorm(n_sims, mu_hat,sd_hat, lower = 0, upper = 1) )/variants_by_country$alleles[i] ) * variants_by_country$cases[i]
}

E_MA_sims <- apply(county_sims, 1, function(x) sum(x, na.rm=T))
quantile(E_MA_sims, c(0.05, 0.5, 0.95))
mean(E_MA_sims)
state_totals[[3]] <- E_MA_sims

E_MA_sims <- data.frame(Total_Cases = E_MA_sims)

p9 <- ggplot(E_MA_sims, aes(x = Total_Cases)) + geom_histogram(aes(y=..density..)) + theme_bw() + 
  xlab("Cases of G26233T") + 
  ylab("Density") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank())

# create plot of 
sum_df <- data.frame(C2416T = state_totals[[1]], C2416T_r1 = state_totals[[4]], C2416T_r2 = state_totals[[5]], G26233T_dom = state_totals[[2]], G26233T_int = state_totals[[3]], G26233T_tot = state_totals[[2]] + state_totals[[3]]) %>% 
  pivot_longer(c("C2416T", "C2416T_r1", "C2416T_r2","G26233T_dom", "G26233T_int", "G26233T_tot"))

p10 <- ggplot(sum_df, aes(x = name, y = value)) + 
  geom_violin(adjust = 1.5, draw_quantiles = c(0.5)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95, angle = 90, size = 14), axis.title=element_text(size = 20), plot.title=element_text(hjust = 0.5)) + 
  xlab("") + ylab("Cases") +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("C2416T (US)","C2416T US (robust 1)", "C2416T US (robust 2)","G26233T (US)", "G26233T (International)", "G26233T (Global Total)"))
p10
plot_grid(p1, p2, p4, p5, labels = c("A", "B", "C", "D"), nrow = 4)
ggsave("data_revision/figures2/S15A-D.pdf", height = 11, width = 12)
ggsave("data_revision/figures2/S15A-D.jpg", height = 11, width = 12)

  
plot_grid(p7, p8, p10, labels = c("E", "F", "G"), nrow = 1)
ggsave("data_revision/figures2/S15_E-G.pdf", height = 5, width = 12)
ggsave("data_revision/figures2/S15_E-G.jpg", height = 5, width = 12 )

lapply(state_totals, function(x) round(quantile(x, probs = c(0.05, 0.5, 0.95)), -3))

plot_grid(s0t, s2t, nrow = 2, labels = c("H", "I"), scale = 0.95)
ggsave("data_revision/figures2/S15_H-I.pdf", height = 12, width = 12)
ggsave("data_revision/figures2/S15_H-I.jpg", height = 12, width = 12 )




# Time-dependent model
mafft_alnDNAss <-readDNAMultipleAlignment("~/Dropbox/COVID/2020_11_02/msa_1102/subset2416_msa_1102.fasta", format="fasta")
variant_association <- as_tibble(data.frame(variant = as.matrix(mafft_alnDNAss)[,1], sample_id = rownames(mafft_alnDNAss)))
date_zero <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
for(dz in date_zero){
  variant_association$sample_id <- gsub(paste("2020-",dz,"-00", sep=""), paste("2020-",dz,"-15", sep=""), variant_association$sample_id)
}
variant_association <- separate(variant_association,col = sample_id, into = c("seq_name", "annotation_num", "date", "location"), remove = FALSE, sep = "\\|")
variant_association <- separate(variant_association,col = sample_id, into = c("organism", "country", "seq_name"), remove = FALSE, sep = "/")
variant_association$date <- as_date(variant_association$date)
to_remove <- read_tsv("~/Dropbox/COVID/SARS-CoV-2/data_revision/metadata/conf_blacklist.tsv")
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

# plot the raw data
vbs_e1 <- variants_by_state[,c(1,6:8, 16, 18)]
vbs_e1 <- cbind(vbs_e1, rep("Epoch 1", nrow(vbs_e1)))
names(vbs_e1) <- c("state","AF", "ci_l", "ci_u", "cases", "linked_cases", "Epoch")
vbs_e2 <- variants_by_state[,c(1,9:11,17, 19)]
vbs_e2 <- cbind(vbs_e2, rep("Epoch 2", nrow(vbs_e2)))
names(vbs_e2) <- c("state","AF", "ci_l", "ci_u", "cases", "linked_cases", "Epoch")
vbs_df <- rbind(vbs_e1, vbs_e2)

p11 <- ggplot(vbs_df, aes(x = state, y = AF, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Allele Frequency") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) + 
  ggtitle("C2416T")
p11

p12 <- ggplot(vbs_df, aes(x = state, y = linked_cases, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = cases*ci_l, ymax = cases*ci_u), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Cases") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  ggtitle("C2416T")
p12
plot_grid(p11, p12, nrow = 2)

names(variants_by_state)[1] <- "STATE"
state_key <- state.map[,c("STATE", "region")]
state_key <- state_key[!duplicated(state_key$STATE),]
s_df <- left_join(variants_by_state, state_key, by = "STATE")
s_df <- s_df[,c("region", "AF1")]
names(s_df)[2] <- "value"
s_df <- rbind(s_df, c("north dakota", "0"))
s_df$value <- as.numeric(s_df$value)
s0 <- state_choropleth(s_df, num_colors = 1, legend = "Allele Frequency")
s0t <- state_choropleth(s_df, num_colors = 1, title = "      C2416T", legend = "Allele Frequency")
s0t


# simulate number of cases based on estimated allele frequency and potential bias
variants_by_state$alleles_epoch1 <- variants_by_state$C0 + variants_by_state$T0
variants_by_state$alleles_epoch2 <- variants_by_state$C1 + variants_by_state$T1

n_sims <- 10000 
county_sims <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
county_sims_param_weighted <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
county_sims_param_weighted_epoch1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
county_sims_param_weighted_epoch2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)


for(i in 1:nrow(variants_by_state)){
  mu_hat.epoch1 <- variants_by_state$AF0[i]
  mu_hat.h.epoch1 <- (variants_by_state$T0[i] + 0.5) / (variants_by_state$alleles_epoch1[i] + 1) #continuity correction
  sd_hat.epoch1 <- sqrt(mu_hat.h.epoch1*(1-mu_hat.h.epoch1)/(variants_by_state$alleles_epoch1[i]+1))
  mu_hat.epoch2 <- variants_by_state$AF1[i]
  mu_hat.h.epoch2 <- (variants_by_state$T1[i] + 0.5) / (variants_by_state$alleles_epoch2[i] + 1) #continuity correction
  sd_hat.epoch2 <- sqrt(mu_hat.h.epoch2*(1-mu_hat.h.epoch2)/(variants_by_state$alleles_epoch2[i]+1))
  mu_hat.pooled <- (variants_by_state$T0[i] + variants_by_state$T1[i]) / (variants_by_state$T0[i] + variants_by_state$T1[i] + variants_by_state$C0[i] + variants_by_state$C1[i])
  mu_hat.h.pooled <- (variants_by_state$T0[i] + variants_by_state$T1[i] + 0.5) / (variants_by_state$T0[i] + variants_by_state$T1[i] + variants_by_state$C0[i] + variants_by_state$C1[i] + 1)
  sd_hat.pooled <- sqrt(mu_hat.h.pooled*(1 - mu_hat.h.pooled)/(variants_by_state$alleles_epoch1[i]+ variants_by_state$alleles_epoch2[i]+1))
  if(variants_by_state$alleles_epoch1[i] < 10){
    mu_hat.epoch1 <- mu_hat.pooled
    sd_hat.epoch1 <- sd_hat.pooled
  }
  if(variants_by_state$alleles_epoch2[i] < 10){
    mu_hat.epoch2 <- mu_hat.pooled
    sd_hat.epoch2 <- sd_hat.pooled
  }
  county_sims_param_weighted_epoch1[,i] <- (rbinom(n = n_sims, size = variants_by_state$alleles_epoch1[i], prob = rtnorm(n_sims, mu_hat.epoch1,sd_hat.epoch1, lower = 0, upper = 1) )/variants_by_state$alleles_epoch1[i] *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1)) * variants_by_state$cases_epoch1[i]
  county_sims_param_weighted_epoch2[,i] <- (rbinom(n = n_sims, size = variants_by_state$alleles_epoch2[i], prob = rtnorm(n_sims, mu_hat.epoch2,sd_hat.epoch2, lower = 0, upper = 1) )/variants_by_state$alleles_epoch2[i] *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1)) * variants_by_state$cases_epoch2[i]
  county_sims_param_weighted[,i] <- county_sims_param_weighted_epoch1[,i] + county_sims_param_weighted_epoch2[,i]
}


E_MA_sims <- apply(county_sims_param_weighted, 1, function(x) sum(x, na.rm=T))
quantile(E_MA_sims, c(0.05, 0.5, 0.95))
E_MA_sims <- apply(county_sims_param_weighted_epoch1, 1, function(x) sum(x, na.rm=T))
quantile(E_MA_sims, c(0.05, 0.5, 0.95))
E_MA_sims <- apply(county_sims_param_weighted_epoch2, 1, function(x) sum(x, na.rm=T))
quantile(E_MA_sims, c(0.05, 0.5, 0.95))

state_totals[[6]] <- apply(county_sims_param_weighted_epoch1, 1, function(x) sum(x, na.rm=T))
state_dists <- apply(county_sims_param_weighted_epoch1, 2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95), na.rm=T))
vbs_e1 <- cbind(vbs_e1, t(state_dists*(1/0.9)))
names(vbs_e1)[8:10] <- c("sim_l_ci", "sim_median", "sim_u_ci")
state_totals[[7]] <- apply(county_sims_param_weighted_epoch2, 1, function(x) sum(x, na.rm=T))
state_dists2 <- apply(county_sims_param_weighted_epoch2, 2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95), na.rm=T))
vbs_e2 <- cbind(vbs_e2, t(state_dists2*(1/0.9)))
names(vbs_e2)[8:10] <- c("sim_l_ci", "sim_median", "sim_u_ci")
state_totals[[8]] <- apply(county_sims_param_weighted, 1, function(x) sum(x, na.rm=T))
vbs_df <- rbind(vbs_e1, vbs_e2)

p13 <- ggplot(vbs_df, aes(x = state, y = sim_median, color = Epoch)) + geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = sim_l_ci, ymax = sim_u_ci), alpha = 0.3, position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("Cases") +
  theme(legend.text = element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
  theme(legend.position = c(0.2, 0.815), legend.title=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95)) +
  ggtitle("C2416T")
p13
plot_grid(p12, p13, nrow = 2)


sum_df <- data.frame(C2416T = state_totals[[1]], C2416T_r1 = state_totals[[4]], C2416T_r2 = state_totals[[5]], 
                     C2416T_t1 = state_totals[[6]], C2416T_t2 = state_totals[[7]], C2416T_t2m = state_totals[[8]], 
                     G26233T_dom = state_totals[[2]], G26233T_int = state_totals[[3]], G26233T_tot = state_totals[[2]] + state_totals[[3]])%>% 
  pivot_longer(c("C2416T", "C2416T_r1", "C2416T_r2", "C2416T_t1", "C2416T_t2", "C2416T_t2m","G26233T_dom", "G26233T_int", "G26233T_tot"))

p10 <- ggplot(sum_df, aes(x = name, y = value)) + 
  geom_violin(adjust = 1.5, draw_quantiles = c(0.5)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95, angle = 90, size = 14), axis.title=element_text(size = 20), plot.title=element_text(hjust = 0.5)) + 
  xlab("") + ylab("Cases") +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("C2416T (US)","C2416T US (robust 1)", "C2416T US (robust 2)","C2416T US (before May 13 2020)", "C2416T US (After May 13 2020)",
                               "C2416T US (time adjusted)", "G26233T (US)", "G26233T (International)", "G26233T (Global Total)"))
p10
#plot_grid(p1, p2, p4, p5, labels = c("A", "B", "C", "D"), nrow = 4)
lapply(state_totals, function(x) quantile(x, c(0.05, 0.5, 0.95)))

