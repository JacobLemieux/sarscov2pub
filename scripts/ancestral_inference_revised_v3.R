### Ancestral inference, revised with updated gisaid sampling
# November 17 2020
# Jacob E. Lemieux

library(devtools)
library(tidyverse)
library(treeio)
library(Biostrings)
library(ggtree)
library(beastio)
library(gt)
library(ape)
library(phytools)
library(jsonlite)
library(MASS)
library(fitdistrplus)
library(lmtest)
library(cowplot)
library(lubridate)
library(DescTools)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

setwd("~/Dropbox/COVID/SARS-CoV-2/")
source("scripts/SARS-CoV-2_functions.R")
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/opt/anaconda3/bin:/opt/anaconda3/condabin:/Users/jy17/bin/google-cloud-sdk/bin:/Users/lemieux/google-cloud-sdk/bin:/Users/lemieux/anaconda2/bin:/usr/local/bin::/Users/lemieux/anaconda2/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:", sep=":"))

# Read in set of nextstrain trees
# samples 154 and 66 failed and were removed from the dataset
sample_set <- read_tsv("data_revision/metadata/gisaid_bootstrap_2_mod.tsv") %>% filter(!is.na(auspice_input_json))
num_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
NA_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
SA_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
Asia_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
Oceania_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
Europe_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
Africa_imports <- matrix(0, nrow = 3, ncol = nrow(sample_set))
sample_size <- rep(0, 0)
poisson_gini <- rep(0,0)
clade_size_gini <- rep(0,0)
nb_gini <- rep(0, 0)
model_list <- list(length = nrow(sample_set))
log_reg_list <- list(length = nrow(sample_set))
lrt_list <- list(length = nrow(sample_set))
next_import_list <- list(length = nrow(sample_set))
for(j in 1:nrow(sample_set)){
  to_download <- sample_set[j,c("time_tree")]
  print(j)
  print(sample_set[j,1])
  at_split <- strsplit(as.character(to_download), "\\/")[[1]] # path of ancestral traits split by delimiter
  print(at_split)
  anc_traits_path <- paste("gs:/", at_split[3], at_split[4], at_split[5], at_split[6], "call-ancestral_traits", "all_samples_aligned.fasta_aligned.filtered.masked_timetree_ancestral_traits.json", sep="/")
  #system(paste("gsutil cp", to_download, "data_revision/sequences_and_trees/nextstrain/trees"))
  #system(paste("mv data_revision/sequences_and_trees/nextstrain/trees/all_samples_aligned.fasta_aligned.filtered.masked_timetree.nwk data_revision/sequences_and_trees/nextstrain/trees/B", sample_set[j,"entity:gisaid_bootstrap_2_id"], ".tree", sep=""))
  #system(paste("gsutil cp", anc_traits_path, "data_revision/sequences_and_trees/nextstrain/ancestral_traits"))
  #system(paste("mv data_revision/sequences_and_trees/nextstrain/ancestral_traits/all_samples_aligned.fasta_aligned.filtered.masked_timetree_ancestral_traits.json data_revision/sequences_and_trees/nextstrain/ancestral_traits/B", sample_set[j,"entity:gisaid_bootstrap_2_id"], ".json", sep=""))
  anc_times_path <- paste("gs:/", at_split[3], at_split[4], at_split[5], at_split[6], "call-refine_augur_tree", "all_samples_aligned.fasta_aligned.filtered.masked_branch_lengths.json", sep="/")
  #system(paste("gsutil cp", anc_times_path, "data_revision/sequences_and_trees/nextstrain/ancestral_times"))
  #system(paste("mv data_revision/sequences_and_trees/nextstrain/ancestral_times/all_samples_aligned.fasta_aligned.filtered.masked_branch_lengths.json data_revision/sequences_and_trees/nextstrain/ancestral_times/B", sample_set[j,"entity:gisaid_bootstrap_2_id"], ".json", sep=""))
  nextstrain_tree <- read.newick(paste("data_revision/sequences_and_trees/nextstrain/trees/B", sample_set[j, "entity:gisaid_bootstrap_2_id"], ".tree", sep=""))
  anc_trait <- read_json(paste("data_revision/sequences_and_trees/nextstrain/ancestral_traits/B", sample_set[j, "entity:gisaid_bootstrap_2_id"], ".json", sep=""), simplifyVector = TRUE)
  anc_times <- read_json(paste("data_revision/sequences_and_trees/nextstrain/ancestral_times/B", sample_set[j, "entity:gisaid_bootstrap_2_id"], ".json", sep=""), simplifyVector = TRUE)
  sample_size <- c(sample_size, sample_set[j, "entity:gisaid_bootstrap_2_id"])
  
  
  #trans_matrix <- anc_trait$models$geoloc_ma$transition_matrix
  #rownames(trans_matrix) <- anc_trait$models$geoloc_ma$alphabet[1:nrow(trans_matrix)]
  #colnames(trans_matrix) <- anc_trait$models$geoloc_ma$alphabet[1:nrow(trans_matrix)]
  #heatmap(trans_matrix, Rowv=NA, Colv=NA)
  
  # flatten list to create a dataframe listing the taxa and trait inference
  
  anc_df <- data.frame(sample_id = names(anc_trait$nodes), country = rep(NA, length(anc_trait$nodes)), geoloc_ma = rep(NA, length(anc_trait$nodes)), 
                       MA = rep(NA, length(anc_trait$nodes)), MA_confNO = rep(NA, length(anc_trait$nodes)), MA_confYES = rep(NA, length(anc_trait$nodes)), 
                       EuropeConf = rep(NA, length(anc_trait$nodes)), AsiaConf = rep(NA, length(anc_trait$nodes)), AfricaConf = rep(NA, length(anc_trait$nodes)),
                       NAConf = rep(NA, length(anc_trait$nodes)), OceaniaConf = rep(NA, length(anc_trait$nodes)), MAConf = rep(NA, length(anc_trait$nodes)),
                       EUSAConf = rep(NA, length(anc_trait$nodes)), SAConf = rep(NA, length(anc_trait$nodes)), NEUSAConf = rep(NA, length(anc_trait$nodes)))
  for(i in 1:length(anc_trait$nodes)){
    #  anc_df$country[i] <- anc_trait$nodes[[i]]$country
    anc_df$geoloc_ma[i] <- anc_trait$nodes[[i]]$geoloc_ma
    anc_df$MA[i] <- anc_trait$nodes[[i]]$mass_binary
    if(!is.null(anc_trait$nodes[[i]]$mass_binary_confidence$NO)){
      anc_df$MA_confNO[i] <- anc_trait$nodes[[i]]$mass_binary_confidence$NO
    }
    if(!is.null(anc_trait$nodes[[i]]$mass_binary_confidence$YES)){
      anc_df$MA_confYES[i] <- anc_trait$nodes[[i]]$mass_binary_confidence$YES
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$Europe)){
      anc_df$EuropeConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$Europe
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$Asia)){
      anc_df$AsiaConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$Asia
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$Oceania)){
      anc_df$OceaniaConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$Oceania
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$Africa)){
      anc_df$AfricaConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$Africa
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$Massachusetts)){
      anc_df$MAConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$Massachusetts
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$`North America`)){
      anc_df$NAConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$`North America`
    }
    if(!is.null(anc_trait$nodes[[i]]$geoloc_ma_confidence$`South America`)){
      anc_df$SAConf[i] <- anc_trait$nodes[[i]]$geoloc_ma_confidence$`South America`
    }
  }
  
  # fill in gaps of confidence
  anc_df$MA_confYES[is.na(anc_df$MA_confYES)] <- 0
  
  #anc_df <- anc_df[ !(anc_df$MA_confNO > 0.1 & anc_df$MA_confNO < 0.9) ,] 
  #anc_df <- anc_df[!is.na(anc_df$sample_id),]
  model_list[[j]] <- anc_trait$models
  
  taxon_tbl <- as_tibble(nextstrain_tree)
  anc_df <- left_join(taxon_tbl, anc_df, by = c("label" = "sample_id"))
  anc_df$istip  <- !(substr(anc_df$label, 1, 4) == "NODE")
  
  # add in ancestral times
  anc_df$numdate <- rep(NA, nrow(anc_df))
  for(i in 1:nrow(anc_df)){
    anc_df$numdate[i] <- anc_times$nodes[[anc_df$label[i]]]$numdate
  } 
  
  # calculate state switches
  next_imports <- calc_MAimports2(anc_df, 0.9, nextstrain_tree)
  next_imports <- next_imports[order(next_imports$numdate),]
  
  
  # count imported nodes, combine with imported tips, produce fig 2D and regressions
  ## calculate import probability over time
  # find all the singleton imports, record their date of collection, record them as "imported"
  # find all non-singleton imports, record the TMRCA for the imported node, record them as "imported"
  # fine all MA samples that are not singleton imports, record their date of collection, record them as "non-imported"
  
  #beast_tibble <- as_tibble(mcc_beast)
  #anc_df is this
  # find the number of tips that are classified as imported
  table(anc_df$MA[nextstrain_tree %>% isTip(1:nrow(anc_df))])
  # OK, this is 791, good.
  
  # pull out tips that are from MA
  MA_df <- anc_df %>% filter(istip == TRUE & MA == "YES" & numdate < 2020.363)
  
  # mark tips as imported or not
  MA_df$imported <- MA_df$node %in% next_imports$node
  
  # pull out the nodes that correspond to imports
  next_imported_nodes <- next_imports$node[next_imports$clade_size > 1]
  
  # subset beast_tibble by these nodes
  next_imported_nodes_tibble <- anc_df[anc_df$node %in% next_imported_nodes,]
  
  next_imported_nodes_tibble$imported <- TRUE
  for(i in 1:length(next_imported_nodes)){
    subset_df <- anc_df[getDescendants(nextstrain_tree, next_imported_nodes[i]),]
    subset_df <- subset_df %>% filter(istip == TRUE & MA == "YES")
    imported_one <- which(subset_df$numdate == min(subset_df$numdate))[1]
    imported_node <- subset_df$node[imported_one]
    # flip one from domestic to imported
    MA_df$imported[MA_df$node == imported_node] <- TRUE
  }
  #next_MA_imports <- rbind(MA_df, next_imported_nodes_tibble)
  next_MA_imports <- MA_df # no longer add in the imported nodes, given that one is flipped
  next_MA_imports$date <- date_decimal(next_MA_imports$numdate)
  next_MA_imports$binary_MA <- ifelse(next_MA_imports$imported == TRUE, 1, 0)
  next_MA_imports$epoch_num <- rep(NA, nrow(next_MA_imports))
  next_MA_imports$epoch_num[next_MA_imports$date  < "2020-03-28"] <- 1
  next_MA_imports$epoch_num[next_MA_imports$date  >= "2020-03-28" & next_MA_imports$date < "2020-04-15"] <- 2
  next_MA_imports$epoch_num[next_MA_imports$date >= "2020-04-15"] <- 3
  log_reg_list[[j]] <- summary(glm(next_MA_imports$binary_MA ~ next_MA_imports$date, family= "binomial"))
  if(j == 1){
    p0 <- ggplot(next_MA_imports, aes(x = date, y = binary_MA)) + geom_jitter(height = 0.02, alpha = 0.1) + theme_bw() + 
      stat_smooth(se=FALSE, alpha = 0.2, span = 1.5, geom = 'line') + 
      stat_smooth(method="glm", method.args = list(family = "binomial"), se=FALSE, alpha = 0.2, color = "red", geom = 'line')
  }else{
    p0 <- p0 + geom_jitter(data = next_MA_imports, aes(x = date, y = binary_MA), width = 0.5, height = 0.02, alpha = 0.1) + 
      stat_smooth(data = next_MA_imports, se=FALSE, alpha = 0.2, color = "blue", geom='line', span = 1.5) + 
      stat_smooth(data = next_MA_imports, method="glm", method.args = list(family = "binomial"), se=FALSE, alpha = 0.2, color = "red", geom = 'line')
    #print(p0)
  }
  
  # Figure 2D
  
  next_imports_by_epoch <- table(next_MA_imports$binary_MA, next_MA_imports$epoch_num)
  next_imports_by_epoch <- data.frame(Imported = next_imports_by_epoch[2,], Total = colSums(next_imports_by_epoch), Proportion_Imported = next_imports_by_epoch[2,]/colSums(next_imports_by_epoch))
  next_imports_by_epoch$upper_ci <- rep(NA, nrow(next_imports_by_epoch))
  next_imports_by_epoch$lower_ci <- rep(NA, nrow(next_imports_by_epoch))
  next_imports_by_epoch$epoch <- c("Early", "Middle", "Late")
  next_imports_by_epoch$epoch <- factor(next_imports_by_epoch$epoch, levels = c("Early", "Middle", "Late"))
  for(i in 1:nrow(next_imports_by_epoch)){
    prop_test <- prop.test(next_imports_by_epoch$Imported[i], next_imports_by_epoch$Total[i])
    next_imports_by_epoch$upper_ci[i] <- prop_test$conf.int[1]
    next_imports_by_epoch$lower_ci[i] <- prop_test$conf.int[2]
  }
  
  next_imports_by_epoch[,3:5] <- round(next_imports_by_epoch[,3:5], 2)
  if(j == 1){
    
    p1 <- ggplot(next_imports_by_epoch, aes(x = epoch, y = Proportion_Imported), alpha = 0.1) + geom_point(size = 4, alpha = 0.2) + 
      #geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, alpha = 0.2) + 
      theme_bw() + 
      theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) +
      theme(legend.text = element_text(size=14))+
      scale_color_manual(values=cbPalette) + 
      ylab("Fraction Imported") + xlab("Time Period") 
  }else{
    p1 <- p1 + geom_jitter(data = next_imports_by_epoch, aes(x = epoch, y = Proportion_Imported), size = 4, alpha = 0.2, width = 0.1) #+ 
    #geom_errorbar(data = next_imports_by_epoch, aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, alpha = 0.2)
  }
  next_imports$date <- date_decimal(next_imports$numdate)
  pois_glm <- glm(next_imports$clade_size ~ next_imports$numdate, family = "poisson")
  nb_glm <- glm.nb(next_imports$clade_size ~ next_imports$numdate)
  lrt_list[[j]] <- lrtest(pois_glm, nb_glm)
  
  next_imports$phat <- predict(pois_glm, type="response")
  next_imports$nbhat <- predict(nb_glm, type="response")
  
  #table(next_imports$clade_size)
  next_imports$imported_before <- cumsum(next_imports$clade_size)
  next_imports$import_num <- 1:nrow(next_imports)
  
  if(j == 1){
    p2 <- ggplot(next_imports, aes(date, clade_size)) + 
      theme_bw() +
      geom_point(data = next_imports, aes(x = date, y = clade_size), alpha = 0.1) + 
      geom_line(data = next_imports, aes(x = date, y = phat), color = 2, alpha = 0.1) +
      geom_line(data = next_imports, aes(x = date, y = nbhat), color = 3, alpha = 0.1)
    #print(p2)
  }else{
    p2 <- p2 + geom_point(data = next_imports, aes(x = date, y = clade_size), alpha = 0.1) + 
      geom_line(data = next_imports, aes(x = date, y = phat), color = 2, alpha = 0.1) +
      geom_line(data = next_imports, aes(x = date, y = nbhat), color = 3, alpha = 0.1)
    #print(p2)
  }
  
  next_poisson_fit <- fitdistr(next_imports$clade_size, "poisson")
  poisson_sample <- rpois(length(next_imports$clade_size), lambda = next_poisson_fit$estimate)
  next_nb_fit <- fitdist(next_imports$clade_size, "nbinom")
  nb_sample <- rnbinom(length(next_imports$clade_size), size = next_nb_fit$estimate["size"], mu = next_nb_fit$estimate["mu"])
  poisson_gini <- c(poisson_gini, Gini(poisson_sample))
  nb_gini <- c(nb_gini, Gini(nb_sample))
  clade_size_gini <- c(clade_size_gini, Gini(next_imports$clade_size))
  
  
  
  #print(p1)
  #px <- plot_grid(p0, p1, labels = c("A", "B"))
  #print(px)
  #ggplot(next_MA_imports, aes(x = epoch_num, y = binary_MA)) + geom_jitter(height = 0.02) + theme_bw() + geom_smooth()+ 
  #  stat_smooth(method="glm", method.args = list(family = "binomial"), se=TRUE)
  
  
  pois_glm <- glm(next_imports$clade_size ~ next_imports$numdate, family = "poisson")
  nb_glm <- glm.nb(next_imports$clade_size ~ next_imports$numdate)
  lrtest(pois_glm, nb_glm)
  
  next_imports$phat <- predict(pois_glm, type="response")
  next_imports$nbhat <- predict(nb_glm, type="response")
  
  # calculate state switches for regional model
  reg_imports <- calc_Regionimports(anc_df, nextstrain_tree)
  
  # calculate imports by epoch
  next_imports$epoch_num <- rep(NA, nrow(next_imports))
  next_imports$epoch_num[next_imports$date  < "2020-03-28"] <- 1
  next_imports$epoch_num[next_imports$date  >= "2020-03-28" & next_imports$date < "2020-04-15"] <- 2
  next_imports$epoch_num[next_imports$date >= "2020-04-15"] <- 3
  
  num_imports[,j] <- table(next_imports$epoch_num)
  next_imports$date <- date_decimal(next_imports$numdate)
  next_imports$imported_before <- cumsum(next_imports$clade_size)
  next_imports$import_num <- 1:nrow(next_imports)
  next_import_list[[j]] <- next_imports
  
  reg_imports <- left_join(reg_imports, anc_df[,c("node", "numdate")], by = "node")
  reg_imports$date <- date_decimal(reg_imports$numdate)
  reg_imports$epoch_num <- rep(NA, nrow(reg_imports))
  reg_imports$epoch_num[reg_imports$date  < "2020-03-28"] <- 1
  reg_imports$epoch_num[reg_imports$date  >= "2020-03-28" & reg_imports$date < "2020-04-15"] <- 2
  reg_imports$epoch_num[reg_imports$date >= "2020-04-15"] <- 3
  if(nrow(table(reg_imports$parent_state == "North America", reg_imports$epoch_num)) == 2){
    NA_imports[,j] <- table(reg_imports$parent_state == "North America", reg_imports$epoch_num)[2,]
  }
  if(nrow(table(reg_imports$parent_state == "South America", reg_imports$epoch_num)) == 2){
    SA_imports[,j] <- table(reg_imports$parent_state == "South America", reg_imports$epoch_num)[2,]
  }
  if(nrow(table(reg_imports$parent_state == "Asia", reg_imports$epoch_num)) == 2){
    Asia_imports[,j] <- table(reg_imports$parent_state == "Asia", reg_imports$epoch_num)[2,]
  }
  if(nrow(table(reg_imports$parent_state == "Europe", reg_imports$epoch_num)) == 2){
    Europe_imports[,j] <- table(reg_imports$parent_state == "Europe", reg_imports$epoch_num)[2,]
  }
  if(nrow(table(reg_imports$parent_state == "Oceania", reg_imports$epoch_num)) == 2){
    Oceania_imports[,j] <- table(reg_imports$parent_state == "Oceania", reg_imports$epoch_num)[2,]
  }
  if(nrow(table(reg_imports$parent_state == "Africa", reg_imports$epoch_num)) == 2){
    Africa_imports[,j] <- table(reg_imports$parent_state == "Africa", reg_imports$epoch_num)[2,]
  }
  
}

reg_import_list <- list(Binary = num_imports, Africa = Africa_imports, Asia = Asia_imports, Europe = Europe_imports, NAmerica = NA_imports, Oceania = Oceania_imports, SAmerica = SA_imports)
reg_median <- matrix(nrow = 3, ncol = length(reg_import_list))
reg_lci <- matrix(nrow = 3, ncol = length(reg_import_list))
reg_uci <- matrix(nrow = 3, ncol = length(reg_import_list))
for(i in 1:ncol(reg_median)){
  reg_median[,i] <- apply(reg_import_list[[i]], 1, median)
  reg_lci[,i] <- apply(reg_import_list[[i]], 1, function(x) quantile(x, prob = c(0.05), type = 1))
  reg_uci[,i] <- apply(reg_import_list[[i]], 1, function(x) quantile(x, prob = c(0.95), type = 1))
}

table1 <- matrix(nrow = 3, paste(reg_median," (", reg_lci,"-", reg_uci, ")", sep=""))
table1 <- data.frame(t(rbind(c("Non-MA", "Africa", "Asia", "Europe", "North America", "Oceania", "South America"), table1)))
colnames(table1) <- c("Region", "Before March 28  ", "March 28 - April 15", "After April 15")
t0 <- table1 %>% gt() %>% tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels(1:4))
t0 %>%  tab_row_group(group = "Regional model", rows = 1:7)%>% tab_row_group(group = "Binary model", rows = 1:1) %>% cols_align(align = "center", columns = 2:4)
#gtsave(t1, "data/tables/T1.pdf")


p0 <- p0 + theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) +
  theme(legend.text = element_text(size=14))+
  scale_color_manual(values=cbPalette) + 
  ylab("Fraction Imported") + xlab("Date") 
p0

log_reg_coeffs <- sapply(log_reg_list, function(x) x$coefficients[6])
log_reg_p <- sapply(log_reg_list, function(x) x$coefficients[8])
round(quantile(log_reg_coeffs, probs = c(0.05, 0.5, 0.95), type = 1), 2)
quantile(log_reg_p, probs = c(0.25, 0.5, 0.75), type = 1)

p2 <- p2 + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) +
  theme(legend.text = element_text(size=12))+
  scale_color_manual(values=cbPalette) + 
  ylab("Imported Clade Size") + xlab("Date")
p2
lrt_p <- sapply(lrt_list, function(x) x$`Pr(>Chisq)`[2])

plot_grid(p0, p2,p3, labels = c("A", "B", "C"), label_size = 16, rel_widths = c(2,2,1), nrow = 1)

ggsave("data_revision/figures2/S9.pdf", height = 5, width = 10)
ggsave("data_revision/figures2/S9.jpg", height = 5, width = 10)



p1
ggsave("data_revision/figures2/Fig2D.pdf", height = 5, width = 5)
ggsave("data_revision/figures2/Fig2D.jpg", height = 5, width = 5)

gini_df <- data.frame(Observed = clade_size_gini, Poisson = poisson_gini, Negative_Binomial = nb_gini) %>% pivot_longer(c("Observed", "Poisson", "Negative_Binomial"), names_to = "Dist_Type", values_to = "Gini")

p3 <- ggplot(gini_df, aes(x = reorder(Dist_Type, Gini, na.rm=T), y = Gini)) + geom_jitter(width = 0.2, alpha = 0.2) + #geom_violin(adjust = 5) + #geom_density(adjust = 2) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14), axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  theme(legend.text = element_text(size=12))+
  scale_color_manual(values=cbPalette) + 
  ylab("Gini Coefficient") + xlab("") + 
  theme(legend.position = c(0.25, 0.85), legend.title=element_blank()) +
  scale_x_discrete(labels = c("Poisson", "Neg Binomial", "Observed"))
p3




# calculate percentage of imports that are due to large events
percent_cases_singletons <- sapply(next_import_list, function(x) sum(x$clade_size[x$clade_size < 2])/sum(x$clade_size))
percent_cases_not_singletons <- sapply(next_import_list, function(x) sum(x$clade_size[x$clade_size > 1])/sum(x$clade_size))
percent_imports_singletons <- sapply(next_import_list, function(x) sum(x$clade_size < 2)/nrow(x))
percent_imports_not_singletons <- sapply(next_import_list, function(x) sum(x$clade_size > 1)/nrow(x))
percent_cases_top10 <- sapply(next_import_list, function(x) sum(x$clade_size[order(x$clade_size, decreasing = T)][1:10])/sum(x$clade_size))
percent_imports_top10 <- sapply(next_import_list, function(x) 10/nrow(x))

# calculate quantiles

case_stats_df <- data.frame(Cases_Linked_to_Singletons = percent_cases_singletons, Cases_non_singletons = percent_cases_not_singletons, Singleton_imports = percent_imports_singletons, Not_singleton_imports = percent_imports_not_singletons)

apply(case_stats_df, 2, function(x) quantile(x, prob = c(0.05, 0.5, 0.95)))

percent_cases <- data.frame(Cases_Linked_to_Singletons = percent_cases_singletons, Cases_non_singletons = percent_cases_not_singletons) %>% 
  pivot_longer(cols = c("Cases_Linked_to_Singletons", "Cases_non_singletons"), values_to = c("Proportion"), names_to = c("Type"))
p4 <- ggplot(percent_cases, aes(x = reorder(Type, Proportion), y = Proportion, color = Type)) + geom_jitter(alpha  = 0.2) +  #+ geom_density(adjust = 2) + 
  theme_bw() + 
  expand_limits(x = c(0,1)) + 
  #coord_cartesian(xlim = c(0, 1)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14), axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  theme(legend.text = element_text(size=12))+
  scale_color_manual(values=cbPalette) + 
  ylab("Proportion of Cases") + xlab("") + 
  scale_x_discrete(labels = c("Singletons", "Non-singletons")) + 
  theme(legend.position = c(0.60, 0.85), legend.title=element_blank()) + theme(legend.position = "none")
p4

percent_imports <- data.frame(Singleton_imports = percent_imports_singletons, Not_singleton_imports = percent_imports_not_singletons) %>% 
  pivot_longer(cols = c("Singleton_imports", "Not_singleton_imports"), values_to = c("Proportion"), names_to = c("Type"))

p5 <- ggplot(percent_imports, aes(x = reorder(Type, Proportion), y = Proportion, color = Type)) + geom_jitter(alpha  = 0.2) +  #+ geom_density(adjust = 2) + 
  theme_bw() + 
  expand_limits(x = c(0,1)) + 
  #coord_cartesian(xlim = c(0, 1)) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14), axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  theme(legend.text = element_text(size=12))+
  scale_color_manual(values=cbPalette[c(2,1)]) + 
  ylab("Proportion of Imports") + xlab("") + 
  scale_x_discrete(labels = c("Non-singletons", "Singletons")) + 
  theme(legend.position = c(0.60, 0.85), legend.title=element_blank()) + theme(legend.position = "none")
p5

plot_grid(p5, p4)
ggsave("data_revision/figures2/2E.pdf", height = 5, width = 6)
ggsave("data_revision/figures2/2E.jpg", height = 5, width = 6)




plot_grid(p4, p3, p2, labels=c("C", "D", "E"), nrow = 1)
ggsave("data_revision/figures2/S8_C-E.pdf", height = 5, width = 12)
ggsave("data_revision/figures2/S8_C-E.jpg", height = 5, width = 12)

# plot model parameters

model_rates <- sapply(model_list, function(x) x$mass_binary$rate)
region_rates <- sapply(model_list, function(x) x$geoloc_ma$rate)
non_MA <- sapply(model_list, function(x) x$mass_binary$transition_matrix[1,2])
africa_MA <- sapply(model_list, function(x) x$geoloc_ma$transition_matrix[4,1])
asia_MA <- sapply(model_list, function(x) x$geoloc_ma$transition_matrix[4,2])
Europe_MA <- sapply(model_list, function(x) x$geoloc_ma$transition_matrix[4,3])
na_MA <- sapply(model_list, function(x) x$geoloc_ma$transition_matrix[4,5])
oceania_MA <- sapply(model_list, function(x) x$geoloc_ma$transition_matrix[4,6])
sa_MA <- sapply(model_list, function(x) x$geoloc_ma$transition_matrix[4,7])

model_rate_df <- data.frame(model_rates = model_rates, region_rates = region_rates) %>% pivot_longer(c(model_rates, region_rates))

model_trans_df <- data.frame(AAnon_MA = non_MA, africa_MA = africa_MA, asia_MA = asia_MA, Europe_MA = Europe_MA, na_MA = na_MA, oceania_MA = oceania_MA, sa_MA = sa_MA) %>% 
  pivot_longer(c(AAnon_MA, africa_MA, asia_MA, Europe_MA, na_MA, oceania_MA, sa_MA))

p6 <- ggplot(model_rate_df, aes (x =  name, y = value, fill = name)) + geom_violin(draw_quantiles = c(0.5)) + theme_bw() + 
  scale_x_discrete(labels = c("Binary Model", "Region Model")) + 
  theme(axis.text.x = element_text(angle = -90, size = 14, vjust =0.5)) + 
  xlab("") + ylab("Rate") + 
  scale_fill_manual(values=cbPalette, labels = c("Binary model", "Regional model")) + theme(legend.position = c(0.2, 0.8)) + labs(fill = "")# + 
  #scale_fill_discrete(labels = c("Binary model", "Region model"))
p6  
p7 <- ggplot(model_trans_df, aes (x =  name, y = value, fill = name)) + geom_violin(draw_quantiles = c(0.5)) + theme_bw() + 
  scale_x_discrete(breaks = c("africa_MA", "asia_MA", "Europe_MA", "na_MA", "oceania_MA", "sa_MA", "AAnon_MA"), labels = c("Africa <-> MA", "Asia <-> MA", "Europe <-> MA", "NA <-> MA", "Oceania <-> MA", "SA <-> MA", "non-MA <-> MA")) + 
  theme(axis.text.x = element_text(angle = -90, size = 14, vjust =0.5)) + 
  xlab("") + ylab("Transition Matrix Value") + 
  scale_fill_manual(values =cbPalette[c(1,2,2,2,2,2,2)]) + theme(legend.position = "none")


plot_grid(p6,  p7, labels = c("C", "D"), label_size = 16, rel_widths = c(1,1,2))
ggsave("data_revision/figures2/S7_C_D.pdf", height = 6, width = 12)
ggsave("data_revision/figures2/S7_C_D.jpg", height = 6, width = 12)


# quickly recalculate the poisson and nb regression coefficients
pr_list <- list(length = nrow(sample_set))
nbr_list <- list(length = nrow(sample_set))

for(i in 1:length(next_import_list)){
  pr_list[[i]] <- glm(next_import_list[[i]]$clade_size ~ next_import_list[[i]]$numdate, family = "poisson")
  nbr_list[[i]] <- glm.nb(next_import_list[[i]]$clade_size ~ next_import_list[[i]]$numdate)
}

pr_reg_coeffs <- sapply(pr_list, function(x) x$coefficients[2])
nbr_reg_coeffs <- sapply(nbr_list, function(x) x$coefficients[2])
round(quantile(pr_reg_coeffs, probs = c(0.05, 0.5, 0.95), type = 1), 2)
round(quantile(nbr_reg_coeffs, probs = c(0.05, 0.5, 0.95), type = 1), 2)