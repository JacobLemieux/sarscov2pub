### Functions to manipulate SARS-CoV-2 data
# Dec 1 2020
# Jacob E. Lemieux


# convert REDCap table by ID to sample table
convert_REDCap_longer <- function(REDCap_df){
  RC_data <- REDCap_df[,c("record_id","death", "gender_1", "ethnicity", "date1", "s1rs", "s1broadid", "confa", "snfa", "zip_code_of_residence", "clinic")]
  names(RC_data) <- c("record_id", "death", "gender_1", "ethnicity", "date", "rs", "sample_id", "confa", "snfa", "zip_code_of_residence", "clinic")
  for(i in 2:6){
    tmp_df <- REDCap_df[,c("record_id","death", "gender_1", "ethnicity", "confa", "snfa","zip_code_of_residence", "clinic",  paste("date", i,sep=""), paste("s", i, "rs", sep=""), paste("s", i, "broadid", sep=""))]
    names(tmp_df) <- c("record_id","death", "gender_1", "ethnicity", "confa", "snfa","zip_code_of_residence", "clinic", "date", "rs",  "sample_id")
    RC_data <- rbind(RC_data, tmp_df)
  }
  RC_data <- RC_data[!is.na(RC_data$date),]
  RC_data$sample_id <- gsub("MAMGH", "MA_MGH_", RC_data$sample_id)
  RC_data$date <- as_date(RC_data$date)
  RC_data
}

# reads redcap data and converts to longer format (each sample is a row)
read_REDCap_data <- function(redcap_csv){
  clin_data <- read_csv(redcap_csv)
  clin_data <- convert_REDCap_longer(clin_data)
  clin_data
}

# convert redcap to metadata

convert_REDCap_to_metadata <- function(clin_data){
  clin_metadata <- clin_data[,c("sample_id", "date", "confa", "snfa", "zip_code_of_residence", "clinic")]
  clin_metadata$exposure <- rep("NoKnown", nrow(clin_metadata))
  clin_metadata$exposure[clin_metadata$confa == 1] <- "ConfA"
  clin_metadata$exposure[clin_metadata$snfa == 1 | clin_metadata$snfa == 2] <- "SNFA"
  clin_metadata <- unite(clin_metadata, label, c(sample_id, exposure, date), sep="|", remove = FALSE)
  clin_metadata[order(clin_metadata$sample_id), c("label", "sample_id", "date", "exposure", "zip_code_of_residence", "clinic")]
}

# read DPH data
read_DPH <- function(dph_csv){
  dph_df <- read_csv(dph_csv, na = "N/A")
  names(dph_df) <- c("sample_id", "DPH_acc", "DPH_provider", "DPH_date","specimen_type", "DPH_N1", "DPH_N2", "DPH_RP")
  dph_df$sample_id <- gsub("SC2-MA-", "MA_DPH_00", dph_df$sample_id)
  dph_df$DPH_provider[dph_df$DPH_provider == "UMASS BOSTON HEALTH SERVICES" | dph_df$DPH_provider == "EMERSON HOSPITAL"] <- "other"
  dph_df$DPH_provider[dph_df$DPH_provider == "BERKSHIRE MEDICAL CENTER"] <- "BerkshireMed"
  dph_df$DPH_provider[dph_df$DPH_provider == "BOSTON HEALTH CARE FOR THE HOMELESS"] <- "BHCHP"
  dph_df <- unite(dph_df, label, c(sample_id, DPH_provider, DPH_date), sep="|", remove=FALSE)
  dph_df
}

read_DPH2 <- function(dph_csv){
  dph_df <- read_csv(dph_csv, na = "N/A")
  names(dph_df) <- c("sample_id", "DPH_date", "DPH_N1", "DPH_N2", "DPH_RP")
  dph_df$sample_id <- gsub("SC2-MA-", "MA_DPH_00", dph_df$sample_id)
  dph_df$DPH_provider <- rep("NoKnown", nrow(dph_df))
  dph_df$case_type <- c(rep("Travel_associated", 3), rep("Western_MA", 8), rep("Homeless", nrow(dph_df)-11))
  dph_df$DPH_acc <- c(1:nrow(dph_df))
  dph_df <- unite(dph_df, label, c(sample_id, case_type, DPH_date), sep="|", remove=FALSE)
  dph_df
}

merge_metadata <- function(dph_metadata, mgh_metadata){
  dph_metadata <- dph_metadata[,c("label", "sample_id", "DPH_date", "case_type")]
  dph_metadata$zip_code_of_residence <- rep(NA, nrow(dph_metadata))
  dph_metadata$clinic <- rep(NA, nrow(dph_metadata))
  names(dph_metadata)[3:4] <- c("date", "exposure")
  mgh_metadata <- mgh_metadata[,1:6]
  rbind(mgh_metadata, dph_metadata)
}

# read compiled metadata file (eventually this should merge with REDcap csv)

read_compiled_metadata <- function(metadata_tsv){
  ss_m <- read_tsv(metadata_tsv)
  ss_m <- unite(ss_m, label, c(sample_id, exposure, date), sep="|", remove=FALSE)
  ss_m
}

# clean up sample set from Terra

clean_sample_set <- function(ss){
  names(ss)[1] <- "sample_id" # rename column 1 to make it easier to merge with other datasets
  ss <- ss[ss$QC == "PASS",] # only consider samples that have passed
  ss <- ss[grep("MA", ss$sample_id),] # remove empty lines 
  ss
}

# Filter for high quality genomes

filter_sample_set <- function(ss, assembly_percent_complete){
  reference_length = 29903
  ss_hq <- ss[ss$assembly_length_unambiguous > (assembly_percent_complete * reference_length),]
  ss_hq
}

# Filter for unique genomes

uniquify_sample_set <- function(ss, redcap_csv){
  ss_unique <- ss[ss$sample_id %in% redcap_csv$s1broadid,] # include only samples that have complete genomes from first sample
  ss_unique <- rbind(ss_unique, ss[ss$sample_id %in% redcap_csv$s2broadid & !(ss$sample_id %in% ss_unique$sample_id),]) # include samples that had complete genomes in subsequent sample but not yet included
  ss_unique <- rbind(ss_unique, ss[ss$sample_id %in% redcap_csv$s3broadid & !(ss$sample_id %in% ss_unique$sample_id),]) # needs to be converted to a for loop...
  ss_unique <- rbind(ss_unique, ss[ss$sample_id %in% redcap_csv$s4broadid & !(ss$sample_id %in% ss_unique$sample_id),])
  ss_unique <- rbind(ss_unique, ss[ss$sample_id %in% redcap_csv$s5broadid & !(ss$sample_id %in% ss_unique$sample_id),])
  ss_unique <- rbind(ss_unique, ss[ss$sample_id %in% redcap_csv$s6broadid & !(ss$sample_id %in% ss_unique$sample_id),])
  ss_unique
}

# Read and clean viral Ct

read_viral_ct <- function(ct_csv){
  v_ct <- read_csv(ct_csv)
  v_ct <- v_ct[!v_ct$is.control,]
  names(v_ct)[1] <- "sample_id"
  v_ct$ViralCT[v_ct$ViralCT  == 45] <- NA
  v_ct
}

# read Pangolin lineage report CSV
read_pangolin <- function(lin_report){
  pang_lin <- read_csv(lin_report)
  #pang_lin <- separate(pang_lin, col = taxon, sep = "\\|", into = c("sample_id", "exposure_group", "sample_date"))
  pang_lin <- pang_lin[,1:2]
  names(pang_lin)[1] <- "sample_id"
  pang_lin
}


# read Pangolin lineage report CSV with version 2 lineages
read_pangolinv2 <- function(lin_report){
  pang_lin <- read_csv(lin_report)
  #pang_lin <- separate(pang_lin, col = taxon, sep = "\\|", into = c("sample_id", "exposure_group", "sample_date"))
  pang_lin <- pang_lin[,1:2]
  names(pang_lin)[1] <- "sample_id"
  pang_lin
}



# download fastas from google cloud

pull_fastas_and_concatenate <- function(ss,fasta_path = "data/fasta", mfa_path = "data/sequences_and_trees/concatenated_unaligned.fasta", do_not_download = TRUE, clean_dir = FALSE){
  if(do_not_download == FALSE){
    if(clean_dir==TRUE){
      system(paste("rm ", fasta_path, "/*.fasta", sep=""))
    }
    for(i in 1:nrow(ss)){
      system(paste("gsutil cp ", ss$assembly_fasta[i], fasta_path))
    }
    system(paste("cat ",fasta_path, "/*.fasta > ", mfa_path, sep=""))
    print(paste("placed concatenated fastas in ", mfa_path, sep=""))
  }
  fa_seqs <- readDNAStringSet(mfa_path, format="fasta")
  fa_seqs
}

# rename alignment with metadata
# for some reason the %in% operator for vectors is not cooperating with the DNAStringSet
relabel_DNAss <- function(DNA_string_set, ss_metadata){
  for(i in 1:length(names(DNA_string_set))){
    if(names(DNA_string_set)[i] %in% ss_metadata$sample_id){
      names(DNA_string_set)[i] <- ss_metadata$label[ss_metadata$sample_id == names(DNA_string_set)[i]]
    }
  }
  DNA_string_set
}

# read and reformat case data

read_and_parse_DPH_cases <- function(cases_csv, ss_metadata, ss_unique_hq){
  case_counts <- read.csv(cases_csv)
  case_counts$Cases <- case_counts$Cases + 1   # correct for MA_1
  MGH_dates_df <- data.frame(table(ss_metadata$date))
  colnames(MGH_dates_df) <- c("Date", "Cases")
  MGH_dates_sequenced <- ss_metadata[ss_metadata$sample_id %in% ss_unique_hq$sample_id,]
  MGH_dates_df_sequenced <- data.frame(table(MGH_dates_sequenced$date))
  colnames(MGH_dates_df_sequenced) <- c("Date", "Cases")
  case_counts$MA_sequenced <- rep(0, nrow(case_counts))
  case_counts$MA_attempted <- rep(0, nrow(case_counts))
  case_counts$MA_attempted[1] <- 1 # correct for MA_1
  case_counts$MA_sequenced[1] <- 1 # correct for MA_1
  for(i in 1:nrow(case_counts)){
    if(sum(MGH_dates_df_sequenced$Date %in% case_counts$Date[i]) > 0){
      case_counts$MA_sequenced[i] <- MGH_dates_df_sequenced$Cases[MGH_dates_df_sequenced$Date %in% case_counts$Date[i]]
    }
    if(sum(MGH_dates_df$Date %in% case_counts$Date[i]) > 0){
      case_counts$MA_attempted[i] <- MGH_dates_df$Cases[MGH_dates_df$Date %in% case_counts$Date[i]]
    }
  }
  case_counts$Proportion <- cumsum(case_counts$MA_sequenced) / case_counts$Cases
  out_df <- rbind(case_counts[,c("Date", "Cases")], data.frame(Date = case_counts$Date, Cases = cumsum(case_counts$MA_attempted)), data.frame(Date = case_counts$Date, Cases = cumsum(case_counts$MA_sequenced)), data.frame(Date = case_counts$Date, Cases = case_counts$Proportion))
  out_df$Site <- c(rep("Reported", nrow(case_counts)), rep("Attempted", nrow(case_counts)), rep("Sequenced", nrow(case_counts)), rep("Proportion", nrow(case_counts)))
  as_tibble(out_df)
}
  
make_lineage_group <- function(lin){
  if(substr(lin, 1, 1) == "A"){
    return("A")
  }
  if(lin == "B.1"){
    return("B.1")
  }else{
    return("B.other")
  }
}

# convert to state abbreviation
state_to_abbrev <- function(x){
  state.abb[match(x,state.name)]
}


# calculate MA imports from nextstrain tree


calc_MAimports <- function(a_df, conf_level, phy_tree){
  import_df <- data.frame(parent = rep(NA, 0), node = rep(NA, 0), MA = rep(NA, 0), numdate = rep(NA, 0))
  for(i in 1:nrow(a_df)){
    if(a_df$MA[i] == "YES" & a_df$MA_confYES[i] > conf_level & a_df$MA[a_df$parent[i]] == "NO" & a_df$MA_confNO[a_df$parent[i]] > conf_level){
      import_df <- rbind(import_df, a_df[i,c("parent", "node", "MA", "numdate")])
    }
  }
  import_df$clade_size <- sapply(import_df$node, function(x) length(grep("MA_", a_df$label[getDescendants(phy_tree, x)])))
  imp_nodes <- rep(0, 0)
  for(i in 1:nrow(import_df)){
    imp_nodes <- c(imp_nodes, getDescendants(phy_tree, import_df$node[i]))
  }
  import_df$reimport <- import_df$parent %in% imp_nodes
  import_df
}

# revised function
calc_MAimports2 <- function(a_df, conf_level, phy_tree){
  import_df <- data.frame(parent = rep(NA, 0), node = rep(NA, 0), MA = rep(NA, 0), numdate = rep(NA, 0))
  for(i in 1:nrow(a_df)){
    if(a_df$MA_confYES[i] > conf_level & a_df$MA_confYES[a_df$parent[i]] < (1-conf_level) & a_df$numdate[i] < 2020.363){
      import_df <- rbind(import_df, a_df[i,c("parent", "node", "MA", "numdate")])
    }
  }
  # count only imports in this study
  nodes_in_this_study <- a_df$node[grep("MA_", a_df$label)]
  nodes_in_this_study <- c(nodes_in_this_study, a_df$node[grep("MA1", a_df$label)])
  import_df$clade_size <- sapply(import_df$node, function(x) length(grep("MA_", a_df$label[getDescendants(phy_tree, x)])) + length(grep("MA-[QU]", a_df$label[getDescendants(phy_tree, x)])) + length(grep("MA1", a_df$label[getDescendants(phy_tree, x)])))
  imp_nodes <- rep(0, 0)
  for(i in 1:nrow(import_df)){
    imp_nodes <- c(imp_nodes, getDescendants(phy_tree, import_df$node[i]))
  }
  import_df$reimport <- import_df$parent %in% imp_nodes
  import_df
}


# calculate region imports from nextstrain tree
calc_Regionimports <- function(a_df, phy_tree){
  import_df <- data.frame(parent = rep(NA, 0), node = rep(NA, 0), MA = rep(NA, 0), parent_state = rep("", 0), branch.length = rep(0, 0))
  for(i in 1:nrow(a_df)){
    if(a_df$geoloc_ma[i] == "Massachusetts" & a_df$geoloc_ma[a_df$parent[i]] != "Massachusetts" & a_df$numdate[i] < 2020.363){
      p_state <- a_df$geoloc_ma[a_df$parent[i]]
      # make vector of confidence, find element-wise maximum
      p_conf_vec <- c(a_df$EuropeConf[a_df$parent[i]], a_df$AsiaConf[a_df$parent[i]], a_df$AfricaConf[a_df$parent[i]], a_df$NAConf[a_df$parent[i]], a_df$OceaniaConf[a_df$parent[i]], a_df$SAConf[a_df$parent[i]])
      child_conf <- a_df$MAConf[i]
      if((sum(p_conf_vec, na.rm=TRUE) > 0.9) & (child_conf > 0.9)){
        import_df <- rbind(import_df, data.frame(a_df[i,c("parent", "node", "MA", "branch.length")], parent_state = p_state))
      }
    }
  }
  import_df$clade_size <- sapply(import_df$node, function(x) length(grep("MA_", a_df$label[getDescendants(phy_tree, x)])) + length(grep("MA-[QU]", a_df$label[getDescendants(phy_tree, x)])) + length(grep("MA1", a_df$label[getDescendants(phy_tree, x)])))
  imp_nodes <- rep(0, 0)
  for(i in 1:nrow(import_df)){
    imp_nodes <- c(imp_nodes, getDescendants(phy_tree, import_df$node[i]))
  }
  import_df$reimport <- import_df$parent %in% imp_nodes
  import_df
}

# calculate MA imports from beast tree

count_importsMAbeast <- function(beast_tree){
  as_df <- as_tibble(beast_tree)
  imports <- data.frame(parent = rep(NA, 0), node = rep(NA, 0), MA = rep(NA, 0), height = rep(0, 0), label = rep(NA, 0))
  for(i in 1:nrow(as_df)){
    if(as_df$MA[i] == "YES" & as_df$MA[as_df$parent[i]] == "NO"){
      imports <- rbind(imports, as_df[i,c("parent", "node", "MA", "height", "label")])
    }
  }
  imports$clade_size <- sapply(imports$node, function(x) sum(as_df$MA[getDescendants(beast_tree@phylo, x)] == "YES" & !is.na(as_df$label[getDescendants(beast_tree@phylo, x)])))
  imports$height <- as.numeric(imports$height)
  #imp_nodes <- rep(0, 0)
  #for(i in 1:nrow(as_df)){
  #  imp_nodes <- c(imp_nodes, getDescendants(beast_tree@phylo, as_df$node[i]))
  #}
  #imports$reimport <- as_df$parent %in% imp_nodes
  imports
}

# calculate imports from parsimony-based reconstruction

calc_MAimports_parsimony <- function(a_df, phy_tree){
  # remove root, which has NA for parsimony reconstruction
  a_df <- a_df[!is.na(a_df$m_pars),]
  import_df <- data.frame(parent = rep(NA, 0), node = rep(NA, 0), MA = rep(NA, 0))
  for(i in 1:nrow(a_df)){
    if(a_df$m_pars[i] == 1 & a_df$m_pars[a_df$parent[i]] == 0){
      import_df <- rbind(import_df, a_df[i,c("parent", "node", "m_pars")])
    }
  }
  import_df$clade_size <- sapply(import_df$node, function(x) length(grep("MA_", a_df$label[getDescendants(phy_tree, x)])) + length(grep("MA-[QU]", a_df$label[getDescendants(phy_tree, x)])))
  imp_nodes <- rep(0, 0)
  for(i in 1:nrow(import_df)){
    imp_nodes <- c(imp_nodes, getDescendants(phy_tree, import_df$node[i]))
  }
  import_df$reimport <- import_df$parent %in% imp_nodes
  import_df
}

# relabel tips of nwk tree for BEAST analysis
relabel_tips <- function(tree, treedb){
  for(i in 1:length(tree$tip.label)){
    tree$tip.label[i] <- treedb$label[treedb$sample_id == tree$tip.label[i]]
  }
  tree
}

# convenience functions for simulating case counts
simRnormProb <- function(number_sims, sim_mean, sim_sd){
  probs <- rnorm(number_sims, mean = sim_mean, sd = sim_sd)
  probs[probs > 1] <- 1
  probs[probs < 0] <- 0
  probs
}

simRnormProportion <- function(number_sims, sim_mean, sim_sd){
  probs <- rnorm(n = number_sims, mean = sim_mean, sd = sim_sd)
  probs[probs > 1] <- 1
  probs
}

transformProb <- function(prob1){
  multiplier <- 1 + rgamma(1, shape = 0.5 , rate = 1)
  if(runif(1, min = 0 , max = 1)> 0.5){
    prob1 <- 1/(1+exp(-log((prob1*multiplier)/(1-prob1))))
  }else{
    multiplier <- 1/multiplier
    prob1 <- 1/(1+exp(-log((prob1*multiplier)/(1-prob1))))
  }
  prob1  
}

# simulate time-dependent allele frequencies

simTimeDep2416 <- function(n_sims, var_by_state_df){
    # set up data structures to store simulation results
    county_sims <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
    napprox <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
    napprox_epoch1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
    napprox_epoch2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
    
    bb <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
    bb_epoch1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
    bb_epoch2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  
    # loop through states  
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
      
      # beta-binomial model
      alpha.epoch1 <- variants_by_state$T0[i] + 0.5
      beta.epoch1 <- variants_by_state$C0[i] + 0.5
      alpha.epoch2 <- variants_by_state$T1[i] + 0.5
      beta.epoch2 <- variants_by_state$C1[i] + 0.5
      alpha.pooled <- variants_by_state$T0[i] + variants_by_state$T1[i] + 0.5
      beta.pooled <- variants_by_state$C0[i] + variants_by_state$C1[i] + 0.5
      
      # replace low count epochs with pooled quantity
      if(variants_by_state$alleles_epoch1[i] < 10){
        mu_hat.epoch1 <- mu_hat.pooled
        sd_hat.epoch1 <- sd_hat.pooled
        alpha.epoch1 <- alpha.pooled
        beta.epoch1 <- beta.pooled
      }
      # replace low count epochs with pooled quantity
      if(variants_by_state$alleles_epoch2[i] < 10){
        mu_hat.epoch2 <- mu_hat.pooled
        sd_hat.epoch2 <- sd_hat.pooled
        alpha.epoch2 <- alpha.pooled
        beta.epoch2 <- beta.pooled
      }
      
      # conduct the simulation
      napprox_epoch1[,i] <- rtnorm(n_sims, mu_hat.epoch1,sd_hat.epoch1, lower = 0, upper = 1) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch1[i]
      napprox_epoch2[,i] <- rtnorm(n_sims, mu_hat.epoch2,sd_hat.epoch2, lower = 0, upper = 1) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch2[i]
      napprox[,i] <- napprox_epoch1[,i] + napprox_epoch2[,i]
      bb_epoch1[,i] <- rbeta(n_sims, alpha.epoch1,beta.epoch1) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch1[i]
      bb_epoch2[,i] <- rbeta(n_sims, alpha.epoch2,beta.epoch2) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch2[i]
      bb[,i] <- bb_epoch1[,i] + bb_epoch2[,i]
  }
  result_list <- list(na1 = napprox_epoch1, na2 = napprox_epoch2, na0 = napprox, bb1 = bb_epoch1, bb2 = bb_epoch2, bb0 = bb)
  result_list
}

simTimeDep26233 <- function(n_sims, var_by_state_df){
  # set up data structures to store simulation results
  county_sims <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox_epoch1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox_epoch2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  
  bb <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  bb_epoch1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  bb_epoch2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  
  # loop through states  
  for(i in 1:nrow(variants_by_state)){
    mu_hat.epoch1 <- variants_by_state$AF0[i]
    mu_hat.h.epoch1 <- (variants_by_state$T0[i] + 0.5) / (variants_by_state$alleles_epoch1[i] + 1) #continuity correction
    sd_hat.epoch1 <- sqrt(mu_hat.h.epoch1*(1-mu_hat.h.epoch1)/(variants_by_state$alleles_epoch1[i]+1))
    mu_hat.epoch2 <- variants_by_state$AF1[i]
    mu_hat.h.epoch2 <- (variants_by_state$T1[i] + 0.5) / (variants_by_state$alleles_epoch2[i] + 1) #continuity correction
    sd_hat.epoch2 <- sqrt(mu_hat.h.epoch2*(1-mu_hat.h.epoch2)/(variants_by_state$alleles_epoch2[i]+1))
    mu_hat.pooled <- (variants_by_state$T0[i] + variants_by_state$T1[i]) / (variants_by_state$T0[i] + variants_by_state$T1[i] + variants_by_state$G0[i] + variants_by_state$G1[i])
    mu_hat.h.pooled <- (variants_by_state$T0[i] + variants_by_state$T1[i] + 0.5) / (variants_by_state$T0[i] + variants_by_state$T1[i] + variants_by_state$G0[i] + variants_by_state$G1[i] + 1)
    sd_hat.pooled <- sqrt(mu_hat.h.pooled*(1 - mu_hat.h.pooled)/(variants_by_state$alleles_epoch1[i]+ variants_by_state$alleles_epoch2[i]+1))
    
    # beta-binomial model
    alpha.epoch1 <- variants_by_state$T0[i] + 0.5
    beta.epoch1 <- variants_by_state$G0[i] + 0.5
    alpha.epoch2 <- variants_by_state$T1[i] + 0.5
    beta.epoch2 <- variants_by_state$G1[i] + 0.5
    alpha.pooled <- variants_by_state$T0[i] + variants_by_state$T1[i] + 0.5
    beta.pooled <- variants_by_state$G0[i] + variants_by_state$G1[i] + 0.5
    
    # replace low count epochs with pooled quantity
    if(variants_by_state$alleles_epoch1[i] < 10){
      mu_hat.epoch1 <- mu_hat.pooled
      sd_hat.epoch1 <- sd_hat.pooled
      alpha.epoch1 <- alpha.pooled
      beta.epoch1 <- beta.pooled
    }
    # replace low count epochs with pooled quantity
    if(variants_by_state$alleles_epoch2[i] < 10){
      mu_hat.epoch2 <- mu_hat.pooled
      sd_hat.epoch2 <- sd_hat.pooled
      alpha.epoch2 <- alpha.pooled
      beta.epoch2 <- beta.pooled
    }
    
    # conduct the simulation
    napprox_epoch1[,i] <- rtnorm(n_sims, mu_hat.epoch1,sd_hat.epoch1, lower = 0, upper = 1) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch1[i]
    napprox_epoch2[,i] <- rtnorm(n_sims, mu_hat.epoch2,sd_hat.epoch2, lower = 0, upper = 1) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch2[i]
    napprox[,i] <- napprox_epoch1[,i] + napprox_epoch2[,i]
    bb_epoch1[,i] <- rbeta(n_sims, alpha.epoch1,beta.epoch1) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch1[i]
    bb_epoch2[,i] <- rbeta(n_sims, alpha.epoch2,beta.epoch2) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases_epoch2[i]
    bb[,i] <- bb_epoch1[,i] + bb_epoch2[,i]
  }
  result_list <- list(na1 = napprox_epoch1, na2 = napprox_epoch2, na0 = napprox, bb1 = bb_epoch1, bb2 = bb_epoch2, bb0 = bb)
  result_list
}

# 2416 non time-dependent simulation 

sim2416 <- function(n_sims, variants_by_state){
  # set up data structures to store results of simulation
  napprox <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox_robust1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox_robust2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  bb <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  
  # conduct simulation
  for(i in 1:nrow(variants_by_state)){
    # beta-binomial model
    alpha.epoch <- variants_by_state$T[i] + 0.5
    beta.epoch <- variants_by_state$C[i] + 0.5
    bb[,i] <- rbeta(n_sims, alpha.epoch,beta.epoch) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases[i]
    
    # normal-approximation model
    mu_hat <- variants_by_state$AF[i]
    mu_hat.h <- (variants_by_state$T[i] + 0.5) / (variants_by_state$alleles[i] + 1)
    sd_hat <- sqrt((mu_hat.h*(1-mu_hat.h))/(variants_by_state$alleles[i] + 1))
    napprox[,i] <- rtnorm(n_sims, mu_hat,sd_hat, lower = 0, upper = 1)  * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases[i]
    napprox_robust1[,i] <- rtnorm(n_sims, mu_hat,2*sd_hat, lower = 0, upper = 1) *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases[i]
    napprox_robust2[,i] <- rtnorm(n_sims, mu_hat,3*sd_hat, lower = 0, upper = 1) *rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases[i]
  }
  result_list <- list(na0 = napprox, nar1 = napprox_robust1, nar2 = napprox_robust2, bb0 = bb)
  result_list
}

# 26233 non time-dependent

sim26233 <- function(n_sims, variants_by_state){
  # set up data structures to store results of simulation
  napprox <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox_robust1 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  napprox_robust2 <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  bb <- matrix(rep(NA, n_sims*nrow(variants_by_state)), nrow = n_sims)
  
  for(i in 1:nrow(variants_by_state)){
    # beta-binomial model
    alpha.epoch <- variants_by_state$T[i] + 0.5
    beta.epoch <- variants_by_state$G[i] + 0.5
    bb[,i] <- rbeta(n_sims, alpha.epoch,beta.epoch) * rtnorm(n_sims, 0.9, 0.05, lower = 0, upper = 1) * variants_by_state$cases[i]
    
    # normal-approximation model
    mu_hat <- variants_by_state$AF[i]
    mu_hat.h <- (variants_by_state$T[i] + 0.5) / (variants_by_state$alleles[i] + 1)
    sd_hat <- sqrt((mu_hat.h*(1-mu_hat.h))/(variants_by_state$alleles[i] + 1))
    napprox[,i] <- rtnorm(n_sims, mu_hat,sd_hat, lower = 0, upper = 1) * variants_by_state$cases[i]
    napprox_robust1[,i] <- rtnorm(n_sims, mu_hat,2*sd_hat, lower = 0, upper = 1) * variants_by_state$cases[i]
    napprox_robust2[,i] <- rtnorm(n_sims, mu_hat,3*sd_hat, lower = 0, upper = 1) * variants_by_state$cases[i]
  }
  result_list <- list(na0 = napprox, nar1 = napprox_robust1, nar2 = napprox_robust2, bb0 = bb)
  result_list
}

