# packages #
#install.packages("adegenet")
#install.packages("haplotypes")
library(adegenet)
library(haplotypes)

# metadata #
group_info<-read.csv("Trees/metadata.tsv",sep="\t")
interest_group<-"BHCHP"
#interest_group<-"SNFA"

# test data #
test_seq<-read.fas("Trees/trimmed_alignment.fasta")
test_seq_df<-as.list(test_seq)
samples<-sapply(names(test_seq_df),function(x) strsplit(x,"[|]")[[1]][1])
of_interest<-group_info$sample_id[group_info$exposure==interest_group]
for_subset<-which(samples %in% of_interest)
test_seq_subset_df<-test_seq_df[for_subset]
test_seq_subset<-as.dna(test_seq_subset_df)

# distance graphs of genomes # (using gengraph from the R package adegenet)
d<-distance(test_seq_subset,indels="missing")
x<-as.matrix(d)
x[x==0]<-1
g<-gengraph(d+1,cutoff=5)
plot(g$graph,vertex.label=NA)
