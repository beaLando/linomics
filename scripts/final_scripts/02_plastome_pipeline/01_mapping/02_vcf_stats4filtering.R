#see https://speciationgenomics.github.io/filtering_vcfs/


library(tidyverse)

# UNFILTERED VCF
var_qual <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

# var_qual <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.lqual", delim = "\t",
#                        col_names = c("chr", "pos", "qual"), skip = 1)

str(var_qual)

ggplot(var_qual, aes(qual)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() +
  xlim(c(0,500)) #>30



var_depth <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# var_depth <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.ldepth.mean", delim = "\t",
#                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

str(var_depth)
summary(var_depth$mean_depth)

ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() #min 100 max 177x2=354



var_miss <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

str(var_miss)
summary(var_miss$fmiss)

ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() #no need to filter, for each site, genotypes are absent in a maximum of 3% individuals



var_freq <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

summary(var_freq$maf) #not sure how this translates to plastome?



ind_depth <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

# ind_depth <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.idepth", delim = "\t",
#                         col_names = c("ind", "nsites", "depth"), skip = 1)

str(ind_depth)

ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()




ind_miss  <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)


str(ind_miss)

ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light() #very very small proportion of missing data for most inds, leave as is

ind_miss$ind[which.max(ind_miss$fmiss)] #maybe I could remove L24



# VCF post-FILTERING
var_qual <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

str(var_qual)

ggplot(var_qual, aes(qual)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() #>30



var_depth <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

str(var_depth)
summary(var_depth$mean_depth)

ggplot(var_depth, aes(mean_depth)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() #min 100 max 177x2=354



var_miss <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

str(var_miss)
summary(var_miss$fmiss)

ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light() #no need to filter, for each site, genotypes are absent in a maximum of 3% individuals



var_freq <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

ggplot(var_freq, aes(maf)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

summary(var_freq$maf) #not sure how this translates to plastome?



ind_depth <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

str(ind_depth)

ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()



ind_miss  <- read_delim("./data/Lb_plastomes_raw/Consensus/stats_vcf_unfiltered/all_plastomes_filtered.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

str(ind_miss)

ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light() #very very small proportion of missing data for most inds, leave as is

ind_miss$ind[which.max(ind_miss$fmiss)] #maybe I could remove L24

