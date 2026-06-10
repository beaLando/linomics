# LIBRARIES
library(ggplot2)
library(tidyverse)
# library(devtools)
# install_github("nbedelman/quiblR")
library(quiblR)
library("ape")
library("hash")
library("ggtree")
library("ggpubr")


# QUIBL
largeSpeciesTree <- read_speciesTree("./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.phyxed.tre")
largeQuiblOutput <- read_csv_quibl("./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/quiblOut.csv")
largeOriginalTrees <- read.tree("./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned190_prunedSPP_completeCases.nwk")

plot(largeSpeciesTree)
largeQuiblOutput

totalTrees <- length(largeOriginalTrees)

largeSpeciesTree$tip.label <- gsub("_nuclear", "", largeSpeciesTree$tip.label)
largeQuiblOutput$triplet <- gsub("_nuclear", "", largeQuiblOutput$triplet)
largeQuiblOutput$outgroup <- gsub("_nuclear", "", largeQuiblOutput$outgroup)

## Assign species topology & Mark which of the topologies' 2-distribution model is a significantly better fit than the ILS only model, and calculate the proportion of loci with a history of introgression for the full data
largeQuiblOutput <- mutate(largeQuiblOutput,
                           isDiscordant = as.integer(! apply(largeQuiblOutput, 1, isSpeciesTree, sTree=largeSpeciesTree)),
                           isSignificant = as.integer(apply(largeQuiblOutput, 1, testSignificance, threshold=10)),
                           totalIntrogressionFraction=(mixprop2*count*isDiscordant)/totalTrees,
                           significantFraction = (isSignificant * count * isDiscordant) / totalTrees)

largeQuiblOutput %>%
  filter(isDiscordant == 1 & isSignificant == 1) %>% View()

largeQuiblOutput <- largeQuiblOutput %>%
  group_by(rownames(.)) %>%
  mutate(tax = gsub(unique(as.character(outgroup)), "", triplet)) %>%
  mutate(tax = gsub("__", "_", tax)) %>%
  mutate(tax = gsub("^_", "", tax)) %>%
  mutate(tax = gsub("_$", "", tax)) %>%
  ungroup() %>%
  separate(tax, sep = "_", into = c("tax1", "tax2")) %>%
  select(-c('rownames(.)')) %>%
  as.data.frame() #%>% str()

write.csv(largeQuiblOutput, "./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/quiblOut_formatted.csv")

## Produce introgression summary & make heatmap
introgressionSummary <- getIntrogressionSummary(largeQuiblOutput,largeSpeciesTree)
#introgressionSummary$value <- as.numeric(ifelse(introgressionSummary$tax1==introgressionSummary$tax2,NA_integer_,introgressionSummary$value))
head(introgressionSummary)

introgressionSummary.BIC <- largeQuiblOutput %>%
  dplyr::filter(isDiscordant == 1) %>%
  droplevels() %>%
  select(tax1, tax2, BIC1Dist, BIC2Dist, count) %>%
  mutate(deltaBIC = BIC1Dist - BIC2Dist) %>%
  group_by(tax1, tax2) %>%
  summarize(
    deltaBIC_wmean     = weighted.mean(deltaBIC, w = count, na.rm = TRUE),
    prop_significant   = sum(count[deltaBIC >= 10], na.rm = TRUE) / sum(count),
    prop_ILS           = sum(count[deltaBIC <= -10],  na.rm = TRUE) / sum(count),
    n_trees            = sum(count)
  ) %>%
  ungroup() %>%
  as.data.frame()

introgressionSummary <- introgressionSummary %>%
  left_join(.,
            introgressionSummary.BIC,
            by = c("tax1", "tax2")) %>%
  as.data.frame()

speciesTreeSubset <- ggtree(extractTripletTree(largeSpeciesTree, unique(introgressionSummary$tax1)))
speciesTreeSubset_down <- ggtree(extractTripletTree(largeSpeciesTree, unique(introgressionSummary$tax1)))+ coord_flip()

speciesTreeSubset + geom_tiplab()
speciesTreeSubset_down + geom_tiplab()


introgressionSummary <- bind_rows(
  introgressionSummary,
  introgressionSummary %>% rename(tax1 = tax2, tax2 = tax1)
) %>%
  distinct()

species <- c("KCPT", "TXMP", "VXOD", "HNCF", "ERR2040380", "AEPI", "BHYC", "D07", "C96", "D03", "D13", "D01", "D02")

diagonal <- data.frame(
  tax1 = species,
  tax2 = species,
  deltaBIC = NA_real_,
  value = NA_real_
)

introgressionSummary <- bind_rows(introgressionSummary, diagonal) %>% #str()
  mutate(
    tax1 = factor(tax1, levels = species),
    tax2 = factor(tax2, levels = rev(species))
  ) #%>% str()


introgressionSummary$prop_significant <- as.character(if_else(round(introgressionSummary$prop_significant, digits = 2)==0, NA, round(introgressionSummary$prop_significant, digits = 2)))
introgressionSummary$deltaBIC_wmean <- if_else(introgressionSummary$deltaBIC_wmean < 10 | is.na(introgressionSummary$deltaBIC_wmean), "", "*")
introgressionSummary$label <- paste0(introgressionSummary$prop_significant, introgressionSummary$deltaBIC_wmean)
introgressionSummary$label <- if_else(introgressionSummary$label == "NA", "", introgressionSummary$label)

summaryMatrix <- ggplot(data = introgressionSummary, aes(tax1, tax2, fill = value))+
  geom_tile(col = "white")+
  #geom_text(aes(label = if_else(is.na(deltaBIC), "", "sig."))) +
  geom_text(aes(label = label)) +
  scale_fill_gradient2(low = "white", mid = "goldenrod", high = "red", na.value = "white",
                       midpoint = 0.07, limit = c(0,max(introgressionSummary$value, na.rm = TRUE)),
                       name="Average introgression fraction") +
  geom_abline(slope = -1, intercept=0)+
  # geom_vline(xintercept=seq(1.5,nrow(introgressionSummary)+0.5,1), alpha=0.6)+
  # geom_hline(yintercept=seq(1.5,nrow(introgressionSummary)+0.5,1), alpha=0.6)+
  labs(x="", y="")+
  #scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "left")

ggarrange(  speciesTreeSubset, summaryMatrix, NULL, speciesTreeSubset_down,
            ncol = 2,nrow=2, heights=c(2,1), widths=c(1,2), align="hv",common.legend = TRUE)
