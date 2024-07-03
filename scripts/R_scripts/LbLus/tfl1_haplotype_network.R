library(tidyverse)
library(ggplot2)
library(phangorn)

#Import population info
pop.info <- read.csv("./data/lblus_all_info.csv")

aln <- ape::read.FASTA("./data/Lb_TFL1/lb_all.mft.tal.fasta")

labs.fixed <- data.frame(inds = str_remove(labels(aln), " 1281 bp")) %>%
  separate(inds, into = c("inds", "inds2"), sep = "L") %>% #View()
  mutate(
    inds = if_else(inds == "", "L", inds),
    inds2 = sprintf("%02d", as.numeric(inds2))
  ) %>%
  mutate(inds.fixed = paste0(inds, inds2)) %>%
  mutate(inds.fixed = str_remove(inds.fixed, "NA")) %>%
  as.data.frame()

names(aln) <- labs.fixed$inds.fixed

as.data.frame(labels(aln)) %>%
  filter(!`labels(aln)` %in% pop.info$seq_code)

pop.info2 <- pop.info %>%
  select(-X) %>% 
  filter(seq_code %in% labels(aln)) %>% #str()
  mutate(purpose = case_when(
    purpose == "wild" ~ "wild",
    purpose == "oil" ~ "oil",
    purpose == "fibre" ~ "fiber",
    purpose == "fiber" ~ "fiber",
    purpose == "unknown" ~ "unkn-interm",
    purpose == "intermediate" ~ "unkn-interm"
  )) %>%
  mutate(
    region_purpose = if_else(spp == "usitatissimum", purpose, region),
    LuTFL1 = if_else(is.na(LuTFL1) == TRUE, "new", LuTFL1)) %>%
  mutate(LuTFL1 = if_else(LuTFL1 == "", "others", LuTFL1)) %>%
  distinct(seq_code, .keep_all = TRUE) %>% #str()
  arrange(factor(seq_code, levels = labels(aln))) %>% #str()
  as.data.frame() #%>% View()

# network to retrieve individuals
h <- pegas::haplotype(aln, strict = FALSE, trailingGapsAsN = TRUE)

hname <- paste("H", 1:nrow(h), sep = "")
rownames(h) <- paste(hname)
net <- pegas::haploNet(h, d = NULL, getProb = TRUE)

ind.hap <- with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=labels(aln))
)

ind.hap <- as.data.frame(ind.hap) %>%
  filter(Freq == 1) %>%
  select(-Freq) %>% #View()
  left_join(.,
            pop.info2,
            by = c("individuals" = "seq_code")) %>%
  rename(seq_code = "individuals",
         TFL1_new = "hap") %>%
  select(seq_code:alt, region_purpose, LuTFL1, TFL1_new, LuTFL2:RADseq) %>%
  as.data.frame()

write.csv(ind.hap, "./output/tfl1/lblus_all_info&haplo.csv")

# network colored by gutaker haplotype
aln.grp <- aln
names(aln.grp) <- pop.info2$LuTFL1

h.grp <- pegas::haplotype(aln.grp, strict = FALSE, trailingGapsAsN = TRUE)

hname <- paste("H", 1:nrow(h.grp), sep = "")
rownames(h.grp) <- paste(hname)
net.grp <- pegas::haploNet(h.grp, d = NULL, getProb = TRUE)

ind.hapg <- with(
  utils::stack(setNames(attr(h.grp, "index"), rownames(h.grp))),
  table(hap=ind, individuals=labels(aln.grp))
)

dev.off()
plot(net.grp, size = attr(net.grp, "freq"), scale.ratio = 10, cex = 0.6, labels = TRUE, pie = ind.hapg, show.mutation = 0, font = 1, fast = TRUE)
legend(x = -100, y = -40, colnames(ind.hapg), fill = rainbow(ncol(ind.hapg)), cex=0.5, ncol=6, x.intersp=0,003, text.width=9)


# network colored by purpose_region
aln.grp <- aln
names(aln.grp) <- pop.info2$region_purpose

h.grp <- pegas::haplotype(aln.grp, strict = FALSE, trailingGapsAsN = TRUE)

hname <- paste("H", 1:nrow(h.grp), sep = "")
rownames(h.grp) <- paste(hname)
net.grp <- pegas::haploNet(h.grp, d = NULL, getProb = TRUE)

ind.hapg <- with(
  utils::stack(setNames(attr(h.grp, "index"), rownames(h.grp))),
  table(hap=ind, individuals=labels(aln.grp))
)

dev.off()
plot(net.grp, size = attr(net.grp, "freq"), scale.ratio = 10, cex = 0.6, labels = TRUE, pie = ind.hapg, show.mutation = 0, font = 1, fast = TRUE)
legend(x = -100, y = -40, colnames(ind.hapg), fill = rainbow(ncol(ind.hapg)), cex=0.5, ncol=6, x.intersp=0,003, text.width=9)
