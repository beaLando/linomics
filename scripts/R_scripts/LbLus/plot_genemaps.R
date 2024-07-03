library(genoPlotR)

# Annotated plastomes for the 4 lineages (I picked samples with concordance between nuclear and chloroplast genomes)
MedSW_p1 <- read_dna_seg_from_genbank("./output/plastome/annotation_denovo/annotation_denovo_L40_GenBank.gb")
AtlNW_pVil <- read_dna_seg_from_genbank("./output/plastome/annotation_denovo/annotation_denovo_L49_GenBank.gb") #lebanon sample
AtlNW_pSUT <- read_dna_seg_from_genbank("./output/plastome/annotation_denovo/annotation_denovo_L41_GenBank.gb")

# Transform to dna_seq objects
MedSW.seg <- dna_seg(MedSW_p1)
AtlNW1.seg <- dna_seg(AtlNW_pVil)
AtlNW2.seg <- dna_seg(AtlNW_pSUT)

list.segs <- list(MedSW.seg, AtlNW1.seg) #, AtlNW2.seg
seg.nms <- c("MedSW", "AtlNW1") #, "AtlNW2"
names(list.segs) <- seg.nms

annots <- lapply(list.segs, function(x){
  mid <- middle(x)
  annot <- annotation(x1=mid, text=x$name, rot=30)
  idx <- grep("^[^B]", annot$text, perl=TRUE)
  annot[idx,] #annot[idx[idx %% 4 == 0],]
  })


# Transform to comparison object

SW_NW1.df <- data.frame(start1 = MedSW.seg$start,
                        end1 = MedSW.seg$end,
                        start2 = AtlNW1.seg$start,
                        end2 = AtlNW1.seg$end)
SW_NW1.com <-comparison(SW_NW1.df)

# SW_NW2.df <- data.frame(start1 = MedSW.seg$start,
#                         end1 = MedSW.seg$end,
#                         start2 = AtlNW2.seg$start,
#                         end2 = AtlNW2.seg$end)
# SW_NW2.com <-comparison(SW_NW2.df)


list.coms <- list(SW_NW1.com)



# Plot

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              #annotations = annots,
              legend = TRUE)

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(0, 20000), c(0, 20000))) #, c(0, 20000)

plot_gene_map(dna_segs = list.segs,
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(20000, 40000), c(20000, 40000)))

###plot_gene_map(dna_segs = list.segs,
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(40000, 60000), c(40000, 60000)))

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(60000, 80000), c(60000, 80000)))

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(80000, 100000), c(80000, 100000)))

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(100000, 120000), c(100000, 120000)))

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(120000, 140000), c(120000, 140000)))

plot_gene_map(dna_segs = list.segs, 
              comparisons = list.coms,
              annotations = annots,
              legend = TRUE,
              xlims = list(c(140000, 157000), c(140000, 157000)))
