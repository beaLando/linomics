library(tidyverse)

getwd()

#parts <- read.table("./data/Angiosperm353/aln/aligned_trimmed_clean_supermatrix/supermx353_partitions.txt", sep = ",", header = FALSE)
#parts <-read.table("./data/Angiosperm353/aln_noPaftol/aligned_trimmed_clean_supermatrix/supermx353_partitions.txt", sep = ",", header = FALSE)
parts <-read.table("./data/Angiosperm353/aln_noPerenne/aligned_trimmed_clean_supermatrix/supermx353_partitions.txt", sep = ",", header = FALSE)
str(parts)

parts1 <- parts %>%
  select(-V1) %>%
  mutate(grp = rep(c("1","2"), nrow(.)/2)) %>%
  filter(grp == "1") %>%
  droplevels() %>%
  bind_cols(.,
            parts %>%
              select(-V1) %>%
              mutate(grp = rep(c("1","2"), nrow(.)/2)) %>%
              filter(grp == "2") %>%
              droplevels()) %>% #str()
  mutate(position = str_replace(V2...3, "s", ""),
         geneName_geneOrder = str_replace(str_replace(str_replace(V2...1, "s", ""), ".mft.tal", ""), "gene", "part")) %>% #str()
  select(geneName_geneOrder, position) %>%
  separate(geneName_geneOrder, sep = "_", into = c("geneName", "geneOrder"), remove = TRUE) %>% #str()
  mutate(type = rep("DNA", nrow(.))) %>%
  select(type, geneOrder, geneName, position) %>% #View()
  mutate(partition = paste0(type, ", ", geneOrder, " = ", position),
         partition_beauti = paste0("charset ", geneOrder, " = ", position, ";")) %>% #str()
  as.data.frame()

#check for errors from https://github.com/Cibiv/IQ-TREE/issues/208:

#1 "Too large site ID" if (upper >= aln->getNSite()) throw "The following line of your partition definitions includes a site with an ID that is larger than your alignment. All site ID's in your partition file must refer to a column number of the alignment between 1 and the length of the alignment. Please fix this and try again.";
#2 "Negative site ID" if (lower < 0) throw "The following line of your partition definitions includes a site with a negative ID, which is impossible. All site ID's in your partition file must refer to a column number of the alignment between 1 and the length of the alignment. Please fix this and try again.";
#3 "Wrong range" if (lower > upper) throw "The following line of your partition definitions includes a range of sites where the lower end of the range (the number before dash) is higher than the upper end of the range (the number after the dash). Please fix this and try again.";
#4 "Wrong step size" if (step < 1) throw "The following line of your partition definitions includes a range of sites where the step size (the number after the \) is less than one. Step sizes have to be at least 1. Please fix this and try again.";

parts1 %>% #View()
  separate(position, sep = "-", into = c("position_lower", "position_upper")) %>%
  mutate(position_lower = as.numeric(as.character(position_lower)),
         position_upper = as.numeric(as.character(position_upper))) %>%
  mutate(test3 = if_else(position_lower > position_upper, "fail", "pass")) %>%
  View()



parts2 <- parts1 %>% select(partition) %>% as.data.frame()
parts3 <- parts1 %>% select(partition_beauti) %>% as.data.frame()  


#write.table(parts1, "./data/Angiosperm353/aln/aligned_trimmed_clean_supermatrix/supermx353_partitions_R.txt")
#write.table(parts2, "./data/Angiosperm353/aln/aligned_trimmed_clean_supermatrix/supermx353_partitions_final.txt", row.names = FALSE, col.names = FALSE)
#write.table(parts1, "./data/Angiosperm353/aln_noPaftol/aligned_trimmed_clean_supermatrix/supermx353_partitions_R.txt")
#write.table(parts2, "./data/Angiosperm353/aln_noPaftol/aligned_trimmed_clean_supermatrix/supermx353_partitions_final.txt", row.names = FALSE, col.names = FALSE)
write.table(parts1, "./data/Angiosperm353/aln_noPerenne/aligned_trimmed_clean_supermatrix/supermx353_partitions_R.txt")
write.table(parts2, "./data/Angiosperm353/aln_noPerenne/aligned_trimmed_clean_supermatrix/supermx353_partitions_final.txt", row.names = FALSE, col.names = FALSE)
write.table(parts3, "./output/Angiosperm353/beast/supermx353_partitions2add.txt", row.names = FALSE, col.names = FALSE)
