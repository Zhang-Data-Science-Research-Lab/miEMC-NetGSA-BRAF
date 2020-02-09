
# Henry Linder 2020 mhlinder@gmail.com

# * Loads all raw miRNA files downloaded using TCGA-Assembler
# * Reads each raw file and records all miRNA IDs
# * Saves a text file with rows containing miRNA IDs

datadir <- "~/local/data/TCGA/mirna"

library(magrittr)
library(dplyr)
library(readr)

all.fn <- list.files(datadir, pattern="BRCA", full.names=TRUE)

all.mirna <- c()
for (d in all.fn) {
    ## Get TCGA cancer files
    fn <- list.files(d, full.names=TRUE)
    ## Only consider a specific reference genome
    fn <- fn[grepl("hg19", fn) &
             grepl("mirbase20", fn) &
             grepl("HiSeq", fn)]
    if (length(fn) == 0) {
        next
    }
    ## Unique miRNAs
    all.mirna <- c(all.mirna, read_tsv(fn) %>% use_series(miRNA_ID)) %>% unique %>% sort
}
all.mirna <- all.mirna[all.mirna != "miRNA_ID"]

data.frame(miRNA=all.mirna) %>%
    write_csv("all-mirna.csv")
