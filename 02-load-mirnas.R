
## Henry Linder 2020 mhlinder@gmail.com

## * Loads mirDIP database
## * Constructs a lookup by miRNA ID of all target genes
## * Plots the scores of all miRNA for each gene in BRAF pathway
## * Subsets to only miRNA targeting BRAF pathway genes, with an 
##   integrated score of "very high" confidence of targeting the
##   genes
## * Plots the BRAF pathway, gene expression only
## * Constructs the adjacency matrix for the full BRAF pathway 
##   with EMC and miRNA
## * Plots the full BRAF pathway
## * Saves the adjacency matrix

lookup.path <- "~/local/data/He_Li_mirDIP_Unidirectional_search_v_4_1_1547654112/mirDIP_Unidirectional_search_v.4.1.txt"

library(magrittr)
library(dplyr)

library(readr)
library(igraph)
library(scales)

source("color.R")

## Lookup of gene symbol by miRNA
lookup.all <- read_tsv(lookup.path,
              col_names=c("GeneSymbol", "miRNA", "SourceNumber",
                          "IntegratedRank", "Sources", "ScoreClass")) %>%
    mutate(miRNA = tolower(miRNA))

## All genes, miRNAs in the lookup
genes.all <- lookup.all$GeneSymbol %>% unique %>% sort
mirna.all <- lookup.all$miRNA %>% unique %>% sort

g.all <- length(genes.all)
m.all <- length(mirna.all)

## miRNA in TCGA dataset
mirna.tcga <- read_tsv("all-mirna.csv")$miRNA
m.tcga <- length(mirna.tcga)

## Genes in BRAF pathway
genes.braf <- readLines("BRAF_E.txt") %>%
    tail(-1) %>%
    strsplit(" ") %>%
    unlist %>%
    unique %>%
    sort
g.braf <- length(genes.braf)

lookup <- lookup.all %>% filter(GeneSymbol %in% genes.braf)

## Lookup all genes for each miRNA
all.lookups <- list()
for (i in 1:m.tcga) {
    if (i %% 100 == 0) print(sprintf("%d of %d", i, m.tcga))
    x <- mirna.tcga[i]
    re.mirna <- sprintf("%s[\"-]", x)
    l <- lookup %>% filter(grepl(re.mirna, miRNA))
    if (nrow(l) > 0) {
        l$miRNA.TCGA <- x
        all.lookups[[i]] <- l
    }
}

## All matches of TCGA miRNAs with BRAF genes
l.all <- bind_rows(all.lookups) %>%
    mutate(ScoreClass = factor(ScoreClass, levels=c("V", "H", "M", "L"))) %>%
    arrange(GeneSymbol, miRNA)

cols_ranks <- cols_all[paste0(c("Blue", "GreenB", "Yellow", "Red"), 1)]

## png("all-ranks.png", width = 1600, height = 900)
setEPS()
postscript("all-ranks.eps", width = 16, height = 10)
par(cex.lab=1.5, cex.main=2, cex.axis=1.4)
plot(l.all$IntegratedRank, type="n",
     axes=FALSE, ylim=c(0, 1),
     xlab="miRNA, nested within gene",
     ylab="Integrative score",
     main="Integrative score predicting miRNA targeting of BRAF pathway genes")
usr <- par("usr")
abline(h=0, lwd=1.5)
## 
ix.gene <- c(which(!duplicated(l.all$GeneSymbol)), nrow(l.all))
abline(v=ix.gene, lwd=1.5, lty=2, col="darkgrey")
abline(h=seq(0, 1, by=.2)[-1], lwd=1.5, lty=2, col="darkgrey")
## Data points
points(l.all$IntegratedRank,
       col=cols_ranks[l.all$ScoreClass],
       pch=19, cex=1.2)
##
ix.lab <- rowMeans(cbind(head(ix.gene, -1),
                         tail(ix.gene, -1)))
axis(1, at=c(ix.gene[1], tail(ix.gene, 1)), labels=c("",""))
axis(1, at=ix.lab, labels=genes.braf)
##
axis(2, at=seq(0, 1, by=.1))
##
box()
legend("topleft", col=cols_ranks, pch=19,
       legend=paste(c("Very high", "High", "Medium", "Low"),
                    "confidence"),
       inset=c(.05, .05), bty="n", cex=1.5)
dev.off()

## Only matches with very high confidence
l <- l.all %>%
    filter(ScoreClass == "V")
mirna.V <- l$miRNA %>% unique %>% sort
m.V <- length(mirna.V)

write_csv(l, path="BRAF-mirna-gene-targets.csv")

A.V11 <- matrix(0, g.braf, g.braf, dimnames=list(genes.braf, genes.braf))
A.V12 <- matrix(0, g.braf, m.V, dimnames=list(genes.braf, mirna.V))

A.V21 <- matrix(0, m.V, g.braf, dimnames=list(mirna.V, genes.braf))
A.V22 <- matrix(0, m.V, m.V, dimnames=list(mirna.V, mirna.V))

A.V <- rbind(cbind(A.V11, A.V12),
             cbind(A.V21, A.V22))

for (g in genes.braf) {
    l.g <- l %>% filter(GeneSymbol == g)
    if (nrow(l.g) > 0) {
        A.V[g, l.g$miRNA] <- 1
    }
}

graph.V <- graph_from_adjacency_matrix(t(A.V))

set.seed(20190116)
lo <- layout_(graph.V, with_fr())

all.vert <- colnames(A.V)
names(all.vert) <- all.vert
all.vert[!all.vert %in% genes.braf] <- ""

all.edges <- as.data.frame(get.edgelist(graph.V), stringsAsFactors=FALSE)
colnames(all.edges) <- c("miRNA", "Gene")

degree.counts <- all.edges %>%
    dplyr::select(miRNA) %>%
    group_by(miRNA) %>%
    summarize(n=n()) %>%
    ungroup %>%
    rename(Shared=n)

all.edges <- all.edges %>% left_join(degree.counts, by="miRNA")

cols.nodes <- c("grey95", cols_all[c("Red2", "Yellow2", "GreenB2", "Blue2", "Purple2")])
cols.edges <- c("grey90",
                "grey70",
                "grey70",
                "grey30",
                "black")

mirna.counts <- (all.edges$Shared %>% table) / 1:5

## w <- 2000
## h <- 1600
## png("graph-braf-mirna.png", width = w, height = h)
w <- 24
h <- 20
setEPS()
postscript("graph-braf-mirna.eps", width = w, height = h,
           fonts="serif")
par(cex.main=3)
plot(graph.V, layout=lo,
     vertex.size=c(rep(12, g.braf), rep(2, m.V)),
     vertex.color=c(rep(cols.nodes[1], g.braf), cols.nodes[-1][degree.counts$Shared]),
     vertex.label=all.vert,
     vertex.label.color="black",
     vertex.label.cex=1.5,
     vertex.label.font=1,
     edge.width=c(3, 5, 5, 7, 10)[all.edges$Shared],
     edge.color=cols.edges[all.edges$Shared],
     edge.arrow.size=0, edge.arrow.width=0,
     asp=1/(w/h))
labs <- sprintf("miRNA targets %d genes (N=%d)", 2:5, mirna.counts[2:5])
labs <- c(matrix(c(labs, rep("", 2*length(labs))), byrow=TRUE, nrow=3))
legend("bottomright",
       legend=c("Gene",
                sprintf("miRNA targets 1 gene   (N=%d)", mirna.counts[1]),
                "",
                "",
                labs),
       lty=c(NA, rep(c(NA, 1, NA), 5)),
       lwd=c(NA, NA, 3, NA, NA, 5, NA, NA, 5, NA, NA, 7, NA, NA, 10),
       col=c("black", c(matrix(c(cols.nodes[-1], cols.edges, rep("", 5)), byrow=TRUE, nrow=3))),
       pch=c(21, rep(c(19, NA, NA), 5)),
       pt.bg=c(cols.nodes[1], rep(NA, 15)),
       cex=2.5, pt.cex=3,
       bty="n")
dev.off()

gene.edges <-
    all.edges %>%
    left_join(all.edges, by="miRNA") %>%  # Gene pairs that are both targeted by the same miRNA
    distinct(Gene.x, Gene.y) %>%  # Only one entry per pair of genes
    filter(Gene.x != Gene.y) %>%  # Remove pairs that are duplicates
    apply(1, sort) %>%  # Sort each row, so that edges are not directed
    t %>% as.data.frame(stringsAsFactors=FALSE) %>% distinct  # Only unique edges
colnames(gene.edges) <- c("from", "to")
    
edgelist.braf <- read.table("BRAF_E.txt", header=TRUE, stringsAsFactors=FALSE)

make_adj_from_df <- function(df, genes) {
    g <- length(genes)
    A <- matrix(0, g, g, dimnames=list(genes, genes))
    for (i in 1:nrow(df)) {
        A[df$to[i],df$from[i]] <- 1
    }
    A
}

A.braf <- make_adj_from_df(edgelist.braf, genes.braf)
A.mirna <- make_adj_from_df(gene.edges, genes.braf)

graph.braf <- graph_from_adjacency_matrix(t(A.braf))
graph.mirna <- graph_from_adjacency_matrix(t(A.mirna), mode="undirected")

## png("gene-networks.png", width = 1600, height = 800)
setEPS()
postscript("gene-networks.eps", width = 16, height = 10,
           fonts="serif")
par(mfrow=c(1, 2))
set.seed(20190116)
lo <- layout_(graph.braf, nicely())
plot(graph.braf, layout=lo,
     vertex.size=20,
     vertex.color="white",
     vertex.weight=3,
     vertex.label.color="black",
     vertex.arrow.width=0,
     edge.width=3,  edge.color="grey20")
plot(graph.mirna, layout=lo,
     vertex.color="white",
     vertex.label.color="black",
     edge.width=3, edge.color="grey20",
     edge.arrow.mode="-")
     ## main="BRAF pathway genes that are co-targets of at least one miRNA")
dev.off()

A.V[rownames(A.braf),colnames(A.braf)] <- A.braf
A <- A.V

save(A, file="A-braf-expr-mirna.Rdata")
