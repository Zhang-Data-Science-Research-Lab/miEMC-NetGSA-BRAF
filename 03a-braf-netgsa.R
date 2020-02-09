
## Henry Linder 2020 mhlinder@gmail.com

data.path <- "~/local/data/TCGA/formatted"
out.path <- "gsaout"

library(magrittr)
library(dplyr)
library(readr)

library(RColorBrewer)

library(ppcor)
library(netgsa)

source("color.R")

innames <- load("A-braf-expr-mirna.Rdata")

## Names of omics features
all.names <- dimnames(A)[[1]]
## Names of miRNA
mirna.A <- all.names[grepl("^hsa-", all.names)]

## Identify miRNA which have both 3p and 5p, keep only one, and remove
## its identifier as being "-3p" or "-5p": TGCA does not include this
## information
rm.col <- gsub("-[35]p", "", mirna.A) %>% duplicated %>% which
rm.col <- which(colnames(A) %in% mirna.A[rm.col])
A <- A[-rm.col,-rm.col]
rownames(A) <- gsub("-[35]p", "", rownames(A))
colnames(A) <- gsub("-[35]p", "", colnames(A))

## Gene and miRNA symbols
all.names <- dimnames(A)[[1]]
mirna.A <- all.names[grepl("^hsa-", all.names)]
genes.A <- all.names[!all.names %in% mirna.A]

## All miRNAs in TCGA data
mirna.tcga <- read_csv("all-mirna.csv") %>% use_series(miRNA)

## Convert mirDIP identifiers to their corresponding TCGA format
mirdip.to.tcga <- function(x) {
    x %>%
        gsub("^hsa-", "", .) %>%
        gsub("^mir-", "MIR", .) %>%
        gsub("^let-", "LET", .) %>%
        toupper
}
## Get the TCGA labels for the mirDIP miRNA tags
mirna.tcga.as.des <- mirdip.to.tcga(mirna.A)

if (!all(mirna.A %in% mirna.tcga)) {
    stop("Something got mismatched!")
}

## All feature names
feature.names <- c(genes.A, mirna.tcga.as.des)

infiles <- list.files(data.path, full.names=TRUE)
nfiles <- length(infiles)

source("functions.R")

for (ix.cancer in 1:length(infiles)) {
    infile <- infiles[ix.cancer]
    cancer <- infile %>% basename %>% strsplit("-") %>% sapply(head, 1)
    print(sprintf("[%d / %d] %s", ix.cancer, nfiles, cancer))
    innames <- load(infile)

    ## miRNA tags
    mirna.des <- des %>% filter(Platform == "miRNAExp") %>% use_series(GeneSymbol)

    if (!all(mirna.tcga.as.des %in% mirna.des)) print(sprintf("Missing features in %s", cancer))

    ## Data is a miRNA or a gene we are interested in
    des <- des %>% filter(GeneSymbol %in% feature.names)
    all.features <- list()
    ## Check which data platforms are available for each feature
    for (fn in feature.names) {
        d <- des %>% filter(GeneSymbol == fn)
        dd <- data.frame(GeneSymbol=fn,
                         Expression="geneExp" %in% d$Platform,
                         Copy.Number="copyNumber" %in% d$Platform,
                         Methylation="methylation" %in% d$Platform,
                         mi="miRNAExp" %in% d$Platform,
                         stringsAsFactors=FALSE)
        all.features[[length(all.features)+1]] <- dd
    }
    all.features <- bind_rows(all.features)

    ## Get the data, name the rows
    dat <- dat[des$ix,]
    dat <- dat[,!apply(dat, 2, . %>% is.na %>% any)]
    dat <- dat[,grepl("-[01]1", colnames(dat))]

    ## Remove missing values
    rm.rows <- dat %>% apply(1, . %>% is.na %>% any) %>% which
    if (length(rm.rows) > 0) {
        rm.counts <- dat %>% apply(1, . %>% is.na %>% sum)
        rm.counts <- rm.counts[rm.rows]
        dat <- dat[-rm.rows,]
        ##
        rm.des <- des %>% filter(tag %in% names(rm.rows))
        des <- des %>% filter(!tag %in% rm.des$tag)
    }

    ## Remove those features with too few unique values in each population
    nunique <- apply(dat[,grepl("-01$", colnames(dat)), drop=FALSE], 1, . %>% unique %>% length)
    rm.rows1 <- which(nunique < 5)

    nunique <- apply(dat[,grepl("-11$", colnames(dat)), drop=FALSE], 1, . %>% unique %>% length)
    rm.rows2 <- which(nunique < 5)

    rm.rows <- c(rm.rows1, rm.rows2) %>% unique

    if (length(rm.rows) > 0) {
        dat <- dat[-rm.rows,]
        keep.des <- des %>% filter(tag %in% rownames(dat))
        des <- des %>% filter(tag %in% keep.des$tag)
    }

    ## Rename the dimensions of A to include information on genomic
    ## platform
    mirna.tcga.labels <- mirdip.to.tcga(mirna.A)
    feature.labels <- c(paste0(genes.A, "-GE"),
                        paste0(mirna.tcga.labels, "-mi"))

    A.cancer <- A
    dimnames(A.cancer) <- list(feature.labels, feature.labels)

    ## number of miRNA per gene
    ## A.cancer[grepl("-GE$", rownames(A.cancer)),grepl("-mi$", colnames(A.cancer))] %>% rowSums
    ##   AKT1-GE   BRAF-GE MAP2K1-GE MAP2K2-GE  MAPK1-GE   MTOR-GE   NRAS-GE PIK3CA-GE 
    ##        14         0        40         0        88        15        61        26 
    ##   PTEN-GE   RAF1-GE 
    ##       130        16 
    ## 
    ## number of genes targeted by each miRNA
    ## A.cancer[grepl("-GE$", rownames(A.cancer)),grepl("-mi$", colnames(A.cancer))] %>% colSums %>% table
    ##    1   2   3   4   5 
    ##  108  56  31  18   1 
    ## g=214 total miRNAs

    ## Expand A to include miRNA, EMC features
    n.genes <- length(genes.A)
    n.mirna <- length(mirna.A)
    emc <- rbind(diag(n.genes),
                 matrix(0, n.mirna, n.genes))
    emc.cn <- emc
    colnames(emc.cn) <- paste0(genes.A, "-CN")
    emc.me <- emc
    colnames(emc.me) <- paste0(genes.A, "-ME")

    A1 <- cbind(A.cancer, emc.cn, emc.me)
    A2 <- matrix(0,
                 nrow=2*n.genes,
                 ncol=ncol(A) + 2*n.genes,
                 dimnames=list(c(colnames(emc.cn), colnames(emc.me)),
                               colnames(A1)))

    A.cancer <- rbind(A1, A2)

    ix.mir <- grepl("-mi$", colnames(A.cancer)) %>% which
    ix.cn <- grepl("-CN$", colnames(A.cancer)) %>% which
    ix.me <- grepl("-ME$", colnames(A.cancer)) %>% which
    ix.ge <- grepl("-GE$", colnames(A.cancer)) %>% which

    ix <- c(ix.mir, ix.cn, ix.me, ix.ge)

    A.cancer <- A.cancer[ix, ix]

    ## Features we removed because of missing values or insufficiently
    ## many unique values
    feat.A <- rownames(A.cancer)
    feat.dat <- rownames(dat)
    rm.feat <- which(!feat.A %in% feat.dat)
    A.cancer <- A.cancer[-rm.feat,-rm.feat]
    if (!all(sprintf("%s-GE", genes.A) %in% rownames(A.cancer))) {
        print("Too many missing values")
        next
    }

    x <- dat[rownames(A.cancer),]
    if (ncol(x) == 0) {
        print("no subjects")
        next
    }

    keep.subj <- grepl("-[01]1$", colnames(x)) %>% which
    x <- x[,keep.subj]

    y <- colnames(x) %>%
        strsplit("-") %>%
        sapply(tail, 1) %>%
        factor(levels=c("11", "01")) %>%
        as.numeric

    y.tab <- y %>% table %>% c
    if ((length(y.tab) == 1) |
        (!all(y.tab>10))) {
        print("not enough subjects")
        next
    }

    estimate_adj <- function(A, x, y) {
        n_tags <- nrow(A)
        node_tags <- rownames(A)
        A1 <- A2 <- A

        for (j in 1:n_tags) {
            tag <- node_tags[j]
            if (sum(A[j,]) > 0) {
                parents_tags <- node_tags[as.logical(A[j,])]
                correlation.data1 <- x[c(tag, parents_tags), y == 1, drop = FALSE] %>%
                    t %>% unname %>% data.frame
                rm_ix1 <-
                    correlation.data1 %>%
                    apply(1, . %>% is.na %>% any) %>%
                    which
                if (length(rm_ix1) > 0)
                    correlation.data1 <- correlation.data1[-rm_ix1,]

                correlation.data2 <- x[c(tag, parents_tags), y == 2, drop = FALSE] %>%
                    t %>% unname %>% data.frame
                rm_ix2 <-
                    correlation.data2 %>%
                    apply(1, . %>% is.na %>% any) %>%
                    which
                if (length(rm_ix2) > 0)
                    correlation.data2 <- correlation.data2[-rm_ix2,]

                cor.out1 <- pcor(correlation.data1)
                cor.out2 <- pcor(correlation.data2)

                partial.corr1 <- cor.out1[["estimate"]]
                partial.corr2 <- cor.out2[["estimate"]]
                rownames(partial.corr1) <- colnames(partial.corr1) <- c(tag, parents_tags)
                rownames(partial.corr2) <- colnames(partial.corr2) <- c(tag, parents_tags)

                A1[j, node_tags %in% parents_tags] <- partial.corr1[tag, parents_tags]
                A2[j, node_tags %in% parents_tags] <- partial.corr2[tag, parents_tags]
            }
        }
        ## A1[is.na(A1)] <- 0
        ## A2[is.na(A2)] <- 0

        return(list(A1 = A1, A2 = A2))
    }

    AA.miEMC <- estimate_adj(A.cancer, x, y)
    subset.miEMC <- list(AA=AA.miEMC,
                         x=x)

    make.subset <- function(regex, negate=TRUE, x) {
        ## Removes rows that match the regex
        ix <- grepl(regex, rownames(x))
        ix.row <- grepl(regex, rownames(A.cancer))
        ix.col <- grepl(regex, colnames(A.cancer))
        if (negate) {
            ix <- not(ix)
            ix.col <- not(ix.col)
            ix.row <- not(ix.row)
        }
        xx <- x[ix,]
        list(AA=estimate_adj(A.cancer[ix.row, ix.col],
                             xx, y),
             x=xx)
    }

    subset.EMC <- make.subset("(^MIR|^LET|-mi$)", x=x)
    subset.E <- make.subset("(^MIR|^LET|-mi$|-CN$|-ME$)", x=x)
    subset.miE <- make.subset("(^MIR|^LET|-GE$)", negate=FALSE, x=x)

    subsets <- list(miEMC=subset.miEMC,
                    miE=subset.miE,
                    EMC=subset.EMC,
                    E=subset.E)
    all.gsa <- list()
    for (i in 1:length(subsets)) {
        print(sprintf("%s", names(subsets)[i]))
        AA <- subsets[[i]]$AA
        xx <- subsets[[i]]$x
        ## A regular expression for any of the gene names, followed by "-GE"---eg,
        ## BRAF-GE
        re.ge <- genes.A %>% paste(collapse="|") %>% sprintf("(%s)-GE", .)
        ## Any gene 
        re.genes <- genes.A %>% paste(collapse="|") %>% sprintf("(%s)-", .)
        B <- matrix(c(1*grepl(re.ge, rownames(AA[[1]])),
                      1*grepl(re.genes, rownames(AA[[1]]))),
                    byrow=TRUE, nrow=2)
        rownames(B) <- c("BRAF-GE", "BRAF-EMC")
        gsa <- NetGSA(AA, xx, y, B=B, directed=TRUE, lklMethod="REML")
        rownames(gsa$teststat) <- rownames(gsa$df) <- rownames(gsa$p.value) <- rownames(B)
        all.gsa[[i]] <- gsa
    }
    names(all.gsa) <- names(subsets)

    p <- do.call(cbind, lapply(all.gsa, function(g) g$p.value))
    colnames(p) <- names(all.gsa)

    T <- do.call(cbind, lapply(all.gsa, function(g) g$teststat))
    colnames(T) <- names(all.gsa)

    s2.e <- sapply(all.gsa, function(g) g$s2.epsilon)
    s2.g <- sapply(all.gsa, function(g) g$s2.gamma)

    beta1 <- list()
    beta2 <- list()
    for (gn in names(all.gsa)) {
        gsa <- all.gsa[[gn]]
        beta1[[gn]] <- gsa$beta[[1]]
        beta2[[gn]] <- gsa$beta[[2]]
    }
    all.est <- list(beta1 = beta1,
                    beta2 = beta2,
                    s2.e = s2.e,
                    s2.g = s2.g,
                    p = p,
                    T = T)

    gsa <- all.gsa
    est <- all.est
    dat <- subsets$miEMC$x

    ## Save output from netgsa
    save(gsa, est, dat,
         file=sprintf("out/%s/out-%s.Rdata", out.path, cancer))

    rm(all.gsa, all.est, subsets,
       gsa, est, dat)
}
