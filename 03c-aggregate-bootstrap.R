
## Henry Linder 2020 mhlinder@gmail.com

library(magrittr)
library(dplyr)

raw.out <- list.files("~/local/data/TCGA/formatted", full.names=TRUE)
infiles <- list.files("out/bootstrap", full.names=TRUE)

all.p <- list()
all.bs.p <- list()
reps <- list()
raw.gsa <- list()
for (i in 1:length(infiles)) {
    infile <- infiles[i]
    cancer <- basename(infile) %>% gsub("out-|\\.Rdata", "", .)
    print(sprintf("%d / %d - %s", i, length(infiles), cancer))
    innames <- load(infile)
    reps[[cancer]] <- bootstrap.out
    ##
    rawfile <- raw.out[grepl(cancer, raw.out)]
    innames <- load(rawfile)
    ##
    raw.gsa[[cancer]] <- list(gsa=gsa,
                              est=est)
    ##
    ## B <- length(bootstrap.out)
    ## bs.T.ge <- sapply(bootstrap.out, function(x) x$est$T["BRAF-GE",])
    ## T.ge <- est$T["BRAF-GE",]
    ## ##
    ## all.bs.p[[cancer]] <- 2*sapply(1:4, function(i) min(sum(bs.T.ge[i,] < T.ge[i]), mean(bs.T.ge[i,] > T.ge[i])))
    bs.p <- sapply(bootstrap.out, function(x) x$est$p["BRAF-GE",])
    bs.p <- apply(bs.p, 1, quantile, probs=c(.05, .95))
    all.bs.p[[cancer]] <- bs.p
    all.p[[cancer]] <- est$p["BRAF-GE",]
}

save(all.bs.p, all.p, file="bootstrapped-p.Rdata")
