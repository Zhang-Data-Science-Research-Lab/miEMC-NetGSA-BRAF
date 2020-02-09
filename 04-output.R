
## Henry Linder 2020 mhlinder@gmail.com

library(magrittr)
library(dplyr)
library(readr)
library(xtable)

source("color.R")

infiles <- list.files("out/gsaout", full.names=TRUE)

cancers <- rep(NA_character_, length(infiles))
all.est <- all.gsa <- all.dat <- list()
for (i in 1:length(infiles)) {
    infile <- infiles[i]
    cancer <- infile %>% basename %>% gsub("out-|\\.Rdata", "", .)
    load(infile)
    all.est[[i]] <- est
    all.dat[[i]] <- dat
    all.gsa[[i]] <- gsa
    cancers[i] <- cancer
}

ss <-
    lapply(all.dat,
           function(d)
               substr(colnames(d), 14, 16) %>%
               table %>%
               c %>%
               t %>% 
               as.data.frame(stringsAsFactors=FALSE)) %>%
    bind_rows %>%
    rename(Tumor=`01`,
           Normal=`11`)
ss$Cancer <- cancers
    
cancernames <- read_tsv("cancers.txt")
ss <- ss %>% left_join(cancernames, by=c("Cancer"="Code"))

ss %>%
    select(Cancer, Name, Tumor, Normal) %>%
    rename(Code=Cancer,
           Cancer=Name) %>%
    mutate(Cancer=tolower(Cancer),
           Cancer=paste0(toupper(substr(Cancer, 1, 1)),
                         substr(Cancer, 2, nchar(Cancer)))) %>%
    select(Code, Cancer, Tumor, Normal) %>% 
    xtable %>%
    print(booktabs=TRUE,
          include.rownames=FALSE,
          table.placement=NULL)

round.out <- function (x, nearest = 5) {
    down <- nearest * floor(x[1]/nearest)
    up <- nearest * ceiling(x[2]/nearest)
    c(down, up)
}
cols <- c(cols_all[c("GreenB1", "Yellow1", "Blue1")], "white")

p <- lapply(1:length(all.est), function(i) {
    est <- all.est[[i]]
    pp <- est$p["BRAF-GE",]
    o <- data.frame(p=pp,
                    stringsAsFactors=FALSE)
    colnames(o) <- cancers[i]
    o
}) %>% bind_cols
rownames(p) <- colnames(all.est[[1]]$p)
p <- as.matrix(p)

p <- apply(p, 2, p.adjust, method="BH")
## p<.05

p <- -log10(p)


setEPS()
postscript("Tbar-pbar.eps", width = 12, height = 10)

## png("pbar.png", width = 800, height = 250)
## setEPS()
## postscript("pbar.eps", width = 12, height = 5)
par(lwd=1.5, mar=c(6.1, 5.1, 2.1, 2.1),
    cex.axis=1.5, cex.lab=2,
    cex.main=1.5,
    mfrow=c(2, 1))
yr <- round.out(range(p), 5)
bp <- barplot(p, beside=TRUE, 
        col=cols, cex.names=1.5,
        ylim=yr, ylab=expression(paste("-log"[10], "(p)")),
        main=expression(paste("-log"[10], "(p)-values")))
abline(h=0)
abline(h=-log10(.05), col="black", lty=2)
splits <- sapply(2:ncol(bp), function(i) mean(c(tail(bp[,i-1], 1),
                                                head(bp[,i], 1))))
diffs <- sapply(2:ncol(bp), function(i) diff(c(tail(bp[,i-1], 1),
                                               head(bp[,i], 1)))) %>%
    unique
splits <- c(bp[1]-diffs/2, splits, bp[length(bp)]+diffs/2)
abline(v=splits, col="grey60", lty=2)
box()
legend("topright", fill=c(cols, NA), bg="white",
       lty=c(rep(NA, 4), 2),
       lwd=c(rep(NA, 4), 2),
       col=c(rep("black", 4), "black"),
       border=c(rep("black", 4), NA),
       c(rownames(p), expression(paste(alpha, "=0.05"))),
       cex=1.5, inset=c(.01, .02))

T <- lapply(1:length(all.est), function(i) {
    est <- all.est[[i]]
    TT <- est$T["BRAF-GE",]
    o <- data.frame(T=TT,
                    stringsAsFactors=FALSE)
    colnames(o) <- cancers[i]
    o
}) %>% bind_cols
rownames(T) <- colnames(all.est[[1]]$T)
T <- as.matrix(T)

yr <- round.out(range(T), 5)

## png("Tbar.png", width = 800, height = 250)
## setEPS()
## postscript("Tbar.eps", width = 12, height = 5)
par(mar=c(3.1, 5.1, 4.1, 2.1))
##     cex.axis=1.5, cex.lab=2,
##     mfrow=c(2, 1)))

bp <- barplot(T, beside=TRUE, 
        col=cols, cex.names=1.5,
        ylim=yr, ylab="T",
        main=expression("Test statistics"))
abline(h=0)
splits <- sapply(2:ncol(bp), function(i) mean(c(tail(bp[,i-1], 1),
                                                head(bp[,i], 1))))
diffs <- sapply(2:ncol(bp), function(i) diff(c(tail(bp[,i-1], 1),
                                               head(bp[,i], 1)))) %>%
    unique
splits <- c(bp[1]-diffs/2, splits, bp[length(bp)]+diffs/2)
abline(v=splits, col="grey60", lty=2)
box()
legend("bottomright", fill=cols, bg="white",
       rownames(T), cex=1.5, inset=c(.01, .02))
dev.off()


