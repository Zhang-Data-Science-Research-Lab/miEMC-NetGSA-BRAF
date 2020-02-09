
## Henry Linder 2020 mhlinder@gmail.com


plot.beta <- function(beta1, beta2, genes.A) {
    all.ge1 <- sapply(beta1, function(b) b[grepl("-GE$", rownames(b))])
    all.ge2 <- sapply(beta2, function(b) b[grepl("-GE$", rownames(b))])
    rownames(all.ge1) <- rownames(all.ge2) <- genes.A

    cols <- cols_all[c("Red2", "Blue2", "GreenB2", "Orange2")]
    yr <- c(min(0, all.ge1, all.ge2), max(all.ge1, all.ge2, 0))

    par(mfrow=c(1, 2))
    barplot(t(all.ge1), beside=TRUE,
            col=cols, ylim=yr,
            main="Healthy tissue",
            las=2)
    abline(h=0)
    barplot(t(all.ge2), beside=TRUE,
            col=cols, ylim=yr,
            main=sprintf("Tumor - %s", cancer),
            las=2)
    abline(h=0)
    legend("topright", bty="n",
           fill=cols, legend=colnames(all.ge1))
}
