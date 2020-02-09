
## Henry Linder 2020 mhlinder@gmail.com

library(magrittr)
library(dplyr)

innames <- load("bootstrapped-p.Rdata")

p <- list()
for (cancer in names(all.p)) {
    p[[cancer]] <- rbind("Observed p-value"=all.p[[cancer]],
                         "95%"=all.bs.p[[cancer]])
}

lapply(p, function(pp) pp < .05)

## $BLCA
##                  miEMC  miE   EMC     E
## Observed p-value  TRUE TRUE  TRUE FALSE
## 5%                TRUE TRUE  TRUE  TRUE
## 95%               TRUE TRUE FALSE FALSE

## $BRCA
##                  miEMC  miE   EMC    E
## Observed p-value  TRUE TRUE FALSE TRUE
## 5%                TRUE TRUE  TRUE TRUE
## 95%               TRUE TRUE FALSE TRUE

## $HNSC
##                  miEMC  miE   EMC     E
## Observed p-value  TRUE TRUE FALSE FALSE
## 5%                TRUE TRUE  TRUE FALSE
## 95%               TRUE TRUE FALSE FALSE

## $KIRC
##                  miEMC  miE  EMC    E
## Observed p-value  TRUE TRUE TRUE TRUE
## 5%                TRUE TRUE TRUE TRUE
## 95%               TRUE TRUE TRUE TRUE

## $KIRP
##                  miEMC  miE  EMC    E
## Observed p-value  TRUE TRUE TRUE TRUE
## 5%                TRUE TRUE TRUE TRUE
## 95%               TRUE TRUE TRUE TRUE

## $LIHC
##                  miEMC  miE   EMC    E
## Observed p-value  TRUE TRUE  TRUE TRUE
## 5%                TRUE TRUE  TRUE TRUE
## 95%               TRUE TRUE FALSE TRUE

## $THCA
##                  miEMC  miE   EMC     E
## Observed p-value  TRUE TRUE  TRUE FALSE
## 5%                TRUE TRUE  TRUE  TRUE
## 95%               TRUE TRUE FALSE FALSE

do.call(rbind, lapply(p, function(pp) {
    apply(pp < .05, 2, function(x) {
        if (x["Observed p-value"]) {
            length(unique(x[c("Observed p-value", "95%")])) == 1
        } else {
            length(unique(x[c("Observed p-value", "5%")])) == 1
        }
    })
}))

##      miEMC  miE   EMC     E
## BLCA  TRUE TRUE FALSE FALSE
## BRCA  TRUE TRUE FALSE  TRUE
## HNSC  TRUE TRUE FALSE  TRUE
## KIRC  TRUE TRUE  TRUE  TRUE
## KIRP  TRUE TRUE  TRUE  TRUE
## LIHC  TRUE TRUE FALSE  TRUE
## THCA  TRUE TRUE FALSE FALSE
