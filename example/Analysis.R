# Created on 16/03/20 by jdutheil

dat <- read.table("nd5-rates.tsv", header = TRUE)
dat$P.value.adjust <- p.adjust(dat$P.value, method = "fdr")
dat <- dat[order(dat$P.value, decreasing = FALSE),]
