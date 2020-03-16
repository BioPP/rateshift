# Created on 16/03/20 by jdutheil

# Gather resulds for all genes, adjust p-values for multiple testing:

require(data.table)
dat.lst <- list()
for (gene in c("atp6", "atp8", "cox1", "cox2", "cox3", "cytb", "nd1", "nd2", "nd3", "nd4", "nd4L", "nd5", "nd6")) {
  dat.gene <- read.table(paste(gene, "-rates.tsv", sep = ""), header = TRUE)
  dat.gene$Gene <- gene
  dat.gene$P.value.adjust <- p.adjust(dat.gene$P.value, method = "fdr")
  dat.lst[[gene]] <- dat.gene
}
dat <- rbindlist(dat.lst)

# Output top results:
dat.signif <- subset(dat, P.value < 0.05)
