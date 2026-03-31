
setwd("E:/HACscript")
library(corrplot)
library(psych)
library(ggplot2)


# 读取数据
state.clean  <- read.table("Data/vb_gene_corr_data.txt", header = TRUE, row.names = 1)
dput(state.clean)

rhizo.spearman.clean<-cor(state.clean,use = "complete.obs",method = "spearman")
res5 <- cor.mtest(state.clean, conf.level = .95,  method = "spearman") 
genes_v <- c("UGDH", "galE", "K16150")
genes_b <- c("rifA", "rifB", "rifC_D", "rifL", "rifN")

corr_sub <- rhizo.spearman.clean[genes_v, genes_b]
p_sub    <- res5$p[genes_v, genes_b]


#figure output
pdf("Results/vb_gene_corr_result.pdf", width = 500/80, height = 350/80)

corrplot(
  corr = corr_sub,
  p.mat = p_sub,
  col = colorRampPalette(
    c("#3B5B92", "#F7F7F7", "#B24C4C")
  )(200),
  method = "color",
  insig = "label_sig",
  sig.level = c(.001, .01, .05),
  pch.cex = 3,
  tl.cex = 1.8,
  cl.cex = 1.5,
  title = "Spearman correlation",
  title.cex = 2
)
dev.off()

#data output
corr_df <- expand.grid(genes_v,genes_b, stringsAsFactors = FALSE)

corr_df$spearman.r <- as.vector(corr_sub)
corr_df$spearman.p <- as.vector(p_sub)

# rename
colnames(corr_df) <- c("virus.gene", "bacteria.gene", "spearman.r", "spearman.p")

#write.table( corr_df,file = "Results/rhizo.spearman.r.txt", sep = "\t", row.names = FALSE, quote = FALSE)
