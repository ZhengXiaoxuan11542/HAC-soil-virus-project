####23redo#####
#read data
library(dplyr)
library(picante)
library(ggplot2)
library(ggsignif)
library(patchwork)

setwd("E:/HACscript")
d=read.table("Data/votu_table.txt",header=TRUE,row.names = 1)


d1<-d[,(1:1267)]
group<-d[,c(1,1268)]
env<-d[,c(1,1269)]

# delete singleton
d1 <- d1[,colSums(d1>0)>1]
d1 <- d1[rowSums(d1)>0,]


#rel_abun and translate values into integer by x 10000
rel_ab <- function(otu) {
  otu/rowSums(otu) * 100
}

d_norm<-d1 %>% rel_ab
rowSums(d_norm)

alpha <- function(x, base = exp(1)) {
  x1 =x
  x1[x1 > 0] <- 1
  Shannon <- diversity(x, index = 'shannon', base = base)
  result <- data.frame(Shannon)
  result
}

alpha_all <- alpha(d_norm)
#write.csv(alpha_all,file = "v_alpha_div.csv",quote = F)
data_ggplot=data.frame(alpha_all, env["env"], group["group"])

p1<-ggplot(data_ggplot,aes(x=group,y=Shannon,fill=group))+geom_boxplot()+ 
  labs(title="Shannon", x="Group", y="Shannon")+ 
  scale_fill_manual(values = c("#c5ccae","#ffc97b", "#167eb2"))+ 
  theme(plot.title=element_text(hjust=0.5), legend.title=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme_classic()

combined_plots <- p1 + plot_annotation(title = "Virus alpha diversity")
combined_plots
ggsave(plot = combined_plots, filename = "Figure/v_alpha.pdf",  width = 120, height = 150, units = "mm", dpi = 300)

#significance test
wilcox.test(Shannon ~ env, data = data_ggplot)
wilcox.test(Shannon ~ group,data = subset(data_ggplot, group %in% c("Con", "HAC")))

########pcoa###########
library(patchwork)
library(plyr)
library(tibble)
library(dplyr)
library(dbplyr)
library(vegan)
library(ggplot2)
library(ade4)


find_hull_pcoa <- function(df) df[chull(df$PCoA1, df$PCoA2), ]

pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

by.dist<-vegdist(d_norm,method = 'bray')
by.dist = as.data.frame(as.matrix(by.dist)) # output the bray matrix
#write.csv(by.dist,file = "v_pcoa_bc_matrix.csv",quote = F)

pcoa_by<-dudi.coa(by.dist,scan = FALSE, nf = 10)
pcoa_by_eig<-(pcoa_by$eig)[1:10] / sum(pcoa_by$eig)
pcoa_by_data<-data.frame(pcoa_by$li)
names(pcoa_by_data)[1:2]<-c('PCoA1','PCoA2')
pcoa_by_data$group = d $group
view(pcoa_by_data) #to check
adonis2(formula = d_norm ~ group, data = d, permutations = 999,method = "bray")
pairwise.adonis(d_norm,  group$group, sim.method="bray", p.adjust.m= "fdr")
pcoa<-ggplot(pcoa_by_data,aes(PCoA1, PCoA2,color = group))+stat_ellipse(level = 0.8, show.legend = F) +
  theme_bw() + labs(title = "Viral beta diversity", subtitle = "R2 = 0.398, P=0.001**")+
  geom_point(aes(color = group, shape = group), size =4, alpha = 0.8) + #可在这里修改点的透明度、大小
  scale_shape_manual(values = c(17,17,17)) +
  scale_color_manual(values = c("#c5ccae","#ffc97b", "#167eb2")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(x = paste('PCoA1:', round(100*pcoa_by_eig[1],2),'%'),y=paste('PCoA2:',round(100*pcoa_by_eig[2],2),'%'))
pcoa

ggsave(plot = pcoa, filename = "Figure/v_pcoa.pdf",  width = 106, height = 94, units = "mm", dpi = 300)
#write.csv(pcoa_by_data,file = "v_pcoa_bc.csv",quote = F)

#rhizo
rhizo = d_norm[c(7:18),]
group.rhizo = group[c(7:18),]
by.dist.rhizo<-vegdist(rhizo,method = 'bray')
by.dist.rhizo = as.data.frame(as.matrix(by.dist.rhizo)) # output the bray matrix
#write.csv(by.dist.rhizo,file = "v_rhizo_pcoa_bc_matrix.csv",quote = F)
pcoa_by.rhizo<-dudi.coa(by.dist.rhizo,scan = FALSE, nf = 10)
pcoa_by_eig.rhizo<-(pcoa_by.rhizo$eig)[1:10] / sum(pcoa_by.rhizo$eig)
pcoa_by_data.rhizo<-data.frame(pcoa_by.rhizo$li)
names(pcoa_by_data.rhizo)[1:2]<-c('PCoA1','PCoA2')
pcoa_by_data.rhizo$group = group.rhizo $group
#view(pcoa_by_data.rhizo) #to check
adonis2(formula = rhizo ~ group, data = group.rhizo, permutations = 999,method = "bray")
pcoa_rhizo<-ggplot(pcoa_by_data.rhizo,aes(PCoA1, PCoA2,color = group))+stat_ellipse(level = 0.8, show.legend = F) +
  theme_bw() + labs(title = "Viral rhizosphere beta diversity", subtitle = "R2 = 0.114, P=0.197")+
  geom_point(aes(color = group, shape = group), size =4, alpha = 0.8) + #可在这里修改点的透明度、大小
  scale_shape_manual(values = c(17,17)) +
  scale_color_manual(values = c("#ffc97b", "#167eb2")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(x = paste('PCoA1:', round(100*pcoa_by_eig.rhizo[1],2),'%'),y=paste('PCoA2:',round(100*pcoa_by_eig.rhizo[2],2),'%'))
pcoa_rhizo
ggsave(plot = pcoa_rhizo, filename = "Figure/v_pcoa_rhizo.pdf",  width = 106, height = 94, units = "mm", dpi = 300)
#write.csv(pcoa_by_data.rhizo,file = "v_rhizo_pcoa_bc.csv",quote = F)


