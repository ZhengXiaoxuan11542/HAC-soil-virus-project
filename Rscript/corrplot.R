library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

setwd("E:/HACscript")
#read data
aa=read.table("Data/VLP_data.txt",header=TRUE,row.names = 1)


### 出图2 
p1<-ggplot(aa,aes(x=lytic_per,y=VLP))+
geom_point(size=4,aes(color=treat))+
scale_color_manual(values=c("#c5ccae","#ffc97b", "#167eb2"))+
geom_smooth(method=lm,level=0.95,color="gray4")+
stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
theme_classic()+
labs(y="VLP",x="Lytic percentage(%)")+
theme(axis.title = element_text(color='black',size=9),
axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
axis.line = element_line(colour = "black"),
axis.title.x=element_text(colour='black', size=11),
axis.title.y=element_text(colour='black', size=11),
axis.text=element_text(colour='black',size=9),
legend.title=element_blank(),
legend.text=element_text(size=9),
legend.key=element_blank(),
legend.background = element_rect(colour = "White"))
#p1
# the pearson result
cor.test(aa$lytic_per,aa$VLP,method = "pearson")


p2<- ggplot(aa, aes(x = lytic_per, y = X16s_copies))+
geom_point(size=4,aes(color=aa$treat))+
scale_color_manual(values=c("#c5ccae","#ffc97b", "#167eb2"))+
geom_smooth(method=lm,level=0.95,color="gray4")+
stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
theme_classic()+
labs(y="16s_copies",x="Lytic percentage(%)")+
theme(axis.title = element_text(color='black',size=9),
axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
axis.line = element_line(colour = "black"),
axis.title.x=element_text(colour='black', size=11),
axis.title.y=element_text(colour='black', size=11),
axis.text=element_text(colour='black',size=9),
legend.title=element_blank(),
legend.text=element_text(size=9),
legend.key=element_blank(),
legend.background = element_rect(colour = "White"))
# the pearson result
cor.test(aa$lytic_per,aa$X16s_copies,method = "pearson")


p3<- ggplot(aa, aes(x = lytic_per, y = VMR))+
geom_point(size=4,aes(color=aa$treat))+
scale_color_manual(values=c("#c5ccae","#ffc97b", "#167eb2"))+
geom_smooth(method=lm,level=0.95,color="gray4")+
stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
theme_classic()+
labs(y="VMR",x="Lytic percentage(%)")+
theme(axis.title = element_text(color='black',size=9),
axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
axis.line = element_line(colour = "black"),
axis.title.x=element_text(colour='black', size=11),
axis.title.y=element_text(colour='black', size=11),
axis.text=element_text(colour='black',size=9),
legend.title=element_blank(),
legend.text=element_text(size=9),
legend.key=element_blank(),
legend.background = element_rect(colour = "White"))
# the pearson result
cor.test(aa$lytic_per,aa$VMR,method = "pearson")

combined_plots <-p1+p2+p3+
  patchwork::plot_layout(design = "
                         A
                         B
                         C
                         ")
combined_plots
ggsave(plot = combined_plots, filename = "Figure/corr_with_lyticphage.pdf",  width = 140, height = 300, units = "mm", dpi = 300)




#read data
bb=read.table("Data/MAG-vOTU_cor_data.txt",header=TRUE,row.names = 1)
### 出图2 
p4<-ggplot(bb,aes(x=vOTUs.trimmean.,y=MAGs.trimmean.))+
  geom_point(size=4,aes(color=treat))+
  scale_color_manual(values=c("#c5ccae","#ffc97b", "#167eb2"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
  theme_classic()+
  labs(y="MAGs.trimmean.",x="vOTUs.trimmean.")+
  theme(axis.title = element_text(color='black',size=9),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=11),
        axis.title.y=element_text(colour='black', size=11),
        axis.text=element_text(colour='black',size=9),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        legend.key=element_blank(),
        legend.background = element_rect(colour = "White"))
p4
ggsave(plot = p4, filename = "Figure/vOTU_MAG_corr.pdf",  width = 140, height = 100, units = "mm", dpi = 300)



