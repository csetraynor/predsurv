library(ggplot2)
library(plyr)
library(grid)
library(extrafont)
library(cowplot)
library(ggpubr)

DCF_Pre <- read.table("//mokey.ads.warwick.ac.uk/User41/u/u1795546/Desktop/DCF_Gem_Pre.txt",head=TRUE)
#Data: ID,  DOSE, TIME, DV



DCF_Pre$TIME <- as.factor(DCF_Pre$TIME)
DCF_Pre$DOSE <- as.factor(DCF_Pre$DOSE)


# unequal variation = default, unpaired = default
compare_means(data=DCF_Pre,DV~DOSE, method="t.test",ref.group = "0") 
#compare_means(data=DCF_Pre,DV~TIME, method="t.test",ref.group = "0")

theme_set(theme_cowplot(font_size=8,font_family = "Times New Roman"))


fig1 <- ggbarplot(data=DCF_Pre,x="DOSE",
                  y="DV",
                  add="mean_se",
                  position=position_dodge(0.8),
                  fill = "lightgrey",
                  width=0.8)+
  stat_compare_means(method="t.test",
                     label="p.signif",
                     label.y=90,
                     hide.ns="TRUE",
                     ref="0")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8)
  )+
  xlab("Gemfibrozil incubation concentration (nmole/ml)")+
  ylab("% of DCF control")


fig2 <- ggdraw(fig1)+
  theme(rect=element_rect(fill="white"))
# fig3 <- add_sub(fig2,"Gemfibrozil incubation concentration (nmoles/ml)",hjust=0.5,vjust=0,size=8) 
# fig4 <- add_sub(fig3,"% of DCF control",hjust=0,vjust=0.5,size=8,angle=90)
ggsave("DCF_Pre.jpg",fig2,width=6,height=3,units="in",dpi=300)        



# file://mokey.ads.warwick.ac.uk/User41/u/u1795546/Desktop/DCF_Gem_Pre.txt
###############################################################################################

res <- ddply(DCF_Pre,.(TIME,DOSE),
             summarise,
             meani = mean(DV),
             se = sd(DV)/2)

res$ymin = res$meani-res$se
res$ymax = res$meani + res$se

compare_means(data=DCF_Pre,DV~DOSE, method="t.test",ref.group = "0")

fig1 <- ggplot(data=res,aes(y=meani,x=DOSE,ymin=ymin,ymax=ymax,fill=TIME))+
        geom_bar(stat="identity",position="dodge",width=0.8)+
        geom_errorbar(width=0.3,position=position_dodge(0.8))+
        scale_color_brewer(palette="Dark2")+
        scale_fill_brewer(palette="Dark2")

fig1
