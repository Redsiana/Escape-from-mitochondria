########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##    PLOT Fig 3 of two thresholds    ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
# 
# Produces Fig. 3 plots of the 2 threshold rates of fusion (polymorphism, and
# fixation) across costs of fusion, for a given number of mitochondria
# Action required: pick numbers of mitochondria (v.Mitograph) to investigate
# 

library(RColorBrewer)
library(ggplot2)
source("multiplot.R")

############### PICK A NUMBER OF MITOCHONDRIA

v.Mitograph <- c(4, 20, 200)


dataset1 <- read.table("titration_fix_df.txt", h=T)
dataset2 <- read.table("titration_poly_df.txt", h=T)

gg_list <- list()
x_list <- list("", "relative reproductive success of fusers (1-c)", "")
y_list <- list("rate of fusion", "", "")


for(i in 1:(length(v.Mitograph))){
  
  df.fix <- dataset1[dataset1$Mito == v.Mitograph[i],]
  df.poly <- dataset2[dataset2$Mito == v.Mitograph[i],]
  df <- data.frame(rate_sex = df.fix$fert_sex, 
                   Mito = df.fix$Mito,
                   rate_fix = df.fix$result,
                   rate_poly = df.poly$result)
  df[df==-Inf] <- 0
  
  mypalette <- colorRampPalette( brewer.pal(9, 'PuBuGn'))
  
  mypalette <- brewer.pal(9, 'PuBuGn')
  newcol <- colorRampPalette(mypalette)
  mypalette <- newcol(100)
  
  gg_list[[i]] = ggplot(df, aes(x=rate_sex)) +
    geom_line(aes(y=rate_poly)) + geom_area( aes(y=rate_poly), fill = "#E8DFEE") +
    geom_line(aes(y=rate_fix))  + geom_area( aes(y=rate_fix), fill = "#014636") +
    scale_x_continuous(breaks = seq(0, 1, 0.1),
                       x_list[[i]]) +
    scale_y_continuous(breaks = seq(0,1,0.1), 
                       y_list[[i]]) +
    geom_text(data=df.fix, aes(x=0.9, y=0.075, label='FUSION\nFIXES'), 
              color = "salmon", fontface = "bold") +
    geom_text(data=df.fix, aes(x=0.75, y=0.25, label='POLYMORPHIC\nEQUILIBRIUM'), 
              color = "salmon", fontface = "bold") +
    geom_text(data=df.fix, aes(x=0.4, y=0.6, label='FUSION DOES NOT\nEVOLVE'), 
              color = "salmon", fontface = "bold") +
    geom_text(data=df.fix, x=0.1, y=1,
                               label=paste("M =", 
                                           v.Mitograph[i]), color = 'black', 
              fontface = "bold") +
    theme(#panel.background = element_rect(fill = 'darkgoldenrod1'), 
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      aspect.ratio=4/3.5)
  
}


dev.new(width=300, height=10)
pdf('Fig3.pdf'  )
multiplot(gg_list[[1]], gg_list[[2]], gg_list[[3]], cols=3)
dev.off()
