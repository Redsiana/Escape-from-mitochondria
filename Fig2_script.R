########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##    PLOT Fig.2 fixation threshold   ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
#
# Produces Fig.2 displaying the maximum rate of sex to reach fixation (fixation 
# threshold), for each combination of mitochondria number and cost of fusion.
# Also produces equivalent graph for the maximum rate of sex to not go extinct
# (polymorphic thresholf)
#


library(plyr)
library(ggplot2)
library(RColorBrewer)



# -------------------------------------------------------
# ########### PART 1
# ########### Generate graph displaying the highest rate of fusion to fix
# -------------------------------------------------------
# 


dataset1 <- read.table("titration_fix_df.txt", h=T)

#### And now the plot

mypalette <- brewer.pal(9, 'PuBuGn')
newcol <- colorRampPalette(mypalette)
mypalette <- newcol(100)
# pie(rep(1, ncols), col = mypalette, border = NA, labels = NA)

c = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
b = log10(c)
v.Mito = sort(unique(dataset1$Mito))
p1 <- ggplot(data=dataset1,aes(x=fert_sex,y=as.factor(Mito)))+
  geom_tile(aes(fill=log10(result)), color='white')+
  scale_fill_gradientn(colours=mypalette,
                       breaks=b,
                       labels=c, 
                       na.value="white") +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     "relative reproductive success of fusers (1-c)") +
  scale_y_discrete(labels=v.Mito, "# mitochondria / cell") +
  labs(fill='Max. rate of\nfusion to fix') +
  theme(#panel.background = element_rect(fill = 'darkgoldenrod1'), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    aspect.ratio=4/3.5)+
  geom_hline(yintercept=10.5)
p1

ggsave('Fig2.pdf', device = "pdf", width = 166, height = 160, unit = "mm")





# -------------------------------------------------------
# ########### PART 2
# ########### Generate graph of displaying highest rate of fusion to not go extinct
# -------------------------------------------------------
# 

dataset2 <- read.table("titration_poly_df.txt", h=T)

#### And now the plot

mypalette <- brewer.pal(9, 'PuBuGn')
newcol <- colorRampPalette(mypalette)
mypalette <- newcol(100)
# pie(rep(1, ncols), col = mypalette, border = NA, labels = NA)

c = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
b = log10(c)

ggplot(data=dataset2,aes(x=fert_sex,y=as.factor(Mito)))+
  geom_tile(aes(fill=log10(result)), color='white')+
  scale_fill_gradientn(colours=mypalette,
                       breaks=b,
                       labels=c, 
                       na.value="white") +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     "relative reproductive success of fusers (1-c)") +
  scale_y_discrete(labels=v.Mito, "# mitochondria / cell") +
  labs(fill='Max. rate of\nfusion to reach\n polym eq') +
  theme(#panel.background = element_rect(fill = 'darkgoldenrod1'), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    aspect.ratio=4/3.5) +
  geom_hline(yintercept=10.5)

ggsave('polymorphic_threshold.pdf', device = "pdf", width = 166, height = 160, unit = "mm")

