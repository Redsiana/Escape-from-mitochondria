########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##   PLOT COMPLEMENTATION FUNCTION    ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
# 
# Produces Fig. 1 plot of the complementation function (if josh = TRUE), 
# and Fig. S1 equivalents (if josh = FALSE)
# Only Ka1 matters for the function presented in Fig. 1 (main text), the other 
# parameters are needed for the function used in the supplementaries
# 

require(ggplot2)
source("PARAM_TABLE.R")


if(josh){
  # fitness.function_p2 = function(x, Ka1_=Ka1){
  #   1-Ka1_*((0.5-x)/0.5)^2
  # }
  
  Mito = 1000
  
  W <- numeric(100+1)
  selection_coefficient1 <- Ka1
  selection_coefficient2 <- Ka2
  expon <- 2
  
  counter <- (Mito/2+1)
  for (ii in 1:(Mito/2+1)){
    W[counter] <- 1 - selection_coefficient1 * (ii/(Mito/2+1))^expon
    counter <- counter - 1
  }
  
  
  for (ii in 1:(Mito/2+1)){
    W[ii+ Mito/2 ] <- 1 - selection_coefficient2 * (ii/(Mito/2+1))^expon
  }
  
  df <- data.frame( x = seq(0,1,0.001), y = W )
  ggplot(df) + 
    geom_line(aes(y=y, x=x)) +
    scale_y_continuous(limits = c(0,1), breaks = c(bsline1, bsline2, 1), "Viability")+
    scale_x_continuous(limits = c(0,1), "% mitotype A") +
    theme_bw()+
    theme(aspect.ratio=4/4, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
} else {
  fitness.function_p2 = function(x, Ka1_=Ka1, Ka2_=Ka2, bsline1_=bsline1, 
                                 bsline2_=bsline2, B1_=B1, B2_=B2, M1_=M1, M2_=M2){
    (Ka1_ /( 1 + 1*exp( -B1_*(x-M1_) )) + bsline1_) *
      ( (Ka2_ /( 1 + 1*exp( B2_*(x-(1-M2_)) )) + bsline2_))
  }
  
  ggplot(data.frame(x=c(0, 1)), aes(x)) +
    scale_y_continuous(limits = c(0,1), breaks = c(bsline1,1), "Viability")+
    scale_x_continuous(limits = c(0,1), "% mitotype A") +
    stat_function(fun = fitness.function_p2) +
    # ggtitle("Mitochondrial threshold function") +
    theme_bw()+
    theme(aspect.ratio=4/4, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


ggsave('Fig1.pdf', device = "pdf", width = 80, height = 80, unit = "mm")

