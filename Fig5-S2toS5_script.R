########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##      FITNESS CHANGE ANALYSIS       ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
#
# Compares the fitness mean and variance, and heteroplasmy mean and variance,
# of a population without fusion, to this population after having evolved fusion.
# The fusion rate is taken to be the maximum fusion rate found to be able to fix.
# The script generates the graphs from Fig. 5, Fig. S2-S5

library(R.utils)
library(ggplot2)
library(RColorBrewer)

source("param_table.R")
source("mitofunction.R")
source("asexual_fitness_function.R")

dataset1 <- read.table("titration_fix_df.txt", h=T)

subset <- dataset1[(dataset1$fert_sex*10)==round(dataset1$fert_sex/0.1) & 
           (dataset1$Mito/2)==round(dataset1$Mito/2) & dataset1$result !=(-Inf), ]


if(josh){
  fitness.function <- function( x, Ka1 = .3) 1-Ka1*((0.5-x)/0.5)^2
} else {
  fitness.function <- function(x) 
  ( Ka1 /( 1 + 1*exp( -B1*(x-M1) )) + bsline1 ) * 
  (Ka2 /( 1 + 1*exp( B2*(x-(1-M2)) )) + bsline2 ) 
}

mean_fitness_before <- numeric( dim(subset)[1] )
variance_fitness_before <- numeric(dim(subset)[1])
mean_fitness_after <- numeric(dim(subset)[1])
variance_fitness_after <- numeric(dim(subset)[1])

mean_mito_before <- numeric( dim(subset)[1] )
variance_mito_before <- numeric(dim(subset)[1])
mean_mito_after <- numeric(dim(subset)[1])
variance_mito_after <- numeric(dim(subset)[1])

# for each combination of parameters, the loop runs the evolutionary script and
# retains the final population composition (or if those runs already exist,
# comment out the function, and the script retrieves it from the datafiles)
# It also let a pop of fusion reach the segregation-selection equilibrium, and
# retains the final population composition

for(i in 1:dim(subset)[1]){

  print(paste(i, "/",dim(subset)[1] ))
  Mito = subset[i, 2]
  cost_sex = subset[i, 1]
  rate_sex = subset[i, 3]
  
  fitness.m.type <-  fitness.function( (0:Mito)/Mito )
  
  namerun = paste("run_M",Mito, "_cost", cost_sex,".RData", sep="" )

  # mitofunction(Ka1 <- Ka1,
  #              Ka2 <- Ka2,
  #              bsline1 <- 1-Ka1,
  #              bsline2 <- 1-Ka2,
  #              B1 <- B1,
  #              B2 <- B2,
  #              M1 <- M1,
  #              M2 <- M2,
  #              Mito <- Mito,
  #              cost_sex <- cost_sex,
  #              rate_sex <- rate_sex,
  #              josh <- josh,
  #              replacement <- replacement)


  homoplasmy_t <- loadToEnv(namerun)[["homoplasmy_t"]];
  mean_fitness_after[i] <- sum( fitness.m.type * homoplasmy_t )
  variance_fitness_after[i] <- sum( homoplasmy_t * fitness.m.type^2 ) - 
    mean_fitness_after[i]^2

  heteroplasmy <- c( seq(0,1, by=1/(Mito/2)), 
                     sort( seq(0,(1-1/(Mito/2)), by=1/(Mito/2)), 
                           decreasing = T ))
  
  mean_mito_after[i] <- sum( homoplasmy_t * heteroplasmy )
  variance_mito_after[i] <-  sum( homoplasmy_t * (heteroplasmy)^2 ) - 
    mean_mito_after[i]^2

  homoplasmy_asex <- asexual.fitness(Ka1 <- Ka1,
                  Ka2 <- Ka2,
                  bsline1 <- 1-Ka1,
                  bsline2 <- 1-Ka2,
                  B1 <- B1,
                  B2 <- B2,
                  M1 <- M1,
                  M2 <- M2,
                  Mito <- Mito,
                  josh <- josh,
                  replacement <- replacement)
  
  mean_fitness_before[i] <- sum( fitness.m.type * homoplasmy_asex )
  variance_fitness_before[i] <- sum( homoplasmy_asex * fitness.m.type^2 ) - 
    mean_fitness_before[i]^2
  
  mean_mito_before[i] <- sum( homoplasmy_asex * heteroplasmy )
  variance_mito_before[i] <-  sum( homoplasmy_asex * heteroplasmy^2 ) - 
    mean_mito_before[i]^2
}


dataset = data.frame(fert_sex = subset[,1], Mito = subset[,2], 
                     freq_sex = subset[,3], mean_before = mean_fitness_before, 
                     mean_after = mean_fitness_after, 
                     var_before = variance_fitness_before, 
                     var_after = variance_fitness_after, 
                     mean_mito_before = mean_mito_before, 
                     mean_mito_after = mean_mito_after, 
                     var_mito_before = variance_mito_before, 
                     var_mito_after = variance_mito_after)

variables = list(dataset$mean_before, dataset$mean_after, 
                 dataset$var_before, dataset$var_after,
                 dataset$mean_mito_before, dataset$mean_mito_after, 
                 dataset$var_mito_before, dataset$var_mito_after)
variable_names = c("mean_before", "mean_after", 
                   "var_before", "var_after",
                   "mean_mito_before", "mean_mito_after", 
                   "var_mito_before", "var_mito_after")
mypalette = list(brewer.pal(9, 'YlGnBu'), brewer.pal(9, 'YlGnBu'),
                 brewer.pal(9, 'YlOrRd'), brewer.pal(9, 'YlOrRd'),
                 brewer.pal(9, 'PuBu'), brewer.pal(9, 'PuBu'),
                 brewer.pal(9, 'PuRd'), brewer.pal(9, 'PuRd'))

mins = c(rep( min(dataset$mean_before, dataset$mean_after), 2),
         rep( min( dataset$var_before, dataset$var_after), 2),
         rep( 0, 2),
         rep( min( dataset$var_mito_before, dataset$var_mito_after), 2))
maxs = c(rep( 1, 2),
         rep( max( dataset$var_before, dataset$var_after), 2),
         rep( 1, 2),
         rep( max( dataset$var_mito_before, dataset$var_mito_after), 2))

write.table(dataset, 'RESULT_before_after.txt')

## plot FIG. S2

for(i in 1:2){
  if(i %% 2 == 0){
    text <- (variables[[i]] - variables[[i-1]] < 0)
    text[text==T] <- " "
    text[text==F] <- "+"
    dataset[,12] <-text
    names(dataset)[12] <- "text"
  } else {
    text <- rep('', length=length(variables[[i]]) )
  }
  print( ggplot(data=dataset,aes(x=fert_sex,y=as.factor(Mito))) +
  geom_tile(aes(fill=variables[[i]]), color='white')+
  scale_fill_gradientn(limits=c((mins[i]-0.1),maxs[i]),
                       colours=mypalette[[i]],
                       breaks=seq(round(mins[i], 2), round(maxs[i]), 0.1),
                       labels=seq(round(mins[i], 2), round(maxs[i]), 0.1),
                       na.value="white") +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     "relative reproductive success of fusers (1-c)") +
  scale_y_discrete(labels=sort(unique( dataset$Mito) ), "# mitochondria / cell") +
  labs(fill='Fitness mean') +
  theme(#panel.background = element_rect(fill = 'darkgoldenrod1'), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    aspect.ratio=4/3.5)+
    geom_hline(yintercept=10.5)+
    ggtitle(variable_names[i]) +
    geom_text(aes(x=fert_sex,y=as.factor(Mito)), label = text)
  )
  graph_title = paste("fitness_", variable_names[i], '.png', sep="")
  ggsave(graph_title)
}


## plot FIG. S4

for(i in 3:4){
  if(i %% 2 == 0){
    text <- (variables[[i]] - variables[[i-1]] < 0)
    text[text==T] <- " "
    text[text==F] <- "+"
    dataset[,12] <-text
    names(dataset)[12] <- "text"
  } else {
    text <- rep('', length=length(variables[[i]]) )
  }
  print( ggplot(data=dataset,aes(x=fert_sex,y=as.factor(Mito))) +
           geom_tile(aes(fill=variables[[i]]), color='white')+
           scale_fill_gradientn(limits=c(mins[i],maxs[i]),
                                colours=mypalette[[i]],
                                # breaks=seq(round(mins[i], 2), round(maxs[i]), 0.005),
                                # labels=seq(round(mins[i], 2), round(maxs[i]), 0.005),
                                na.value="white") +
           scale_x_continuous(breaks = seq(0, 1, 0.1),
                              "relative reproductive success of fusers (1-c)") +
           scale_y_discrete(labels=sort(unique( dataset$Mito) ), "# mitochondria / cell") +
           labs(fill='Fitness variance') +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             panel.background = element_blank(),
             panel.border = element_rect(fill = NA),
             aspect.ratio=4/3.5)+
           geom_hline(yintercept=10.5)+
           # ggtitle(variable_names[i]) +
           geom_text(aes(x=fert_sex,y=as.factor(Mito)), label = text)
  )
  graph_title = paste("fitness_", variable_names[i], '.png', sep="")
  ggsave(graph_title)
}

## plot FIG. S3

for(i in 5:6){
  
  if(i %% 2 == 0){
    text <- (variables[[i]] - variables[[i-1]] < 0)
    text[text==T] <- " "
    text[text==F] <- "+"
    dataset[,12] <-text
    names(dataset)[12] <- "text"
  } else {
    text <- rep('', length=length(variables[[i]]) )
  }
  print( ggplot(data=dataset,aes(x=fert_sex,y=as.factor(Mito))) +
           geom_tile(aes(fill=variables[[i]]), color='white')+
           scale_fill_gradientn(limits=c((mins[i]-0.1),maxs[i]),
                                colours=mypalette[[i]],
                                breaks=seq(round(mins[i], 2), round(maxs[i]), 0.2),
                                labels=seq(round(mins[i], 2), round(maxs[i]), 0.2),
                                na.value="white") +
           scale_x_continuous(breaks = seq(0, 1, 0.1),
                              "relative reproductive success of fusers (1-c)") +
           scale_y_discrete(labels=sort(unique( dataset$Mito) ), "# mitochondria / cell") +
           labs(fill='Mean prop type A mito') +
           theme(#panel.background = element_rect(fill = 'darkgoldenrod1'), 
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             panel.background = element_blank(),
             panel.border = element_rect(fill = NA),
             aspect.ratio=4/3.5)+
           geom_hline(yintercept=10.5)+
           ggtitle(variable_names[i]) +
           geom_text(aes(x=fert_sex,y=as.factor(Mito)), label = text)
  )
  graph_title = paste("heteroplasmy_", variable_names[i], '.png', sep="")
  ggsave(graph_title)
}

## plot FIG. S5

for(i in 7:8){
  
  if(i %% 2 == 0){
    text <- (variables[[i]] - variables[[i-1]] < 0)
    text[text==T] <- " "
    text[text==F] <- "+"
    dataset[,12] <-text
    names(dataset)[12] <- "text"
  } else {
    text <- rep('', length=length(variables[[i]]) )
  }
  
  print( ggplot(data=dataset,aes(x=fert_sex,y=as.factor(Mito))) +
           geom_tile(aes(fill=variables[[i]]), color='white')+
           scale_fill_gradientn(limits=c((mins[i]-0.1),maxs[i]),
                                colours=mypalette[[i]],
                                # breaks=seq(round(mins[i], 2), round(maxs[i]), 0.1),
                                # labels=seq(round(mins[i], 2), round(maxs[i]), 0.1),
                                na.value="white") +
           scale_x_continuous(breaks = seq(0, 1, 0.1),
                              "relative reproductive success of fusers (1-c)") +
           scale_y_discrete(labels=sort(unique( dataset$Mito) ), "# mitochondria / cell") +
           labs(fill='Heteroplasmy variance') +
           theme(#panel.background = element_rect(fill = 'darkgoldenrod1'), 
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             panel.background = element_blank(),
             panel.border = element_rect(fill = NA),
             aspect.ratio=4/3.5)+
           geom_hline(yintercept=10.5) +
           # ggtitle(variable_names[i]) +
           geom_text(aes(x=fert_sex,y=as.factor(Mito)), label = text)
  )
  graph_title = paste("heteroplasmy_", variable_names[i], '.png', sep="")
  ggsave(graph_title)
}
