########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##           Fig. 4, Fig. S6          ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
# 
# the script lets asexual populations evolve until segregation-selection limit
# and for the resulting stable population structure, allows to 1) plot the pop
# composition (as in Fig. 4), and 2) calculate the fitness advantage the offspring
# of a fuser would get, compared to the offspring of a non-fuser.
# The script can generate Fig. S6, that is the initial advantage of fusion,
# for a range of mitochondria number
# 

require(ggplot2)

source('PARAM_TABLE.R')

# the range of mito number to explore in Fig. S6
v.Mito = c(seq(2,20), seq(30, 200, 10))

cost_sex = 1

advantage = numeric(length(v.Mito))

plot_for <- c(2,10,50)
homoplasmy_list <- list()
advantage_list <- list()

## begin loop for each value of mitochondrial number explored
count = 1
for(run in 1:length(v.Mito)){ 
  
  Mito = v.Mito[run]
  
  if(josh){
    W <- numeric(Mito+1)
    selection_coefficient <- Ka1
    expon <- 2
    
    counter <- (Mito/2+1)
    for (ii in 1:(Mito/2+1)){
      W[counter] <- 1 - selection_coefficient * (ii/(Mito/2+1))^expon
      counter <- counter - 1
    }
    
    for (ii in 1:(Mito/2+1)){
      W[ii+ Mito/2 ] <- 1 - selection_coefficient * (ii/(Mito/2+1))^expon
    }
    fitness.diagm <- W * diag(1, Mito+1, Mito+1)
  } else {
    
    fitness.function <- function(x) 
      ( Ka1 /( 1 + 1*exp( -B1*(x-M1) )) + bsline1 ) * 
      (Ka2 /( 1 + 1*exp( B2*(x-(1-M2)) )) + bsline2 ) 
    fitness.diagm <- fitness.function( (0:Mito)/Mito ) * diag(1, Mito+1, Mito+1)
  }
  
  # make matrix of transition from the composition of a pop of non-fusers
  # to a pop formed of their offspring (thx Josh for code snippet)
  
  non_fuser_transition_matrix_with_replacement <- function(Mito) {
    transition <- matrix(numeric( (Mito + 1) * 2 ), Mito + 1, Mito + 1)
    for (ii in 0 : Mito) {
      for (jj in 0 : Mito) {
        transition[ii + 1, jj + 1] <- dbinom(jj, Mito, ii / Mito)
      }
    }
    return(transition)
  }
  non_fuser_transition_matrix <- function(Mito) {
    transition <- matrix(numeric( (Mito + 1) ^ 2 ), Mito + 1, Mito + 1)
    for (ii in 0 : Mito) {
      for (jj in 0 : Mito) {
        transition[ii + 1, jj + 1] <- dhyper(jj,
                                             2 * ii,
                                             (2 * Mito) - (2 * ii),
                                             Mito)
      }
    }
    return(transition)
  }
  

  ## does fission happens with sampling with or without replacement
  
  if(replacement){
    transition <- non_fuser_transition_matrix_with_replacement( Mito ) 
  } else {
    transition <- non_fuser_transition_matrix( Mito ) 
  }
  
  # make matrix of transition from the composition of the pop after fusion
  # (big cells), to the pop after the ensuing fission
  
  fuser_transition_matrix <- function(Mito) {
    transition <- matrix( numeric(), ( 2*Mito + 1 ), Mito + 1 )
    for (ii in 0 : (2*Mito)) {
      for (jj in 0 : Mito) {
        transition[ii + 1, jj + 1] <- dhyper( jj, ii, (2 * Mito) - ii, Mito )
      }
    }
    return(transition)
  }
  transition2 <- fuser_transition_matrix( Mito )
  
  ## begin the simulation to reach population content equilibrium
  
  tpsBI = 1000 # max time to run if segregation-selection equilibrium not reached yet
  t = 1
  homoplasmy <- matrix(nrow = tpsBI, ncol = Mito+1 )
  homoplasmy_start <- numeric(Mito+1)
  homoplasmy_start[(Mito/2+1)] <- 1
  homoplasmy[t,] <- homoplasmy_start
  
  for(t in 1:tpsBI){
    
    # cells undergo fission
    non_fusers2 <- homoplasmy[t,] %*% transition
    
    # cells undergo selection based on their cytotypes
    non_fusers3 <- non_fusers2 %*% fitness.diagm
    
    # gathering results after scaling to sum up to 1
    homoplasmy[t+1,] <- (non_fusers3) / sum(non_fusers3)
    if(t>1){
      if( t/10 == round(t/10)) print(t)
      if( sum( abs( homoplasmy[t,] - homoplasmy[(t-1),] ) ) < 0.00001){
        break
      }
    }
  }
  # barplot(homoplasmy[t,])
  
  fitness.asex.m.type <- numeric() 
  fitness.sex.m.type <- numeric() 
  relative.fitness.m.type <- numeric()
  
  transition.fitness <- transition %*% fitness.diagm
  
  for(m in 1:(Mito+1)){
    
    # a cell with m mitochondria gives that distrib of daughters, 
    # with those fitnesses
    fitness.asex.m.type[m] <- sum( transition.fitness[m,] )
    
    fitness.prospects.m.type <- numeric()
    for( i in 1:(Mito+1)){ # a cell with m mitochondria fuses with a cell with i
      
      # daughter gets m + i mitos (-1 because of indexing with the 0 class)
      fitness.prospects.m.type[i] <- sum( 
        transition2[  (m+i-1), ] * rowSums(transition.fitness) ) 
    }
    
    # fitness prospect of each cross weighted by its probability
    fitness.sex.m.type[m] <- sum( fitness.prospects.m.type * homoplasmy[t,] ) * 
      cost_sex
    
    relative.fitness.m.type[m] <- fitness.sex.m.type[m] / fitness.asex.m.type[m]
  }
  
  
  relative.fitness <- sum( relative.fitness.m.type * homoplasmy[t,] )
  
  if(Mito %in% plot_for ){
    homoplasmy_list[[count]] <- homoplasmy[t,]
    advantage_list[[count]] <- relative.fitness.m.type
    count <- count + 1
  }
  
  advantage[run] = sum( relative.fitness.m.type * homoplasmy[t,] )

  }

## plot Fig. 4

par(mfrow=c(2,length(plot_for)))
for(i in 1:length(plot_for)){
  Mito = plot_for[i]
  names(homoplasmy_list[[i]]) <- 0:(Mito)/Mito
  
  barplot(homoplasmy_list[[i]], ylab="frequency",
          main = paste("M =",  plot_for[i]), xlab="% mitochondria of type A")
}
for(i in 1:length(plot_for)){
  Mito = plot_for[i]
  names(advantage_list[[i]]) <- 0:(Mito)/Mito
  
  barplot(advantage_list[[i]], ylim=c(0, 1.8), col="white", ylab="relative fitness",
          main = '', xlab="% mitochondria of type A")
  abline(h=1, col='red')
}

advantage[v.Mito%in%plot_for]




## plot Fig. S6

res = cbind( v.Mito, advantage )
colnames(res) = c("Mito", 'advantage')
res = na.omit( as.data.frame(res) )

p <- ggplot(res, aes(Mito, advantage))
p + geom_line() + theme_bw() + scale_y_continuous(limits = c(1, 1.125))+ 
  theme(aspect.ratio=4/3.5) +
  ggtitle("advantage of fusing if no cost")

ggsave("2. invasion potential.jpg")

res2 = res[res$Mito<=20,]
p <- ggplot(res2, aes(Mito, advantage))
p + geom_line() + theme_bw() + 
  theme(aspect.ratio=4/3.5) +
  ggtitle("advantage of fusing if no cost")

ggsave("2. zoom invasion potential.jpg")
