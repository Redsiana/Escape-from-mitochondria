########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##     ASEXUAL EVOLUTION FUNCTION     ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
# 
# Starts a population of maximally heteroplasmic asexuals, let the population
# evolve until it reaches segregation-selection equilibrium
# Returns the composition of the population at segregation-selection equilibrium
# 

asexual.fitness <- function(Ka1,
                            Ka2,
                            bsline1,
                            bsline2,
                            B1 ,
                            B2 ,
                            M1 ,
                            M2 ,
                            Mito,
                            josh = F,
                            replacement){
  
  fitness.function <- function(x) 
    ( Ka1 /( 1 + 1*exp( -B1*(x-M1) )) + bsline1 ) * 
    (Ka2 /( 1 + 1*exp( B2*(x-(1-M2)) )) + bsline2 ) 
  
  fitness.diagm <- fitness.function( (0:Mito)/Mito ) * diag(1, Mito+1, Mito+1)
  
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
  } 
  
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
  
  if(replacement){
    transition <- non_fuser_transition_matrix_with_replacement( Mito )
  } else {
    transition <- non_fuser_transition_matrix( Mito ) 
  }
  
  fuser_transition_matrix <- function(Mito) {
    transition <- matrix( numeric( ), ( 2*Mito + 1 ), Mito + 1 )
    for (ii in 0 : (2*Mito)) {
      for (jj in 0 : Mito) {
        transition[ii + 1, jj + 1] <- dhyper(jj, ii, (2 * Mito) - ii, Mito)
      }
    }
    return(transition)
  }
  transition2 <- fuser_transition_matrix( Mito )
  
  
  tpsBI = 1000 # how long should it run if doesn't reach equilibrium before
  t = 1
  homoplasmy <- matrix(nrow = tpsBI+1, ncol = Mito+1 )
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
  return(homoplasmy[t,])
  
}