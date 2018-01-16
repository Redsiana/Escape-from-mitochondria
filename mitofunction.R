########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##         SIMULATION FUNCTION        ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
# 
# Run this script to generate a compiled function that can run a simulation 
# for the specified parameters 
#
# The function returns:
# - the final frequency of the fusing allele in the pop (f_t)
# - and the time it took to reach it (t)
# It also saves the entire function environment in the working directory,
# including the distribution of cytoplasmic content at the last generation
# (heteroplasmy_t)
#
# A simulation can end in 4 ways: 
# - the fusing allele disappears ( f < 10^-8 )
# - the fusing allele fixes ( f > .99 )
# - the fusing allele reaches a polymorphic equilibrium ( f_t - f_tminus1 < 10^-7 for 1000 generations )
# - none of the above is reached, but more that 4 million generations have passed
# 

mitofunction <- function(Ka1,
                          Ka2,
                          bsline1,
                          bsline2,
                          B1 ,
                          B2 ,
                          M1 ,
                          M2 ,
                          Mito ,
                          cost_sex ,
                          rate_sex,
                          replacement = T,
                          josh = F) {

  tps <- 4000009                                                                  # a simulation runs for a max of tps generations
  
  # -------------------------------------------------------
  # ###########          I.
  # ########### Definition of functions
  # -------------------------------------------------------
  
  ## defining the complementation function
  
  fitness.function <- function(x) 
    ( Ka1 /( 1 + 1*exp( -B1*(x-M1) )) + bsline1 ) * 
    (Ka2 /( 1 + 1*exp( B2*(x-(1-M2)) )) + bsline2 )
  
  ## create matrix of transition from pre- to post-selection pop composition
  ## is a diagonal matrix, matching a cytotype with a fitness
  
  fitness.diagm <- fitness.function( (0:Mito)/Mito ) * diag(1, Mito+1, Mito+1)
  
  if(josh){
    W <- numeric(Mito+1)
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
    fitness.diagm <- W * diag(1, Mito+1, Mito+1)
  } 
  
  ## create matrix of transition from the composition of a pop of non-fusers
  ## to a the pop composed of their offspring (thx Josh for code)
  
  non_fuser_transition_matrix_with_replacement <- function(Mito) {                # mother -> daugther: sampling with replacement
    transition <- matrix(numeric( (Mito + 1) * 2 ), Mito + 1, Mito + 1)
    for (ii in 0 : Mito) {
      for (jj in 0 : Mito) {
        transition[ii + 1, jj + 1] <- dbinom(jj, Mito, ii / Mito)
      }
    }
    return(transition)
  }
  
  non_fuser_transition_matrix <- function(Mito) {                                 # mother -> daughter: sampling without replacement
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
  
  # -------------------------------------------------------
  # ###########          II.
  # ########### Define starting conditions 
  # ########### ACTION REQUIRED
  # -------------------------------------------------------
  

  ## CHOOSE MITOCHONDRIAL SAMPLING METHOD FOR ASEXUAL REPRODUCTION (FISSION) 
  ## comment/uncomment 1 of 2 lines below

  if(replacement){
    transition <- non_fuser_transition_matrix_with_replacement( Mito )              # WITH REPLACEMENT (MORE VARIANCE)
  }else{
    transition <- non_fuser_transition_matrix( Mito )                             # WITHOUT REPLACEMENT (LESS VARIANCE)
  }
  
  
  ## create matrix of transition from the composition of the pop after fusion 
  ## (big cells), to the pop after the ensuing fission
  
  fuser_transition_matrix <- function(Mito) {
    transition <- matrix(numeric( ), ( 2*Mito + 1 ), Mito + 1)
    for (ii in 0 : (2*Mito)) {
      for (jj in 0 : Mito) {
        transition[ii + 1, jj + 1] <- dhyper(jj, ii, (2 * Mito) - ii, Mito)
      }
    }
    return(transition)
  }
  
  transition2 <- fuser_transition_matrix( Mito )
  
  
  ## CHOSE WHETHER TO START WITH TWO HOMOPLASMIC LINES, OR A STEADY-STATE  
  ## DISTRIBUTION OF NON-FUSING, HETEROPLASMIC CELLS
  ## BY COMMENTING a) OR b)
  ## MODIFY c) FOR ORIGINAL CYTOPLASMIC BACKGROUND OF f MUTANT
  
  # ### a) SWITCH ON TO START AT STEADY STATE FOR A POP OF NON-FUSERS WITH MIXED CYTO (BURN-IT PHASE)
  # tpsBI = 200
  # t = 1
  # 
  # # homoplasmy <- matrix(nrow = tps+1, ncol = Mito+1 )
  # homoplasmy_start <- numeric(Mito+1)
  # homoplasmy_start[(Mito/2+1)] <- 1
  # homoplasmy_t <- homoplasmy_start
  # 
  # for(t in 1:tpsBI){
  #   
  #   # cells undergo fission
  #   non_fusers2 <- homoplasmy_t %*% transition 
  #   
  #   # cells undergo selection based on their cytotypes
  #   non_fusers3 <- non_fusers2 %*% fitness.diagm
  #   
  #   # gathering results after scaling to sum up to 1
  #   homoplasmy_tplus1 <- (non_fusers3) / sum(non_fusers3)
  #   
  #   if(t>1){
  #     if( t/10 == round(t/10)) print(t)
  #     if( sum( abs( homoplasmy_t - homoplasmy_tminus1 ) ) < 0.00001){
  #       break
  #     }
  #   }
  # }
  # 
  # homoplasmy_start <- homoplasmy_t
  
  
  #### b) SWITCH THOSE 3 LINES ON TO START WITH TWO HOMOPLASMIC LINES
  
  t = 1
  homoplasmy_start <- numeric(Mito+1)
  homoplasmy_start[c(1,(Mito+1))] <- c(0.5, 0.5)
  
  
  #### c) IN WHAT CYTOPLASMIC BACKGROUND DOES THE f ALLELE ORIGINATE?
  
  fuser_allele_start <- numeric(Mito+1)
  fuser_allele_start[1] <- 0.01                                                   # fusing allele starts at freq 0.01 in cytoplasmic class 1
  
  
  # -------------------------------------------------------
  # ###########          IV.
  # ########### Simulation 
  # -------------------------------------------------------
  
  ## define vector of fusing allele frequency among each cytotype at time t
  
  fuser_allele_t <- fuser_allele_start
  
  ## define vector of pop cytoplasmic composition at time t
  
  homoplasmy_t <- homoplasmy_start
  
  ## define total fraction of f allele in pop
  
  f_t <- as.numeric( fuser_allele_start %*% homoplasmy_start )#                   # at time t
  f_tminus2 <- numeric()                                                          # at time t-2
  f_tminus1 <- numeric()                                                          # at time t-1
  
  ## initialize count variable: increments when change in f_t and f_tminus1 is very small (to detect equilibrium)

  count <- 0

  
  # --- BEGIN LOOP --- #
  
  for( t in 1:tps ){
    
    ## cells carrying n (non-fusing) allele
    
    non_fusers <- (1-fuser_allele_t) * homoplasmy_t 
    
    ## cells carrying f (fusing) allele
    
    fusers <- fuser_allele_t * homoplasmy_t
    
    
    # --- SOME CELLS FUSE, SEPARATE, THEN DO MITOSIS --- #
    
    ## probability of each cross
    
    ff_x_ff <-  (fusers * rate_sex) %*% t( fusers * rate_sex )                    # fusing fuser x fusing fuser
    ff_x_nff <- (fusers * rate_sex) %*% t( fusers * (1-rate_sex) )                # fusing fuser x non-fusing fuser
    ff_x_nf <- (fusers * rate_sex) %*% t( non_fusers )                            # fusing fuser x non-fuser
    
    ## cytotype distribution of double cells generated by fusion for each type of cross
    
    fused_cells_ff_x_ff <- numeric( 2*Mito+1 )                                    # fusing fuser x fusing fuser
    
    for(i in 1:(Mito+1)){
      for(j in 1:(Mito+1)){
        fused_cells_ff_x_ff[i+j-1] <- fused_cells_ff_x_ff[i+j-1] + ff_x_ff[i,j]   # probability of inheriting (i+j-1) mitos from a ff x ff cross
      }
    }
    
    fused_cells_ff_x_nf <- numeric( 2*Mito+1 )                                    # fusing fuser x non-fuser
    
    for(i in 1:(Mito+1)){
      for(j in 1:(Mito+1)){
        fused_cells_ff_x_nf[i+j-1] <- fused_cells_ff_x_nf[i+j-1] +  ff_x_nf[i,j]
      }
    }
    
    fused_cells_ff_x_nff <- numeric( 2*Mito+1 )                                   # fusing fuser x non-fusing fuser
    
    for(i in 1:(Mito+1)){
      for(j in 1:(Mito+1)){
        fused_cells_ff_x_nff[i+j-1] <- fused_cells_ff_x_nff[i+j-1] +  ff_x_nff[i,j]
      }
    }
    
    ### double then cells separate into single cells, THEN undergo normal mitosis
    ### because that takes time, their contrib to next generation is cost_sex, not 1
    
    offspring_f_x_f <- fused_cells_ff_x_ff %*%                                    # offspring of fusing fuser x fusing fuser crosses
      transition2  %*% 
      transition * 
      cost_sex
    offspring_f_x_nf <- fused_cells_ff_x_nf %*%                                   # offspring of fusing fuser x non-fuser crosses
      transition2  %*% 
      transition * 
      cost_sex
    offspring_f_x_nff <- fused_cells_ff_x_nff %*%                                 # offspring of fusing fuser x non-fusing fuser crosses
      transition2  %*% 
      transition * 
      cost_sex

    # --- SOME CELLS HAVE ESCAPED FUSION AND UNDERWENT MITOSIS RIGHT AWAY --- #
    
    offspring_nf <- ( non_fusers * (1-f_t*rate_sex) ) %*% transition              # cells with the n allele ('non-fusers') that didn't get grabbed by a fuser
    offspring_nff <- ( fusers*(1-rate_sex) * (1-f_t*rate_sex) ) %*% transition    # cells with f allele, didn't want to fuse this generation ('non-fusing fusers'), didn't get grabbed 
    
    
    # --- POPULATION POST-FUSION, POST-MITOSIS, PRE-SELECTION --- #
    
    fusers2 <- offspring_f_x_f +                                                  # cytotype distribution among fusers after repro, before selection
      offspring_f_x_nf + 
      2*offspring_f_x_nff + 
      offspring_nff
    non_fusers2 <- offspring_nf + offspring_f_x_nf                                # cytotype distribution among non-fusers after repro, before selection

    
    # --- SELECTION ROUND COMES HERE --- #
    
    fusers3 <- fusers2 %*% fitness.diagm                                          # cytotype distribution among fusers after repro and selection
    non_fusers3 <- non_fusers2 %*% fitness.diagm                                  # cytotype distribution among non-fusers after repro and selection
    
    
    # --- GATHERING POP COMPOSITION, SCALING SO THAT IT SUMS TO 1 --- #
    
    homoplasmy_tplus1 <- as.vector(                                               # cytotype distribution in the population of the new generation
      (non_fusers3 + fusers3) / sum(non_fusers3 + fusers3) )
    fuser_allele_tplus1 <- as.vector(                                             # frequency of the fusing allele in each cytoplasmic class of the new generation
      (fusers3 / sum(non_fusers3 + fusers3)) / homoplasmy_tplus1 )
    f_tplus1 <- as.vector( fuser_allele_tplus1 %*% homoplasmy_tplus1 )            # frequency of the fusing allele in the total population of the new generation
    
    
    # --- WRAPPING UP THE GENERATION --- #
    
    ## print the number of generations elapsed (every 1000 generation)
    
    if(t>1){
      if( t/1000 == round(t/1000)) print(t)
    }
    
    ## stop simulation if one mitotype has gone to fixation
    
    if(homoplasmy_t[1]>0.999 | homoplasmy_t[Mito+1]>0.999 ) break
    
    ## stop simulation if the fusing allele disappears
    
    if(f_t < 10^(-8)) break
    
    ## to detect if the frequency of the fusing allele is reaching an equilibrium
    
    if(t > 3){
      delta_over_2g <- 
        na.omit( abs( (f_tminus1-f_tminus2) /f_tminus1 ) < 10^-7 ) *              # is the change in freq below 10^-7 for at least 2 generations?
        na.omit( abs( (f_t-f_tminus1) /f_t ) < 10^-7 )
      
      if( delta_over_2g == 0 ){ count <- 0 } else { count <- count + 1 }          # if yes, increment count; if no, count returns to zero
    }
    
    ## initialize next generation
    
    homoplasmy_tminus1 <- homoplasmy_t
    homoplasmy_t <- homoplasmy_tplus1
    homoplasmy_tplus1 <- numeric()
    
    fuser_allele_t <- fuser_allele_tplus1
    fuser_allele_tplus1 <- numeric()
    
    f_tminus2 <- f_tminus1
    f_tminus1 <- f_t
    f_t <- f_tplus1
    f_tplus1 <- numeric()
    
    ## if freq of f did not change by more that 10^-7 for 1000 generations, stop simulation
    
    if( count > 1000 ) break  
  }
  
  # -------------------------------------------------------
  # ###########          V.
  # ########### The function returns data
  # -------------------------------------------------------
  
  ## return final freq of fusing allele, and the number of generations needed to reach it
  namerun = paste("run_M",Mito, "_cost", cost_sex,".RData", sep="" )
  save( list = ls( all.names = TRUE ), file = namerun, envir = environment() )
  return(c( f_t, t, sum(homoplasmy_t[1], homoplasmy_t[Mito+1])))

}


## compile the function for faster execution
library(compiler)
mitofunction <- cmpfun(mitofunction)

#.