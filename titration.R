########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##         RUNNING SIMULATIONS        ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Script description
# -------------------------------------------------------
# To achieve more precision in the estimation of 
# - the highest rate of fusion able to reach fixation, and
# - the higher rate of fusion that does not go extinct, but reaches a polymorphic equilibrium 
# we use a "titration" procedure in 3 iterations.
# 
# Procedure for determination of the highest rate of fusion to fix:
#
# The max rate of fusion to fix in previous simulations (RESULT.txt) is the low starting bound.
# The min rate of fusion tested, that didn't fix, is the high starting bound.
# The threshold fusion rate value to achieve fixation is somewhere in between.
# First iteration runs a simulation with the average of the two starting bounds (log10-scale) as rate of fusion:
# - if it leads to fixation, this rate of fusion becomes the new low bound, and next iteration will use
# the average of this and the high starting bound.
# - if it doesn't lead to fixation, this rate of fusion becomes the new high bound, and next iteration will use
# the average of this and the low starting bound.
# Repeat.
#
# The procedure is similar for determination of the higher rate of fusion that does not go extinct


source('mitofunction.R')


# -------------------------------------------------------
# ########### PART 1
# ########### Estimate highest rate of fusion to fix (fixation threshold)
# -------------------------------------------------------
# 
.param <- read.table('RESULT.txt', h=T)

## recover parameter set for complementation function

v.cost_sex = unique(.param$cost_sex)[order(unique(.param$cost_sex))]
v.Mito = unique(.param$Mito)[order(unique(.param$Mito))]
v.rate_sex = unique(.param$rate_sex)[order(unique(.param$rate_sex))]
Ka1 <- .param$Ka1[1]
Ka2<- .param$Ka2[1]
bsline1 <- 1-Ka1
bsline2 <- 1-Ka2
B1 <- .param$B1[1]
B2 <- .param$B2[1]
M1 <- .param$M1[1]
M2 <- .param$M2[1]

## add boolean column: does that rate of fusion leads to fixation?

f_fix <- .param$f > .99
.param <- cbind(.param, f_fix)
.param$f_fix[ is.na(.param$f_fix) ] <- TRUE 

## prepare parameter / result table

estimate_higher_rate <- expand.grid(
  it = c("min0", "max0", "it1", "it2", "it3"), 
  cost_sex = v.cost_sex, 
  Mito = v.Mito)

## determine higher and lower bounds for fusion rate to start from, for each parameter set

for( i in 1:length(v.cost_sex) ){
  for( j in 1:length(v.Mito) ){
    df = subset( .param[ .param$cost_sex == v.cost_sex[ i ] & 
                           .param$Mito == v.Mito[ j ] , ] )
    estimate_higher_rate$rate_sex[ 
      estimate_higher_rate$it == "min0" &
        estimate_higher_rate$Mito == v.Mito[ j ] & 
        estimate_higher_rate$cost_sex == v.cost_sex[ i ] ] <- 
      max( df$rate_sex[ df$f_fix == TRUE ] )
  }
}



scale_rate_sex <- c(-Inf, v.rate_sex)

estimate_higher_rate$rate_sex[ estimate_higher_rate$it == "max0" ] <- 
  scale_rate_sex[ sapply( 
    estimate_higher_rate$rate_sex[ estimate_higher_rate$it == "min0"], 
    function(x) (which(scale_rate_sex == x)+1) ) ] 

estimate_higher_rate$lograte <- log10( estimate_higher_rate$rate_sex )

estimate_higher_rate$f <- numeric( length = length(estimate_higher_rate$it))
estimate_higher_rate$f[ estimate_higher_rate$it == "min0"] <- 1
estimate_higher_rate$f[ estimate_higher_rate$it == "max0"] <- 0

iterations <- c("it1", "it2", "it3")


for(i in 1:length(v.cost_sex)){
  for(j in 1:length(v.Mito)){
    # rm(list= c(minratelog, maxratelog, it, res, lograte_try))
    # gc()
    print(paste( count, "Mitos:", v.Mito[j], "cost sex:", v.cost_sex[i]))
    
    minratelog <- estimate_higher_rate$lograte[ 
      estimate_higher_rate$it == 'min0' &
        estimate_higher_rate$Mito == v.Mito[j] &
        estimate_higher_rate$cost_sex == v.cost_sex[i] ]
    
    maxratelog <- estimate_higher_rate$lograte[ 
      estimate_higher_rate$it == 'max0' &
        estimate_higher_rate$Mito == v.Mito[j] &
        estimate_higher_rate$cost_sex == v.cost_sex[i] ]
    
    if( !is.na( maxratelog) & !is.nan(minratelog)){  # if the min maximum rate to evolve is 1, no need to compute
      
      for(it in 1:3){
        lograte_try <- mean( c( minratelog, maxratelog) )
        
        res <- mitofunction(Ka1 = .param$Ka1[1],
                            Ka2 = .param$Ka2[1],
                            bsline1 = 1-Ka1,
                            bsline2 = 1-Ka2,
                            B1 = .param$B1[1],
                            B2 = .param$B2[1],
                            M1 = .param$M1[1],
                            M2 = .param$M2[1],
                            Mito = v.Mito[j],
                            cost_sex = v.cost_sex[i],
                            rate_sex = (10^lograte_try),
                            replacement = .replacement,
                            josh = .josh)
        
        estimate_higher_rate$lograte[ estimate_higher_rate$it == iterations[it] &
                                        estimate_higher_rate$Mito == v.Mito[j] &
                                        estimate_higher_rate$cost_sex == v.cost_sex[i] ] <- lograte_try
        estimate_higher_rate$f[ estimate_higher_rate$it == iterations[it] &
                                  estimate_higher_rate$Mito == v.Mito[j] &
                                  estimate_higher_rate$cost_sex == v.cost_sex[i] ] <- res[1]
        
        if( res[1] > 0.99){  minratelog <- lograte_try } else { maxratelog <- lograte_try }
        
      }
    }
  }
}


estimate_higher_rate$estimate <- round( 10^estimate_higher_rate$lograte, 4 )
write.table(estimate_higher_rate, "estimate_higher_rate.txt")

## extracts and treats results to be ready for plotting

results_titration <- matrix(nrow= length(v.cost_sex)*length(v.Mito), ncol=5)

count=1
for(i in 1:length(v.cost_sex)){
  for(j in 1:length(v.Mito)){
    
    results_titration[count, 1] <- v.Mito[j]
    results_titration[count, 2] <- v.cost_sex[i]
    results_titration[count, 3] <- max( estimate_higher_rate$lograte[ estimate_higher_rate$Mito == v.Mito[j] &
                                                                        estimate_higher_rate$cost_sex == v.cost_sex[i] &
                                                                        estimate_higher_rate$f > 0.99
                                                                      ])
    results_titration[count, 4] <- min( estimate_higher_rate$lograte[ estimate_higher_rate$Mito == v.Mito[j] &
                                                                        estimate_higher_rate$cost_sex == v.cost_sex[i] &
                                                                        estimate_higher_rate$f < 0.99
                                                                      ])
    count = count+1
  }
}
results_titration[!is.na(results_titration[,4]),5] <- 
  10^( (results_titration[!is.na(results_titration[,4]),3] + 
          results_titration[!is.na(results_titration[,4]),4] ) /2 )
results_titration[results_titration[,2] == 1,5] <- 1

colnames( results_titration ) <- c('Mito', "rate_sex", "min", "max", 'mean')

### 
dataset1 <- as.data.frame( results_titration[,c(2,1,5)] )
colnames(dataset1) <- c("fert_sex", "Mito", "result")
dataset1$result[ is.na(dataset1$result )] <- -Inf

write.table(dataset1, "titration_fix_df.txt", col.names = T, row.names = F)



# -------------------------------------------------------
# ########### PART 2
# ########### Estimate highest rate of fusion to reach polymorphism 
# ########### (polymorphism threshold)
# -------------------------------------------------------
# 


.param <- read.table('RESULT.txt', h=T)

v.cost_sex = unique(.param$cost_sex)[order(unique(.param$cost_sex))]
v.Mito = unique(.param$Mito)[order(unique(.param$Mito))]
v.rate_sex = unique(.param$rate_sex)[order(unique(.param$rate_sex))]
f_poly <- ( .param$f > 10^(-8) )
Ka1 <- .param$Ka1[1]
Ka2<- .param$Ka2[1]
bsline1 <- 1-Ka1
bsline2 <- 1-Ka2
B1 <- .param$B1[1]
B2 <- .param$B2[1]
M1 <- .param$M1[1]
M2 <- .param$M2[1]
.param <- cbind(.param, f_poly)
.param$f_poly[ is.na(.param$f_poly) ] <- FALSE # it has a true for f between 0.01 and 1 

estimate_lower_rate <- expand.grid(
  it = c("min0", "max0", "it1", "it2", "it3"), 
  cost_sex = v.cost_sex, 
  Mito = v.Mito)


# for each combination of mito and cost, range for the lowest rate that gives a polymorphism
estimate_lower_rate$f <- NA
for(i in 1:length( v.cost_sex )){
  for(j in 1:length( v.Mito )){
    df = subset( .param[ .param$cost_sex == v.cost_sex[i] & 
                           .param$Mito == v.Mito[j], ])
    estimate_lower_rate$rate_sex[ 
      estimate_lower_rate$it == "min0" &
        estimate_lower_rate$Mito == v.Mito[j] & 
        estimate_lower_rate$cost_sex == v.cost_sex[i] ] <- 
      max( df$rate_sex[ df$f_poly == TRUE ]) # returns Inf if there is no TRUE (no polymorphism) 
    
    # get the polymorphic frequency
    if( max( df$rate_sex[ df$f_poly == TRUE ]) != -Inf){
      estimate_lower_rate$f[ 
        estimate_lower_rate$it == "min0" &
          estimate_lower_rate$Mito == v.Mito[j] & 
          estimate_lower_rate$cost_sex == v.cost_sex[i] ] <- 
        df$f[ df$rate_sex == max( df$rate_sex[ df$f_poly == TRUE ])]
    }
  }
}

scale_rate_sex <- c(-Inf, v.rate_sex)

# upper limit of polymorphism (at that rate fusion was found to disappear)
estimate_lower_rate$rate_sex[ estimate_lower_rate$it == "max0" ] <- 
  scale_rate_sex[ sapply( 
    estimate_lower_rate$rate_sex[ estimate_lower_rate$it == "min0"], 
    function(x) (which(scale_rate_sex == x)+1) ) ] 

# dealing with cases where a rate of fusion of 1 is polymorphic
index_pb <- which( estimate_lower_rate$it == "min0" & estimate_lower_rate$rate_sex == 1)
estimate_lower_rate[index_pb+1,5 ] <- 1

estimate_lower_rate$lograte <- log10( estimate_lower_rate$rate_sex )

estimate_lower_rate$f[ estimate_lower_rate$it == "max0"] <- 0

iterations <- c("it1", "it2", "it3")
count <- 0

for(i in 1:length(v.cost_sex)){
  for(j in 1:length(v.Mito)){
    # rm(list= c(minratelog, maxratelog, it, res, lograte_try))
    # gc()
    
    minratelog <- estimate_lower_rate$lograte[ 
      estimate_lower_rate$it == 'min0' &
        estimate_lower_rate$Mito == v.Mito[j] &
        estimate_lower_rate$cost_sex == v.cost_sex[i] ]
    
    maxratelog <- estimate_lower_rate$lograte[ 
      estimate_lower_rate$it == 'max0' &
        estimate_lower_rate$Mito == v.Mito[j] &
        estimate_lower_rate$cost_sex == v.cost_sex[i] ]
    
    if( !is.na( maxratelog) & !is.nan(minratelog)){  # if the min maximum rate to evolve is 1, no need to compute
      
      for(it in 1:3){
        lograte_try <- mean( c( minratelog, maxratelog) )
        
        res <- mitofunction(Ka1 = .param$Ka1[1],
                            Ka2 = .param$Ka2[1],
                            bsline1 = 1-Ka1,
                            bsline2 = 1-Ka2,
                            B1 = .param$B1[1],
                            B2 = .param$B2[1],
                            M1 = .param$M1[1],
                            M2 = .param$M2[1],
                            Mito = v.Mito[j],
                            cost_sex = v.cost_sex[i],
                            rate_sex = (10^lograte_try), 
                            replacement = .replacement,
                            josh = .josh)
        
        estimate_lower_rate$lograte[ estimate_lower_rate$it == iterations[it] &
                                       estimate_lower_rate$Mito == v.Mito[j] &
                                       estimate_lower_rate$cost_sex == v.cost_sex[i] ] <- lograte_try
        estimate_lower_rate$f[ estimate_lower_rate$it == iterations[it] &
                                 estimate_lower_rate$Mito == v.Mito[j] &
                                 estimate_lower_rate$cost_sex == v.cost_sex[i] ] <- res[1]
        
        if( res[1] > 10^(-8)){  minratelog <- lograte_try } else { maxratelog <- lograte_try }
        print(paste( count, "Mitos:", v.Mito[j], "cost sex:", v.cost_sex[i]))
        count <- count + 1
      }
    }
  }
}

estimate_lower_rate$estimate <- round( 10^estimate_lower_rate$lograte, 4 )
write.table(estimate_lower_rate, "estimate_lower_rate.txt")


## extracts and treats results to be ready for plotting

results_titration <- matrix(nrow= length(v.cost_sex)*length(v.Mito), ncol=5)

count=1
for(i in 1:length(v.cost_sex)){
  for(j in 1:length(v.Mito)){
    
    results_titration[count, 1] <- v.Mito[j]
    results_titration[count, 2] <- v.cost_sex[i]
    results_titration[count, 3] <- max( na.omit( estimate_lower_rate$lograte[ estimate_lower_rate$Mito == v.Mito[j] &
                                                                                estimate_lower_rate$cost_sex == v.cost_sex[i] &
                                                                                estimate_lower_rate$f > 10^(-8)
                                                                              ]))
    results_titration[count, 4] <- min(na.omit( estimate_lower_rate$lograte[ estimate_lower_rate$Mito == v.Mito[j] &
                                                                               estimate_lower_rate$cost_sex == v.cost_sex[i] &
                                                                               estimate_lower_rate$f < 10^(-8)
                                                                             ]))
    # case when a rate of sex of 1 leads to polymorphism:
    if( estimate_lower_rate$rate_sex[ estimate_lower_rate$it == 'min0' &
                                      estimate_lower_rate$Mito == v.Mito[j] &
                                      estimate_lower_rate$cost_sex == v.cost_sex[i]] == 1 ){
      results_titration[count, 5] <- 1
    }
    count = count+1
  }
}
results_titration[!is.na(results_titration[,4]),5] <- 
  10^( (results_titration[!is.na(results_titration[,4]),3] + 
          results_titration[!is.na(results_titration[,4]),4] ) /2 )
results_titration[results_titration[,2] == 1,5] <- 1



colnames( results_titration ) <- c('Mito', "rate_sex", "min", "max", 'mean')

### 
dataset2 <- as.data.frame( results_titration[,c(2,1,5)] )
colnames(dataset2) <- c("fert_sex", "Mito", "result")
dataset2$result[ is.na(dataset2$result )] <- -Inf

write.table(dataset2, "titration_poly_df.txt", col.names = T, row.names = F)
