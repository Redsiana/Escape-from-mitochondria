########################################
##                                    ##
##      EVOLUTION OF CELL FUSION      ##         _~
##           PARAMETERS USED          ##     _~ )_)_~
##          last edit: 01/18          ##    )_))_))_)
##        author: Anais Tilquin       ##    _!__!__!_ /
##          redsiana@orange.fr        ##    \______t/
##                                    ##  ~~~~~~~~~~~~~
########################################

# -------------------------------------------------------
# ########### Creates the table of parameters to run simulations with
# -------------------------------------------------------

## The fitness function parameters are set constant, but the rate of fusion
## coded by the fusing allele, the cost of fusion, and the number of mitochondria
## per cell, are varied - a table is created with all combinations to run
## simulations for

## Fitness function parameters:

## Which function to use? josh = T is the main paper example (Eqn 1), 
## josh = F the one in ## the supplementary material 

josh  <- TRUE

## fitness difference btw fully homo- and fully heteroplasmic cells
## (Ka1 is the only parameter taken by the main text complementation function, Eqn1)
Ka1 <- 0.3
Ka2 <- 0.3
bsline1 <- 1-Ka1
bsline2 <- 1-Ka2
## shape parameter for the threshold steepness (flatness of the hat)
B1 <- 25
B2 <- 25
## shape parameter for the threshold midpoint (width of the hat)
M1 <- 0.2
M2 <- 0.2

## the numbers of mitochondria per cell to explore
v.Mito = c(seq(2,20,2), 50, 160, 200)

## the rates at which a mutant cell attempts fusion to explore
v.rate_sex <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)

## the costs of fusion to explore, between 0 and 1
## (engaging in fusion multiplies offspring number by cost_sex)
v.cost_sex <- seq(.05, 1, .05)

## create parameter table .param

.param <- expand.grid(Ka1, Ka2, B1, B2, M1, M2, v.Mito, v.rate_sex, v.cost_sex)
.param <- cbind(.param, rep( 'NA', dim(.param)[1] ), rep( 'NA', dim(.param)[1] ), rep( 'NA', dim(.param)[1] ))
names(.param) <- c('Ka1', 'Ka2', 'B1', 'B2', 'M1', 'M2', 
                   'Mito', 'rate_sex', 'cost_sex', 'f', 'time', 'h')
.param$f <- as.vector(.param$f)                                                   # final freq of fusing allele will be stored here
.param$time <- as.vector(.param$time)                                             # time to reach this freq will be stored here
.param$h <- as.vector(.param$h)
.param <- .param[with(.param, order(Mito, cost_sex, -rate_sex)), ]
