#####################################################################
##
##        Basic regression analysis with BCSaplingData
##
##        Using the likelihood package to run optim
##
##   NOTE:  this code assumes that you have already run the
##          code for "Basic regression with anneal.R"
##
####################################################################

#  The setup for the optim package can be cumbersome
#  If you are already used to setting up your models for annealing
#    using the likelihood package, you can use the "likeli_4_optim"
#    function to "wrap" the call to optim in terms you are already
#    familiar with...

#  The only real change in setup is to separate your list of parameters
#    into first a vector (not a list) with just initial values
par<-c(5, 0.1, 1, 1)

#  and then create another vector containing the names of the parameters
names <- c("a","s","d","sigma")

var <- list(X1 = "global", X2 = "size")

#  There isn't any real point in specifying lower and upper bounds on the
#    parameters, since the local optimization methods in optim generally don't
#    let you do bounded searches

#  par_lo<-list(a = -10, s = 0.00001, d = 0, sigma = 0)
#  par_hi<-list(a = 100, s = 500, d = 10, sigma = 10)

var$x<-"meanrg"

var$mean<-"predicted"

var$log<-TRUE

## Set your choice of optim controls - pass the other likeli_4_optim arguments
## by name so optim knows they are for likeli_4_optim
## Remember to set the fnscale option of optim to a negative value to perform
## a maximization rather than a minimization

optim(par, likeli_4_optim, method = "Nelder-Mead",
      control = list(fnscale = -1), model = size.MM.model, par_names = names,
      var = var, source_data = data, pdf = normal.with.power.variance)


#######################################################
#
#  optim with the simple normal pdf
#
######################################################

par<-c(5, 0.1, 1, 1)

#  and then create another vector containing the names of the parameters
names <- c("a","s","d","sd")

var <- list(X1 = "global", X2 = "size")

#  There isn't any real point in specifying lower and upper bounds on the
#    parameters, since the local optimization methods in optim generally don't
#    let you do bounded searches

#  par_lo<-list(a = -10, s = 0.00001, d = 0, sigma = 0)
#  par_hi<-list(a = 100, s = 500, d = 10, sigma = 10)

var$x<-"meanrg"

var$mean<-"predicted"

var$log<-TRUE

## Set your choice of optim controls - pass the other likeli_4_optim arguments
## by name so optim knows they are for likeli_4_optim
## Remember to set the fnscale option of optim to a negative value to perform
## a maximization rather than a minimization

optim(par, likeli_4_optim, method = "Nelder-Mead",
      control = list(fnscale = -1), model = size.MM.model, par_names = names,
      var = var, source_data = data, pdf = dnorm)

####################################################################
#
#  combining global and local optimization
#
####################################################################

#  For simple problems and datasets, optim will often work, but you should
#    always be aware of the effects of initial starting conditions

#  Since there are no formal methods for guaranteeing "convergence" of simulated
#    annealing (or any other global optimization routines, for that matter), one
#    option is to do simulated annealing first, then when you are confident you have
#    found the region of the global maximum likelihood, use optim to climb up to the
#    peak (i.e. run annealing, and pass the best parameter estimates to optim).

#  Remember that the object you created with your call to "anneal" (i.e. "results") is a list
#    with many elements:  one of them is a sub-list called "best_pars" that stores the
#    best estimate of the values of each of the parameter in "par".
#  So, pass these values as the starting values for the local optimization:
par<-c(results$best_pars$a, results$best_pars$s, results$best_pars$d,
       results$best_pars$sd)

#  and then create another vector containing the names of the parameters
names <- c("a","s","d","sd")

var <- list(X1 = "global", X2 = "size")
var$x<-"meanrg"
var$mean<-"predicted"
var$log<-TRUE

## Set your choice of optim controls - pass the other likeli_4_optim arguments
## by name so optim knows they are for likeli_4_optim
## Remember to set the fnscale option of optim to a negative value to perform
## a maximization rather than a minimization

optim(par, likeli_4_optim, method = "Nelder-Mead",
      control = list(fnscale = -1), model = size.MM.model, par_names = names,
      var = var, source_data = data, pdf = dnorm)