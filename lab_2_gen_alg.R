#################################################################
#
#       Basic regression analysis with BCSaplingData
#
#       Optimization using a genetic algoritm (RGENOUD)
#
#
#################################################################

# load the genetic algorithm library
library(rgenoud)

# load the BC Sapling Growth datafile
BCdata <- read.table("BC Sapling Growth Data.txt",header=T,sep="\t",as.is=T)

# display list of variables in the data file
str(BCdata)

# create a working dataset for one of the species
#  in one of the subzones
working.data <- subset(BCdata, spp == "HW" & subzone == "ICHMC2")

# drop any observations containing missing values for the independent or
#  dependent variables
#  NOTE:  by first "attaching" the dataframe, you don't have to preface
#    all of the variable names by the dataframe name...
attach(working.data)
data <- na.omit(data.frame(radius,global,meanrg))
detach(working.data)

data$size <- data$radius - 5*data$meanrg  #create a variable with size (radius
#   in mm) at the beginning of the measurement period

# examine scatterplots
plot(data$global,data$meanrg)
windows()
plot(data$size,data$meanrg)


#####################################################
#
#    Define the Function to be optimized (in GENOUD the 
#         function combines both the scientific model
#         (to calculate the predicted values) and the
#         likelihood function
#
#####################################################

#  NOTE:  the first parameter of the function must be a vector (not a list)
#    containing the vector of parameters to be optimized

# define a size-dependent Michaelis Menton model
size.MM.model <- function(a,s,d,x1,x2)
{  (a* x1*x2^d)/((a/s)+x1) }

# create a "wrapper" to hold the function to be optimized
ga.model <- function(par)
{  a <- par[1]  # need to assign the elements of par to the parameter names used within the function
   s <- par[2]
   d <- par[3]
   sd <- par[4]
   # now define the scientific model that calculates predicted values, given the parameters and variables
   pred <- size.MM.model(a,s,d,X1,X2)    # calculate predicted, given the parameters and the variables
   # now calculate (and return) the summed log likelihood
   return(sum(dnorm(x,pred,sd,log=T)))    # return summed log-likelihood
}



#######################################################
##
##  SETTING UP GENOUD
##
#########################################################


##  Define the vector (not a list) using the "c" function
##    for the parameters to be optimized
par<-c(a = 5, s = 1, d = 1, sd = 1)

## identify the independent variables needed by the function
X1 <- data$global
X2 <- data$size

## identify the dependent variable (observed)
x <- data$meanrg

## You can set bounds using the "domains" option in GENOUD
##   the bounds should be in a 2 column matrix, with the first column the
##   lower bound, and the order of rows should match the parameter vector

##   NOTE:  in the model used in the example, the s and sd parameters 
##     have to be greater than zero (algebraically)
bounds <- matrix(0,4,2)     # initialize the bounds matrix
bounds[,1]<-c(-10,0.00001,0,0.00001)
bounds[,2]<-c(100000,500,10,1000)

##  Now, call Genoud
result <- genoud(fn=ga.model,nvars=4,max=TRUE,pop.size=2000,
                 Domains=bounds,boundary.enforcement=2,hessian=T,print.level=1,
                 output.path = "GENOUD output.txt")

## now, display the results:
result

#  $value is the likelihood at the MLEs of the parameters

##  compute and display AIC
AIC <- -2*result$value + 2 * length(par)     # length(par) is the number of parameters
AIC

##  Other options to consider:
#
#   -- increase pop.size  (default is 1000)
#   -- make the domain limits strictly observed (boundary.enforcement = 2)
#   -- save output to a file (using "output.path")

#  The variance-covariance matrix is the inverse of the hessian, so use the "solve" function
#    to calculate it
solve(-1*result$hessian)

## the square roots of the diagonals of the inverted negative Hessian are the standard errors of the 
##   ML parameter estimates
sqrt(diag(solve(-1*result$hessian)))