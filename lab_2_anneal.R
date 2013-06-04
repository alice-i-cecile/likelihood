#################################################################
#
#       Basic regression analysis with BCSaplingData
#
#       I  - Basic Michaelis Menton model with size dependency
#
#################################################################

# load the simulated annealing library
library(likelihood)

setwd("C:\\Users\\canhamc\\Documents\\Likelihood Course\\Course_Schedule_2012\\Day_1")
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

#######################################################
#
#    PDFs
#
######################################################

my.dnorm <- function(x,mean,sd) {dnorm(x,mean,sd,log=T)}

# define an alternative to the normal PDF in which the variance is
#   a power function of the mean
normal.with.power.variance <- function(x,mean,sigma)
{  sd <- mean^sigma
   dnorm(x,mean,sd,log=T)
}

#####################################################
#
#    Alternate scientific models
#
#####################################################

# define a basic Michaelis Menton model
MM.model <- function(a,s,X1)
{ (a* X1)/((a/s)+X1) }

# define a size-dependent Michaelis Menton model
size.MM.model <- function(a,s,d,X1,X2)
{  (a* X1*X2^d)/((a/s)+X1) }

# define a basic linear model
linear.model <- function(a,b,X1) { a + b*X1 }

# define a basic exponential regression model
exp.model <- function(a,b,X1) {  a + exp(b*X1) }

# define a basic power function model
power.model <- function(a,b,X1) { a*X1^b }

#######################################################
##
##  SETTING UP ANNEAL
##    The "anneal" function in the likelihood package requires:
##    1.  a scientific "model"
##    2.  a pdf to use to calculate likelihood
##    3.  a list of the parameters to estimate (par), with initial starting values
##    4.  a list of the variables in the model (var)
##
##    In addition, you can define the upper and lower limits for the
##      range to search for each parameter (par_lo and par_hi)
##
#########################################################

##  Set up par and var for the size-dependent MM model
##    NOTE:  par will include both the parameters of the scientific model
##    plus any parameters needed for the PDF...
par<-list(a = 5, s = 1, d = 1, sd = 1)

var <- list(X1 = "global", X2 = "size")       # you're telling anneal what variables in the dataset
#   to plug in as X1 and X2 in the scientific model

## Set bounds and initial search ranges within which to search for parameters
##   NOTE:  the s and sd parameters have to be greater than zero (algebraically)
par_lo<-list(a = -10, s = 0.00001, d = 0, sd = 0.00001)
par_hi<-list(a = 100000, s = 500, d = 10, sd = 1000)

## Specify the dependent variable, and name it using whatever
##    argument name is used in the pdf ("x" in the normal.pdf function above)
var$x<-"meanrg"

## predicted value in your PDF should be given the reserved name "predicted"
##   Anneal will use your scientific model and your data and calculate the
##   predicted value for each observation and store it internally in an object
##   called "predicted"
var$mean<-"predicted"

## Have it calculate log likelihood
var$log<-TRUE

##  now call the annealing algorithm, choosing which model to use
results<-anneal(size.MM.model,par,var,data,par_lo,par_hi,dnorm,"meanrg",
                hessian = F, max_iter=5000)

## now write the results to a file...
write_results(results,"Absolute MM model for BL.txt")

## It is also useful to save the "results" object (called "results" in this example, but you can name it anything you want)
##   -- If you save it with a filetype ".Rdata" it will be clear (to you and the operating system) that this file contains
##   an R data object
##   -- Once saved, you can reload it later if you need to do more analysis of the results (looking at residuals, etc.)
save(results,file="Absolute MM model for BL.Rdata")

## display some of the results in the console
results$best_pars;
results$max_likeli;
results$aic_corr ;
results$slope;
results$R2

#################################################################
##
##  using an alternate PDF (variance proportional to the mean)
##
################################################################

par<-list(a = 5, s = 0.1, d = 1, sigma = 1)

var <- list(X1 = "global", X2 = "size")

par_lo<-list(a = -10, s = 0.00001, d = 0, sigma = 0)
par_hi<-list(a = 100, s = 500, d = 10, sigma = 10)

var$x<-"meanrg"

var$mean<-"predicted"

var$log<-TRUE

next.results<-anneal(size.MM.model,par,var,data,par_lo,par_hi,normal.with.power.variance,"meanrg",max_iter=5000,
                     hessian = F)

## now write the results to a file...
write_results(next.results,"Absolute MM model for HW with heteroscedasticity.txt")

## and save the Rdata object
save(next.results,file="Absolute MM model for HW with heterscedasticity.Rdata")

## display some of the results in the console
next.results$best_pars;
next.results$max_likeli;
next.results$aic_corr ;
next.results$slope;
next.results$R2
