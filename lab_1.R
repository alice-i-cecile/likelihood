#################################################################
#
#                     Likelihood Methods in Ecology
#
#                      Lab 1 - MLE and Regression
#
#################################################################


#################################################################
#
#   Section 1:  Grid search for MLEs of the mean and
#      variance of a sample
#
##################################################################

# Libraries
library(ggplot2)

#  generate a sample of observations with a normal distribution

mean <- 15
sd <- 2.5
n <- 100

x <- rnorm(n,mean,sd)    # check the R help for "rnorm" if you're not familiar with the function

#  check the traditional (method of moments) estimates of the mean and standard deviation
mean(x)
sd(x)

#   You have to have a PDF to calculate likelihood, so take a first look at your data to make an
#     initial choice of the PDF
ggplot(data=data.frame(x=x), aes(x=x)) + geom_histogram()

#   Now, set up a grid search to compute the likelihoods of different combinations of estimates of
#     the mean and variance.

#   define the dimensions of the grid
range <- 0.25
increment <- 0.01

# these next lines just create two vectors of increments to use for the grid search of 
#   the mean and variance, using the "seq" (sequence) function in R - they will search plus and
#   minus 25% of the values for mean and sd you specified above
muhat <- seq(mean-(range*mean),mean+(range*mean),by=mean*increment)
sdhat <- seq(sd-(range*sd),sd+(range*sd),by=sd*increment)

# Make a data.frame for all possible combinations of these choices
search_space <- expand.grid(muhat, sdhat)
names(search_space) <- c("mu", "sd")


#   It's not called the normal distribution for nothing...
#     Given the histogram of observed data, let's assume the data are normally distributed,
#     and use dnorm (the R function for the normal PDF)
#     NOTE: remember that ultimately it's the distribution of residuals that matters, not the
#     distribution of the observations!


# Create a likelihood function for the normal function
ll_norm <- function(x, muhat, sdhat)
{
  log_likelihood <- sum(dnorm(x,muhat,sdhat,log=T))
  return (log_likelihood)
} 

# Evaluate the log-likelihood at each data point
ll_df <- search_space

ll_df$ll <- mapply(ll_norm, muhat=search_space$mu, sdhat=search_space$sd, MoreArgs=list(x=x))


#  what is the maximum likelihood calculated for the set of points
max_lh <- max(ll_df$ll)
max_lh

#  Now, display a contour plot of the summed log likelihoods
ggplot(data=ll_df, aes(x=mu, y=sd, z=ll)) + geom_contour()

# Or use a heat map
ggplot(data=ll_df, aes(x=mu, y=sd, fill=ll)) + geom_tile() + scale_fill_gradient(low="red", high="white")

# Find which row has the greatest log-likelihood
max_row <- which.max(ll_df$ll)

# check that you actually did find the right maximum
# maximum likelihood estimate of the mean (muhat)
# maximum likelihood estimate of the standard deviation
print (ll_df[max_row, ])

#  Now, calculate residuals (errors), and see if they are indeed normally distributed with a mean of zero
#   and a standard deviation equal to the MLE
res <- x - ll_df[max_row, "mu"]

res_plot <- ggplot(data=data.frame(x=res), aes(x=x)) + geom_histogram()
print(res_plot)

#  We can superimpose a normal PDF over this, given the MLEs of the
pdfx <- seq(from = min(res), to = max(res), by = 0.5)
pdfy <- dnorm(pdfx,0,ll_df[max_row, "sd"])
pdf_df <- data.frame(x=pdfx, y=pdfy)

# Rescale the y values
plot_info <- ggplot_build(res_plot)
ymean <- mean(plot_info$data[[1]]$count)

pdf_df$y <- pdf_df$y/mean(pdf_df$y)*ymean


# Plot the histogram and rescaled normal curve
res_plot + geom_line(data=pdf_df, aes(x=x, y=y))


points(pdfx,pdfy,type="p",pch=19,col="Blue",cex=2)
lines(pdfx,pdfy,type="l",col="Blue",lwd=2)

# Or use a quantile-quantile normal plot
qqnorm(res)

#  "Profile" likelihoods are plots of the change in likelihood as a particularl parameter changes
# (but with other parameters held at their MLE

#  We'll cover "support intervals" in a later lecture, but one way of summarizing the information in
#    the profile is to look at the range of parameter estimates that is within a certain number of units of
#   log likelihood (this is known as "support").  As you'll see, the range of parameter values within 2 units of the
#    maximum likelihood is roughly equivalent to a 95% confidence interval (although a true likelihoodist shudders
#    at that analogy).

#  Superimpose a line on the plot that is 2 units down (less than) the maximum log likelihood

#  Profile likelihood for estimate of the mean
max_sd_df <- ll_df[ll_df$sd==max(ll_df$sd[which.max(ll_df$ll)]),]

mean_profile <- ggplot(data=max_sd_df, aes(x=mu, y=ll)) + geom_line() + geom_abline(intercept=max(ll_df$ll)-2, slope=0)
print(mean_profile)


#  Profile likelihood for estimate of the standard deviation
max_mean_df <- ll_df[ll_df$mu==max(ll_df$mu[which.max(ll_df$ll)]),]

sd_profile <- ggplot(data=max_mean_df, aes(x=sd, y=ll)) + geom_line() + geom_abline(intercept=max(ll_df$ll)-2, slope=0)
print(sd_profile)