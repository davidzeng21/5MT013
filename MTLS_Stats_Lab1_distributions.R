# Lab 1: Distributions and sampling

# For this lab session we’ll be simulating data using R. We strongly advise
# you to attempt “Lab 2: Intro to R” before you attempt this material, unless
# you are an experienced R user.

# All the code you need is included in the R script file
# “MTLS_Stats_Lab2_distributions.R” available on Canvas. We encourage you to
# modify the code and discuss the results with your colleagues.

# load useful R packages
library(ggplot2)

# set the random seed to any number to get the same sequence of random numbers
# every time you run this script
set.seed(2378910)

# 1.
# a) Using R, produce a random sample of 10 000 from a Normal distribution with
# mean 20 and standard deviation 4. You should change the seed before drawing
# your sample. Plot a histogram and convince yourself that the data are indeed
# normally distributed.
sample1 <- data.frame(x=rnorm(10000, mean = 20, sd = 4))
ggplot(sample1) + geom_histogram(aes(x), bins = 30)
#
# b) Now draw a random sample of 10 from the same Normal distribution, and again
# plot a histogram.
sample2 <- data.frame(x=rnorm(10, mean = 20, sd = 4))
ggplot(sample2) + geom_histogram(aes(x), bins = 30)
#
# c) Without knowing how the data were generated, what would you be prepared to
# say about the underlying distribution of the data? Discuss with your neighbor
# – compare your histogram with theirs. Are they similar? Would they lead you to
# the same conclusion about the underlying population?

# 2.
# a) Using R again, draw a random sample of 10 from a Bernoulli distribution
# with p=0.5
# (note that the Bernoulli distribution is a special case if the binomial
# distribution, so we use an R function for the latter)
sample3 <- rbinom(10, 1, 0.5)

#
# b) What proportion of successes do you have in your sample? How far away is
# this from the true population proportion of 0.5?
mean(sample3)

#
# c) Imagine you didn’t know the true proportion – based on your survey alone,
# what range of values around your estimate do you think would capture the true
# proportion. Compare your answer with some colleagues – how much do the
# estimates vary?
#
# d) Repeat the above, but using a sample size of 100, and then 1000.
sample4 <- rbinom(100, 1, 0.5)
mean(sample4)

sample5 <- rbinom(1000, 1, 0.5)
mean(sample5)

# 3. Repeat question 2, but with a true proportion for the Bernoulli
# distribution of 0.96. What do you conclude?
mean(rbinom(10, 1, 0.96))
mean(rbinom(100, 1, 0.96))
mean(rbinom(1000, 1, 0.96))

# 4.
# a) Draw a sample of 1000 from a Poisson distribution with a mean of 1. Plot the data.
sample6 <- data.frame(x = rpois(1000, lambda = 1))
ggplot(sample6) + geom_bar(aes(x))

#
# b) Repeat with a mean of 5, then 10 and 100. What happens to the shape of the
# distribution as the mean increases?
ggplot(data.frame(x = rpois(1000, lambda = 5))) + geom_bar(aes(x))
ggplot(data.frame(x = rpois(1000, lambda = 10))) + geom_bar(aes(x))
ggplot(data.frame(x = rpois(1000, lambda = 100))) + geom_bar(aes(x)) + geom_histogram(aes(x), binwidth = 4)
#
# c) What happens if the sample size is much smaller, say 10?
ggplot(data.frame(x = rpois(10, lambda = 5))) + geom_bar(aes(x))

#5.What happens to the negative binomial distribution NegBin(n, p) as p changes?
#First, try N=1 and investigate p=0.2, 04, 0.6, 0.8
#Then repeat, but with N=2, N=10, N=100?
#What do you find? Discuss with your colleagues.
ggplot(data.frame(x = rnbinom(1000, size=1, prob = 0.2))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=1, prob = 0.4))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=1, prob = 0.6))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=1, prob = 0.8))) + geom_bar(aes(x))

ggplot(data.frame(x = rnbinom(1000, size=2, prob = 0.2))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=2, prob = 0.4))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=2, prob = 0.6))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=2, prob = 0.8))) + geom_bar(aes(x))

ggplot(data.frame(x = rnbinom(1000, size=10, prob = 0.2))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=10, prob = 0.4))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=10, prob = 0.6))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=10, prob = 0.8))) + geom_bar(aes(x))

ggplot(data.frame(x = rnbinom(1000, size=100, prob = 0.2))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=100, prob = 0.4))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=100, prob = 0.6))) + geom_bar(aes(x))
ggplot(data.frame(x = rnbinom(1000, size=100, prob = 0.8))) + geom_bar(aes(x))

# 6. (Advanced)

# There are lots of other distributions. Using R, simulate data to investigate
# some other distributions, and how they change shape as the parameters increase
# or decrease. Here are some you might like to try:

# • Student’s t-distribution
# • Negative binomial
# • Geometric
# • Exponential
# • Beta
# • Gamma

# We’ll leave it to you to work out how many parameters each has, and how to
# generate data using R.

# NOTE: see functions rt, rnbinom, rgeom, rexp, rbeta, rgamma for generating
# data from the distributions above

ggplot(data.frame(x = rt(1000, 999))) + geom_histogram(aes(x), binwidth = .1)
