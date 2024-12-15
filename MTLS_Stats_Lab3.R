# Lab 3: Inference and confidence intervals

# For the first half of this session we would like you to calculate standard
# errors and confidence "by hand". You are welcome to use a calculator.
#
# In the second half of the session you will use R to estimate standard errors
# and construct confidence intervals.
#
# Please note that all data are fictitious and have been created for educational
# purposes

# PART B: Using R
#
# Open the BeatAML2 dataset. R commands to complete these questions are given in
# “MTLS_Stats_Lab3_inference.R”, available on Canvas. We will be working with
# the clinical_data data.

# load useful R packages
library(ggplot2)
library(dplyr)

# load course data
# NOTE: path will be different for students
load("E:/Xuexi/5MT013/5MT013/BeatAML2_clean.RData")

#
# 4.
# a) Plot a histogram to show the distribution of creatinine levels.
ggplot(clinical_data) + geom_histogram(aes(x = creatinine), bins = 30)

# b) What is the mean level among people who provided a valid measurement?
mean(clinical_data$creatinine, na.rm = TRUE)
# c) What is the standard deviation of these measurements?
summary_cre <- clinical_data |>
    filter(!is.na(creatinine)) |>
    summarize(mean_cre = mean(creatinine),
              sd_cre = sd(creatinine),
              n = n())
print(summary_cre)
# d) Calculate (by hand) the standard error of the estimate of the mean
#    creatinine level.
sd(clinical_data$creatinine, na.rm = TRUE)/sqrt(262)
# e) Use R to confirm your calculation of the standard error.
summary_cre <- summary_cre |> mutate(se = sd_cre/sqrt(n))
print(summary_cre)

# f) Use R to construct a 95% confidence interval
# NOTE: here a couple of possible ways to do this

# 1. You manually do: mean ± z * SE. Here z = 1.96 for 95% CI
summary_cre$mean_cre + c(-1,1)*qnorm(0.975)*summary_cre$se

#This is functionally identical to the first method, but uses the with() function
#so you don’t have to repeatedly refer to summary_cre$....
#Inside with(), mean_cre and se are directly accessible.
#The logic is the same: a direct, manual calculation of the CI using a normal approximation.
with(summary_cre, mean_cre + c(-1,1)*qnorm(0.975)*se)

# 2. Use the confint() function on the output of lm() or t.test()
# fitting a linear model (lm) with creatinine ~ 1 (which means an intercept-only model, effectively estimating the mean of creatinine).
# By default, confint() on lm objects uses the t-distribution (assuming normality of residuals) and the appropriate degrees of freedom.
confint(lm(creatinine ~ 1, data = clinical_data))

t.test(creatinine ~ 1, data = clinical_data)$conf.int

#
# 5. Now we will look at the variable “specimentype”
# a) What proportion of people provided a bone marrow aspirate sample?
clinical_data |>
    count(specimenType)
count(clinical_data, specimenType)
count(clinical_data, specimenType) |> mutate(prop = n/sum(n))
# b) Among those for whom we have a valid creatinine measurement, what
#    proportion of people provided a bone marrow aspirate sample?
prop_bm <- clinical_data |>
    filter(!is.na(creatinine)) |>
    with(table(marrow = (specimenType == "Bone Marrow Aspirate")))

prop_bm/sum(prop_bm)

# c) Use R to calculate three confidence intervals for this estimate: 90%, 95%
#    and 99%
# NOTE: a couple of different approaches here too
cl <- c(ci_90 = 0.90, ci_95 = 0.95, ci_99 = 0.99)
# 1. Using prop.test() with the correct number of successes and failures
sapply(cl, \(x) prop.test(182, 182 + 80, conf.level = x)$conf.int)
# 2. Using prop.test() with the raw data
sapply(cl, \(x) prop.test(prop_bm[2:1], conf.level = x)$conf.int)

# d) Compare the widths of these intervals. Are they what you expected? Is the
#    estimate in the center of each interval?





