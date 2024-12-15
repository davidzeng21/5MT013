#L2.1 This question concerns the example used in the lecture. Suppose that in a random sample of
# 178 girls without any copies of the mutated allele of FTO, the mean BMI is 21.7 and the standard
#deviation is 6.67. Suppose also that in a random sample of 100 girls with one or more copies of
#the mutated allele of FTO, the mean BMI is 24.0 and the standard deviation is 6.40. Go through
#the steps below for carrying out a t-test (which does not assume equal variances in the two
#subpopulations).
# 1. What is the null hypothesis? What is the alternative hypothesis?
# The null hypothesis is there is no difference in the mean BMI between the two groups.

#The alternative hypothesis is that there is a difference in the mean BMI between the two groups.

# 2. Estimate the standard error of the difference in means of BMI (under the null hypothesis):

6.67/(178^0.5) # standard error for wt group = 0.5
6.40/(100^0.5) # standard error for mut group = 0.64
(0.5^2+0.64^2)^0.5 # standard error of the difference in means = 0.81

# 3. Now calculate the test statistic value
(21.7-24.0)/0.81 # t-statistic value = -2.84 < -1.96 (critical value for 5% significance level)

# 4. For Welch’s t-test the degrees of freedom in this case is equal to 212.44 (this is a detail you
# don’t need to focus on). These are large samples so the p-value can be calculated either using
# a t-distribution with 212.44 degrees of freedom or a standard normal distribution. They give
#similar results - both yield a p-value of approximately 0.005. What do you conclude?

# The p-value is less than 0.05, so we reject the null hypothesis and conclude that there is a
# significant difference in the mean BMI between the two groups.

#5. For the same samples, in the lecture notes we compare the proportions of girls that have a
#BMI>30 between the two groups, by testing for a difference in proportions. Check that you
#understand the steps of this test by looking carefully through slides 18-20. Do you understand
#why the standard error for the difference in proportions for this test is different for that which
#you would need to compute if you were to construct a confidence interval for the difference in
#proportions?

# When testing a null hypothesis of no difference in proportions,
# the standard error is computed under the assumption that the proportions
# are equal under the null. In that scenario, a pooled estimate of the common proportion is used.
# This pooled estimate replaces the individual group proportions in the calculation
# of the standard error, thereby yielding a different formula than you would use
# if you were constructing a confidence interval (where you do not assume equality of
# the two proportions and do not pool).

# For proportions, the difference in how the standard error is computed for a hypothesis test
#vs. a confidence interval is due to the pooled estimate used under the null hypothesis
#(for the test) versus separate estimates for each group (for the confidence interval).


#L2.2 Return to the BeatAML2 data and the normalised expression of gene JOSD1. Now compare
#the two groups of patients with and without mutations in the FLT3 gene.
#First read in the data, then type the commands below to re-create the variables.

# Load the data
load("E:/Xuexi/5MT013/5MT013/BeatAML2_clean.RData")

# select the gene of interest and the FLT3-ITD mutation carriership variable
JOSD1<-rnaNorm[2222,]
FLT3<-clinical_data$'FLT3-ITD'

# 1. Visualise the distributions of normalised expression levels of gene JOSD1 in the two groups.
par(mfrow=c(1,2))
hist(JOSD1[which(FLT3=="negative")],main="FLT3 negative")
hist(JOSD1[which(FLT3=="positive")],main="FLT3 positive")

# 2. Carry outWelch’s t-test (which does not assume equal variances) at significance level α = 0.05:
t.test(JOSD1~FLT3)

# 3. Discuss the assumptions of the test and their validity. What do you conclude?

# The assumptions of the t-test are that the data are normally distributed
# The validity of these assumptions should be checked before interpreting the results.
# The t-test results indicate that there is not a significant difference in the mean expression of JOSD1
# between the two groups of patients with and without FLT3 mutations.

# Check the normality
shapiro.test(JOSD1[which(FLT3=="negative")])
shapiro.test(JOSD1[which(FLT3=="positive")])
# qq plot
library(car)
qqPlot(JOSD1[which(FLT3=="negative")])
qqPlot(JOSD1[which(FLT3=="positive")])

# 4. Check that you understand how the test statistic (above in point 2) is calculated by typing
# and understanding the commands below, i.e. how it is constructed from the data, and that
# you understand how the p-value is calculated (what distribution the test statistic is compared
# to – and how it is compared).

# The t-test statistic is calculated as the difference in means divided by
#the standard error of the difference in means.
# The p-value is calculated by comparing the t-test statistic to a t-distribution with
# degrees of freedom equal to the Welch-Satterthwaite approximation.

x1<-mean(JOSD1[which(FLT3=="negative")])
x2<-mean(JOSD1[which(FLT3=="positive")])
s1<-sd(JOSD1[which(FLT3=="negative")])
s2<-sd(JOSD1[which(FLT3=="positive")])
n1<-length(JOSD1[which(FLT3=="negative")])
n2<-length(JOSD1[which(FLT3=="positive")])
se<-(s1^2/n1+s2^2/n2)^0.5
T<-((x1-x2)-0)/se
print(T)
df<-(s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
print(df)
# How is df calculated?
# The degrees of freedom for Welch's t-test are calculated using the Welch-Satterthwaite approximation.
# This formula is based on the variances and sample sizes of the two groups.

# Interpret the results.
# The p-value is greater than 0.05, so we do not reject the null hypothesis and
#conclude that there is not a significant difference in the mean expression of
#JOSD1 between the two groups of patients with and without FLT3 mutations.

2*(1-pt(T,df))


# Check that you get the same number as obtained in point 2 and that you understand which
# area of the t-distribution corresponds to the p-value.

# The p-value is the area under the t-distribution curve to the left and the right of the t-test statistic.

#L2.3 Take a subset of (40) observations of the data in L2.1. Examine the distribution of the
#expression levels in one of the groups. Carry out Welch’s t-test and Wilcoxon’s rank sum test.

# Type the following code to take a sample of 40 patients and store their gene
# expression values of JOSD1 as JOSD1_sub and their mutation status of FLT3 as
# FLT3_sub
id<-sample(c(1:length(JOSD1)),40,replace=FALSE)
JOSD1_sub<-JOSD1[id]
FLT3_sub<-FLT3[id]
# Construct a Q-Q plot to examine the distribution of JOSD1 expression values in
# the 40 patients with FLT3 mutations
library(car)
qqPlot(JOSD1_sub[which(FLT3_sub=="positive")])
qqPlot(JOSD1_sub[which(FLT3_sub=="negative")])

# Carry out Welch’s t-test and Wilcoxon’s rank sum test
t.test(JOSD1_sub~FLT3_sub)
wilcox.test(JOSD1_sub~FLT3_sub)

# Is there any difference in the results of the two tests?
# The t-test and Wilcoxon rank sum test yield different results.
# The t-test indicates a significant difference in the mean expression of JOSD1 between the two groups,
# while the Wilcoxon rank sum test does not.

#L2.4 Use the BeatAML2 cohort data to study response to treatment of the drug Venetoclax. There
#are 242 persons in the study that are given this drug. We will use drug response values based on
#area under the curve (AUC).

# First, identify the individuals receiving this treatment:
partic<-drug_data$sample[which(drug_data$inhibitor=="Venetoclax")]
# And store their AUC values:
auc4partic<-drug_data$auc[which(drug_data$inhibitor=="Venetoclax")]
# Now link this information to mutation status:
mutstatus<-matrix(,length(partic),1)
for (a in 1:length(partic)) {
  mutstatus[a]<-FLT3[which(clinical_data$sample==partic[a])] }

#We will analyse AUC as a binary variable (below/above 156 (156 happens to be the median value))

auct<-cut(auc4partic,breaks=c(-Inf,156,+Inf),labels=c("0","1"))
table(mutstatus,auct)


# Now carry out a two-sample test of proportions. What are the null and alternative hypotheses?
# What is the test result?

prop.test(x = c(99, 22), n = c(178, 64))

# The null hypothesis is that there is no difference in the proportion of patients with AUC values
# above 156 between the two groups.
# The alternative hypothesis is that there is a difference in the proportion of patients with AUC values
# above 156 between the two groups.
# The test results indicate that there is a significant difference in the proportion of patients with AUC
# values above 156 between the two groups.












