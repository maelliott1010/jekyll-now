---
layout: post
title: A Factor Analytic Approach to Understanding Adult ADHD: A Walkthrough of the CFA performed in Park et al. (2018)
---

Files associated with this post can be found [here](https://github.com/maelliott1010/CFA-ADHD-Park_et_al.).

**Before we get started, some background on ADHD symptoms...**
The expected factor structure of ADHD symptoms has undergone many changes recently. Historically, the creation and validation of this factor structure was based on primarily juvenile male populations. Although it is typically assumed that the structure of symptoms holds across the lifespan, empirical study of the factor structure of ADHD symptoms in adults is quite recent. 

Our analysis systematically examines the one-factor, two-factor, three-factor, and bifactor models of adult ADHD symptoms. 
- First, we evaluate the fit statistics of each of the models. 
- Next, we obtain estimates of model-based reliability (e.g., coefficient omega and coefficient omega hierarchical; Rodriguez, Reise, & Haviland, 2015).
- We then estimate the stability and construct replicability for each factor (e.g., H-factor; Hancock & Mueller, 2001).
- For the bifactor models, we examine estimates for the unidimensionality of the scale (i.e., explained common variance [ECV]; Rodriguez et al., 2015). 
- As exploratory analyses, we examine the invariance of the models across both gender and datasets. 
- Finally, we explore model validity by regressing variables with known associations with ADHD symptoms, such as education, depression, hostility, and both positive and negative parenting on the latent factors as specified in the best-fitting models. 

**In this post, I will walk through the R code used for our 2018 paper, which takes a confirmatory factor analysis (CFA) approach to addressing the factor structure of ADHD symptoms in adults.**

**Data Cleaning** <br/>
<br/>
Before I get into the analysis, here is some information about our dataset:
- We started with a dataset from three samples from parents of children with and without ADHD (n = 673, 430 females, and 243 males). 
- Getting the right subset of data for the analyses was a bit tricky. We had to worry about a couple important factors: *independence* and *equal n for males and females* For instance, in some cases, only mothers participated; in others, the men and women were couples parenting the same child.
- To tackle these issues, we performed the following data cleaning: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- To ensure *independence* of the data, male and female partners from within the same family were not included (i.e., only one parent from each two-parent family was randomly chosen to be included). <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- In order to create a final dataset with an *equal n for males and females*, we randomly selected 215 males out of the 243 males, and then included the remaining females who were not coupled with the 215 males. 
- The resulting sample included 215 males and 215 females.

**Scale Conversion** <br/>
<br/>
Self-report of ADHD symptoms was assessed using the Current Symptom Scale-Self-Report (CSS; Barkley & Murphy, 2006) in Study 1 and 2, and the Barkley Adult ADHD Rating Scale-IV (BAARS; Barkley, 2011) in Study 3. <br/>
Some quick notes about using two different scales:
- Both the BAARS and the CSS contain 18 items that assess the presence of DSM-5 symptoms for ADHD. The CSS is rated on a scale from 0 to 3, while BAARS is rated from 0 to 4.
- The BAARS is an updated version of the CSS, with the only differences being additional wording to make symptoms more applicable to adults. 
- Since the CSS is rated on a scale from 0 to 3, items on the BAARS had to be converted to a scale from 0 to 3 also. 

**The Analyses** <br/>
<br/>
Tools:
- All analyses were conducted in R with R Studio
- We used the *lavaan* 0.5-22  package for our main CFA analyses, more about *lavaan* [here](http://lavaan.ugent.be).
- We chose to run two different estimators in this analysis. We'll start this tutorial with robust maximum likelihood estimator (MLR), which produces standard errors and a chi-square test statistic that are robust to non-normality for incomplete data. 
- *Note*: Our data is on a 4-point scale, which is on the threshold of recommendations to either utilize a categorical or continuous estimation method (Rhemtulla, Brosseau-Liard, & Savalei, 2012). *Therefore, we also ran our analyses utilizing the unweighted least squares estimator (ULSMV) in order to account for the categorical nature of our data. Given that the pattern of results was highly similar using either the MLR or the ULSMV estimators, we only reported findings with the MLR estimator in our paper. Code for the ULSMV estimator is provided at the bottom half of this post.*  
<br/>

First, we load the **Hancock Mueller** function. *Its utility in this analysis is explained below*.
```r
#Hancock Mueller#
HancockMueller<-function(x){
  loads<-inspect(x,"std")$lambda
  facs<-colnames(loads)
  out<-rep(NA,length(facs))
  names(out)<-facs
  
  for(i in facs){
    out[i]<-(1+(1/sum((loads[,i]^2)/(1-loads[,i]^2))))^-1
  }
  out
}
```
We start with the simple One-Factor Model. We load the data in and specify our model in the code, then run the CFA and ask for a summary:
```r
##One-Factor Model 

Modsa1 <- '
ADHD =~ NA*CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18
ADHD~~1*ADHD
'
Runsa1 <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, missing="FIML", estimator="MLR")
summary(Runsa1,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE)
```
We use the comparative fit index (CFI), root mean square error of approximation (RMSEA), and standardized root mean square residual (SRMR) to evaluate model fit. We use the full information maximum likelihood (FIML), or direct maximum
likelihood, to handle missing data. FIML does not estimate missing values, and instead estimates relationships using all available data, assuming multivariate normality and MCAR or MAR. This function may need to be modified in replication attempts to meet the unique nature of missing observations in new samples. A common alternative approach is mulitple imputations (MI), which is more flexible than FIML, but comes with several costs. More on missing data techniques [here](https://statisticalhorizons.com/wp-content/uploads/Allison-2003-JAP-Special-Issue.pdf). <br/>
<br/>
Interpretation: 
- A value close to or greater than .95 for the TLI and the CFI, and a value close to or less than .06 for the RMSEA (Hu & Bentler, 1998) and SRMR less than .08 (Kline, 2016) indicate a good fit between the model and the observed data. 
Chi-square difference tests were also utilized to determine which model provided a significantly better fit to the data.

```r
#reliability
reliability(Runsa1) 
```
Here, we compute the reliabilities for each of the factors. This is equivalent to the proportion of variance in the indicators of each factor that were accounted for by that factor specifically (ωh), as well as the proportion of variance in all items that were accounted for by all factors together (ω; e.g., Reise, Bonifay, & Haviland, 2013). <br/>
<br/>
Interpretation: 
- Recommendations by Reise et al. suggest a minimual acceptable value of ωh of .50, though a value of .75 is more acceptable. 
- Comparison of ω and ωh values can indicate how much reliable variance could be attributed to general vs. specific factors.

```r
#Hancock
HancockMueller(Runsa1)
```
Here, calculate the H-index (Hancock & Mueller, 2001) to examine construct replicability of each factor, indicating how well a set of items represents a latent variable.  <br/>
<br/>
Interpretation:
- High H values (>.70) indicate the factor is stable and is less likely to fluctuate between various samples. 

```r
#Measurement invariance: dataset
Runsa1A <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Dataset", missing="FIML", estimator="MLR", std.lv=TRUE)
Runsa1B <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Dataset", group.equal=c("loadings"), missing="FIML", estimator="MLR",std.lv=TRUE)
Runsa1C <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Dataset", group.equal=c("loadings", "intercepts"), missing="FIML", estimator="MLR", std.lv=TRUE)
lavTestLRT(Runsa1A,Runsa1B)
lavTestLRT(Runsa1A,Runsa1C)
```
Here, we assess Measurement invariance to examine whether models were equivalent across gender and across datasets, which would provide evidence that these models hold in different subpopulations. To do this, we tested weak invariance first by comparing the fit of a configural model (i.e., a model that was fit to data for males and females, but no constraints were imposed) with that of the same model but with all factor loadings constrained to be equal across gender.  <br/>
<br/>
Interpretation:
- A significant chi-squared difference test of invariance, in the case of our data, would suggest that a hypothesis that the model is invariant across a given variable (i.e., gender, dataset) can be rejected. *Note*: this was only done for nested models that were tested using the NET procedure

```r
#Measurement invariance: Gender
Runsa1AA <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Gender",estimator="MLR", std.lv=TRUE, missing="FIML")
Runsa1BB <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Gender", group.equal=c("loadings"), estimator="MLR", std.lv=TRUE, missing="FIML")
Runsa1CC <- cfa(Modsa1, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Gender", group.equal=c("loadings", "intercepts"), estimator="MLR", std.lv=TRUE, missing="FIML")
lavTestLRT(Runsa1AA,Runsa1BB)
lavTestLRT(Runsa1AA,Runsa1CC)
```
If weak invariance held, we then tested strong invariance by comparing the configural model with the same model but with all factor loadings as well as intercepts constrained to be equal across gender. The same analysis was performed across datasets. However, it is necessary to note that these analyses are exploratory given the limited sample size.

Finally, we  examined the concurrent validity of the models demonstrating the best fit by regressing participants’ education, depression, hostility, positive parenting, and lax parenting on each of the latent factors as specified in their respective configurations.

**Repeat these procedures for Two- and Three- Factor models:**
```r
##Two-Factor Model

Modsa2<- 
  'INATT=~ NA*CuSS_S1 + CuSS_S3 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17

HYP=~ NA*CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 + CuSS_S14 + CuSS_S16 + CuSS_S18

INATT~~1*INATT
HYP~~1*HYP

INATT~~HYP
'

Runsa2 <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, missing="FIML", estimator="MLR")

summary(Runsa2,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

#reliability
reliability(Runsa2)
HancockMueller(Runsa2)

#Measurement invariance: dataset
Runsa2A <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset")
Runsa2B <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset", group.equal=c("loadings"))
Runsa2C <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset", group.equal=c("loadings", "intercepts"))
lavTestLRT(Runsa2A,Runsa2B)
lavTestLRT(Runsa2A,Runsa2C)

#Measurement invariance: Gender
Runsa2AA <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender")
Runsa2BB <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender", group.equal=c("loadings"))
Runsa2CC <- cfa(Modsa2, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender", group.equal=c("loadings", "intercepts"))
lavTestLRT(Runsa2AA,Runsa2BB)
lavTestLRT(Runsa2AA,Runsa2CC)

##Three-Factor Model

Modsa3<- 
  'INATT=~ NA*CuSS_S1 + CuSS_S3 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17

IMP=~ NA*CuSS_S14 + CuSS_S16 + CuSS_S18

HYP=~ NA*CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12
INATT~~1*INATT
IMP~~1*IMP
HYP~~1*HYP
INATT~~HYP
IMP~~INATT
HYP~~IMP
'
Runsa3 <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, missing="FIML", estimator="MLR")
summary(Runsa3,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

#reliability
reliability(Runsa3)
HancockMueller(Runsa3)

#Measurement invariance: dataset
Runsa3A <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset")
Runsa3B <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset", group.equal=c("loadings"))
Runsa3C <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset", group.equal=c("loadings", "intercepts"))
lavTestLRT(Runsa3A,Runsa3B)
lavTestLRT(Runsa3A,Runsa3C)

#Measurement invariance: Gender
Runsa3AA <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender")
Runsa3BB <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender", group.equal=c("loadings"))
Runsa3CC <- cfa(Modsa3, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender", group.equal=c("loadings", "intercepts"))
lavTestLRT(Runsa3AA,Runsa3BB)
lavTestLRT(Runsa3AA,Runsa3CC)
```
**Bifactor Models**

```r
##Bifactor model (INATT/HYP)

Modsa4 <- '
ADHD =~ CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18
INATT=~CuSS_S3 + CuSS_S1 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17
HYP=~CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 + CuSS_S14 + CuSS_S16 + CuSS_S18

ADHD~~0*INATT
ADHD~~0*HYP
INATT~~0*HYP
'
Runsa4 <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, estimator="MLR", missing="FIML",std.lv=TRUE)
summary(Runsa4,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

#Reliability
reliability(Runsa4) 
HancockMueller(Runsa4)

#Measurement invariance: dataset
Runsa4A <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR",  std.lv=TRUE, missing="FIML", group="Dataset")
Runsa4B <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset", group.equal=c("loadings"))
Runsa4C <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Dataset", group.equal=c("loadings", "intercepts"), missing="FIML", estimator="MLR", std.lv=TRUE)
lavTestLRT(Runsa4A,Runsa4B)
lavTestLRT(Runsa4A,Runsa4C)

#Measurement invariance: Gender
Runsa4AA <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR",  std.lv=TRUE, missing="FIML", group="Gender")
Runsa4BB <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender", group.equal=c("loadings"))
Runsa4CC <- cfa(Modsa4, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Gender", group.equal=c("loadings", "intercepts"))
lavTestLRT(Runsa4AA,Runsa4BB)
lavTestLRT(Runsa4AA,Runsa4CC)
```
*Note:* For the bifactor models, we also computed the explained common variance (ECV) to examine the unidimensionality of the substantive ADHD factor (Rodriguez et al., 2015).
```r
#Explained common variance
L<-inspect(Runsa4, "coef")$lambda
lsq<-L*L
ECV.Runsa4<-sum(lsq[,1])/sum(lsq)
```
Here, we compute this using a function to estimate factor loadings of the general and secondary dimensions of a bifactor model. <br/>
<br/>
Interpretation:
- An ECV of greater than .80 indicates that the factor loadings of the general ADHD factor are very similar to those that might be obtained by the estimate of a one-dimensional model. 
- An ECV of greater than .70 indicates that the factor loadings of the general ADHD factor are increasingly similar to those that might be obtained by the estimate of a one-dimensional model (Rodriguez et al., 2015). 

```r
##Bifactor Model-3#

Modsa5 <- '
ADHD =~ CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18
INATT=~CuSS_S3 + CuSS_S1 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17
HYP=~CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 
IMP=~CuSS_S14+CuSS_S16+CuSS_S18


ADHD~~0*INATT
ADHD~~0*HYP
ADHD~~0*IMP
INATT~~0*HYP
INATT~~0*IMP
HYP~~0*IMP

'
Runsa5 <- cfa(Modsa5, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", missing="FIML",std.lv=TRUE)
summary(Runsa5,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

# Reliability

reliability(Runsa5) 
HancockMueller(Runsa5)

#Explained common variance
M<-inspect(Runsa5, "coef")$lambda
lsq<-M*M
ECV.Runsa5<-sum(lsq[,1])/sum(lsq)


#Measurement invariance: dataset
Runsa5A <- cfa(Modsa5, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR",  std.lv=TRUE, missing="FIML", group="Dataset")
Runsa5B <- cfa(Modsa5, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR,estimator="MLR", std.lv=TRUE, missing="FIML", group="Dataset", group.equal=c("loadings"))
Runsa5C <- cfa(Modsa5, dat=Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR, group="Dataset", group.equal=c("loadings", "intercepts"), missing="FIML", estimator="MLR", std.lv=TRUE)
lavTestLRT(Runsa5A,Runsa5B)
lavTestLRT(Runsa5A,Runsa5C)

#Measurement invariance: Gender
Runs5AA <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender")
Runs5BB <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings"))
Runs5CC <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings", "intercepts"))
lavTestLRT(Runs5AA,Runs5BB)
lavTestLRT(Runs5AA,Runs5CC)
```

**NET Procedure**

Below, we provide code for the NET procedure, which allows us to determine whether or not the 3-factor model is nested within Bifactor model 2 and Bifactor model 3.

```r

#Model 1 - This is the 3-factor model, to test whether this model is nested within Bifactor 2 and Bifactor 3#

Model1<- 

  'INATT=~ NA*CuSS_S1 + CuSS_S3 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17



IMP=~ NA*CuSS_S14 + CuSS_S16 + CuSS_S18



HYP=~ NA*CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12

INATT~~1*INATT

IMP~~1*IMP

HYP~~1*HYP

INATT~~HYP

IMP~~INATT

HYP~~IMP

'

Run1 <- cfa(Model1, dat=dataset, estimator="MLR", missing="FIML",std.lv=TRUE)

summary(Run1,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 





#Model 2 - This is the bifactor 2 model#



Model2 <- '

ADHD =~ CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18

INATT=~CuSS_S3 + CuSS_S1 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17

HYP=~CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 + CuSS_S14 + CuSS_S16 + CuSS_S18



ADHD~~0*INATT

ADHD~~0*HYP

INATT~~0*HYP

'

Run2 <- cfa(Model2, dat=dataset, estimator="MLR", missing="FIML",std.lv=TRUE)

summary(Run2,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 









# Model 3 - This is the bifactor 3 model#



Model3 <- '

ADHD =~ CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18

INATT=~CuSS_S3 + CuSS_S1 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17

HYP=~CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 

IMP=~CuSS_S14+CuSS_S16+CuSS_S18


ADHD~~0*INATT

ADHD~~0*HYP

ADHD~~0*IMP

INATT~~0*HYP

INATT~~0*IMP

HYP~~0*IMP


'

Run3 <- cfa(Model3, dat=dataset,estimator="MLR", missing="FIML",std.lv=TRUE)

summary(Run3,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 



#----NET TEST 1: Is Model 2 nested within Model 1?-----------------


#Model-implied covariance matrix from Model 1 run: 

Sigma.Model1<-fitted.values(Run1)$cov   



#Model 2 fitted to the model-implied covariance matrix from Model 1 run:

Run2.nettest21 <- cfa(Model2, sample.cov=Sigma.Model1, sample.nobs=430, std.lv=TRUE) #N doesn't matter

summary(Run2.nettest21) 



#You see that the fit is perfect--these models are nested.



# > summary(Run2.nettest21)

# lavaan (0.5-23.1097) converged normally after  45 iterations

# 

# Number of observations                           430

# 

# Estimator                                         ML

# Minimum Function Test Statistic                0.000

# Degrees of freedom                               117

# P-value (Chi-square)                           1.000



#----NET TEST 2: Is Model 3 nested within Model 1?-----------------



#Model-implied covariance matrix from Model 1 run: 

Sigma.Model1<-fitted.values(Run1)$cov   



#Model 3 fitted to the model-implied covariance matrix from Model 1 run:

Run3.nettest31 <- cfa(Model3, sample.cov=Sigma.Model1, sample.nobs=430, std.lv=TRUE) #N doesn't matter

summary(Run3.nettest31) 



#You see that the fit is perfect--these models are also nested.

# > summary(Run3.nettest31)

# lavaan (0.5-23.1097) converged normally after  53 iterations

# 

# Number of observations                           430

# 

# Estimator                                         ML

# Minimum Function Test Statistic                0.000

# Degrees of freedom                               117

# P-value (Chi-square)                           1.000

```

**ULSMV Estimator**
As mentioned in the beginning of this post, our data is on a 4-point scale, which is on the threshold of recommendations to either utilize a categorical or continuous estimation method (Rhemtulla, Brosseau-Liard, & Savalei, 2012). The unweighted least squares estimator (ULSMV) can be used to account for the categorical nature of a dataset. *Given that the pattern of results was highly similar using either the MLR or the ULSMV estimators, we only reported findings with the MLR estimator in our paper. Nevertheless, we provide code for the ULSMV estimator here:*

```r
#Hancock Mueller#
HancockMueller<-function(x){
  loads<-inspect(x,"std")$lambda
  facs<-colnames(loads)
  out<-rep(NA,length(facs))
  names(out)<-facs
  
  for(i in facs){
    out[i]<-(1+(1/sum((loads[,i]^2)/(1-loads[,i]^2))))^-1
  }
  out
}

#Reducing dataset to 3 indicators#
tempdata2<-Bifactor_Dataset_Nov_8_17_recodedtocorrectCIHR
sapply(tempdata2[,6:24],table)

tempdata2[which(tempdata2$CuSS_S1==3),]$CuSS_S1<-2
tempdata2[which(tempdata2$CuSS_S2==3),]$CuSS_S2<-2
tempdata2[which(tempdata2$CuSS_S3==3),]$CuSS_S3<-2
tempdata2[which(tempdata2$CuSS_S4==3),]$CuSS_S4<-2
tempdata2[which(tempdata2$CuSS_S5==3),]$CuSS_S5<-2
tempdata2[which(tempdata2$CuSS_S6==3),]$CuSS_S6<-2
tempdata2[which(tempdata2$CuSS_S7==3),]$CuSS_S7<-2
tempdata2[which(tempdata2$CuSS_S8==3),]$CuSS_S8<-2
tempdata2[which(tempdata2$CuSS_S9==3),]$CuSS_S9<-2
tempdata2[which(tempdata2$CuSS_S10==3),]$CuSS_S10<-2
tempdata2[which(tempdata2$CuSS_S11==3),]$CuSS_S11<-2
tempdata2[which(tempdata2$CuSS_S12==3),]$CuSS_S12<-2
tempdata2[which(tempdata2$CuSS_S13==3),]$CuSS_S13<-2
tempdata2[which(tempdata2$CuSS_S14==3),]$CuSS_S14<-2
tempdata2[which(tempdata2$CuSS_S15==3),]$CuSS_S15<-2
tempdata2[which(tempdata2$CuSS_S16==3),]$CuSS_S16<-2
tempdata2[which(tempdata2$CuSS_S17==3),]$CuSS_S17<-2
tempdata2[which(tempdata2$CuSS_S18==3),]$CuSS_S18<-2

##One-Factor Model 

Mods1 <- '
ADHD =~ NA*CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18
ADHD~~1*ADHD
'
Runs1 <- cfa(Mods1, dat=tempdata2, missing="pairwise", estimator="WLSMV", zero.add="default", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"))
summary(Runs1,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE)

#reliability
reliability(Runs1) 

#Hancock
HancockMueller(Runs1)

#Measurement invariance: dataset
Runs1A <- cfa(Mods1, dat=tempdata2, group="Dataset", estimator="ULSMV", zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise")
Runs1B <- cfa(Mods1, dat=tempdata2, group="Dataset", group.equal=c("loadings"), missing="pairwise", estimator="ULSMV", zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE)
Runs1C <- cfa(Mods1, dat=tempdata2, group="Dataset", group.equal=c("loadings", "intercepts"), missing="pairwise", estimator="ULSMV", zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE)
lavTestLRT(Runs1A,Runs1B, method="satorra.bentler.2010")
lavTestLRT(Runs1A,Runs1C, method="satorra.bentler.2010")

#Measurement invariance: Gender
Runs1AA <- cfa(Mods1, dat=tempdata2, group="Gender",estimator="ULSMV", zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise")
Runs1BB <- cfa(Mods1, dat=tempdata2, group="Gender", group.equal=c("loadings"), estimator="ULSMV",zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise")
Runs1CC <- cfa(Mods1, dat=tempdata2, group="Gender", group.equal=c("loadings", "intercepts"),zero.add=c(0.5,0.5), estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise")
lavTestLRT(Runs1AA,Runs1BB, method="satorra.bentler.2010")
lavTestLRT(Runs1AA,Runs1CC, method="satorra.bentler.2010")

##Two-Factor Model

Mods2<- 
  'INATT=~ NA*CuSS_S1 + CuSS_S3 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17

HYP=~ NA*CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 + CuSS_S14 + CuSS_S16 + CuSS_S18

INATT~~1*INATT
HYP~~1*HYP

INATT~~HYP
'

Runs2 <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), zero.add=c(0.5,0.5), std.lv=TRUE, missing="pairwise")

summary(Runs2,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

#reliability
reliability(Runs2)
HancockMueller(Runs2)

#Measurement invariance: dataset
Runs2A <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), zero.add=c(0.5,0.5), std.lv=TRUE, missing="pairwise", group="Dataset")
Runs2B <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset", zero.add=c(0.5,0.5), group.equal=c("loadings"))
Runs2C <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset",zero.add=c(0.5,0.5), group.equal=c("loadings", "intercepts"))
lavTestLRT(Runs2A,Runs2B, method="satorra.bentler.2010")
lavTestLRT(Runs2A,Runs2C, method="satorra.bentler.2010")

#Measurement invariance: Gender
Runs2AA <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", zero.add=c(0.5,0.5))
Runs2BB <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings"), zero.add=c(0.5,0.5))
Runs2CC <- cfa(Mods2, dat=tempdata2, estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings", "intercepts"), zero.add=c(0.5,0.5))
lavTestLRT(Runs2AA,Runs2BB, method="satorra.bentler.2010")
lavTestLRT(Runs2AA,Runs2CC, method="satorra.bentler.2010")

##Three-Factor Model

Mods3<- 
  'INATT=~ NA*CuSS_S1 + CuSS_S3 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17

IMP=~ NA*CuSS_S14 + CuSS_S16 + CuSS_S18

HYP=~ NA*CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12
INATT~~1*INATT
IMP~~1*IMP
HYP~~1*HYP
INATT~~HYP
IMP~~INATT
HYP~~IMP
'
Runs3 <- cfa(Mods3, dat=tempdata2, estimator="ULSMV",ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18",std.lv=TRUE, missing="pairwise"),zero.add=c(0.5,0.5))
summary(Runs3,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

#reliability
reliability(Runs3)
HancockMueller(Runs3)

#Measurement invariance: dataset
Runs3A <- cfa(Mods3, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset",zero.add=c(0.5,0.5))
Runs3Ba <- cfa(Mods3, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset", group.equal=c("loadings"), zero.add=c(0.5,0.5))
Runs3C <- cfa(Mods3, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset", group.equal=c("loadings", "intercepts"),zero.add=c(0.5,0.5))
lavTestLRT(Runs3A,Runs3Ba, method="satorra.bentler.2010")
lavTestLRT(Runs3A,Runs3C, method="satorra.bentler.2010")

#Measurement invariance: Gender
Runs3AA <- cfa(Mods3, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender",zero.add=c(0.5,0.5))
Runs3BB <- cfa(Mods3, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings"),zero.add=c(0.5,0.5))
Runs3CC <- cfa(Mods3, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings", "intercepts"),zero.add=c(0.5,0.5))
lavTestLRT(Runs3AA,Runs3BB,method="satorra.bentler.2010")
lavTestLRT(Runs3AA,Runs3CC,method="satorra.bentler.2010")

##Bifactor model (INATT/HYP)

Mods4 <- '
ADHD =~ CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18
INATT=~CuSS_S3 + CuSS_S1 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17
HYP=~CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 + CuSS_S14 + CuSS_S16 + CuSS_S18

ADHD~~0*INATT
ADHD~~0*HYP
INATT~~0*HYP
'
Runs4 <- cfa(Mods4, dat=tempdata2, estimator="ULSMV",ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"),std.lv=TRUE,zero.add=c(0.5,0.5))
summary(Runs4,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

#Reliability
reliability(Runs4) 
HancockMueller(Runs4)

#Measurement invariance: dataset
Runs4A <- cfa(Mods4, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset",zero.add=c(0.5,0.5))
Runs4B <- cfa(Mods4, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset", group.equal=c("loadings"),zero.add=c(0.5,0.5))
Runs4C <- cfa(Mods4, dat=tempdata2, group="Dataset", group.equal=c("loadings", "intercepts"), missing="pairwise", estimator="ULSMV", zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE,zero.add=c(0.5,0.5))
lavTestLRT(Runs4A,Runs4B,method="satorra.bentler.2010")
lavTestLRT(Runs4A,Runs4C,method="satorra.bentler.2010")

#Measurement invariance: Gender
Runs4AA <- cfa(Mods4, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender",zero.add=c(0.5,0.5))
Runs4BB <- cfa(Mods4, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings"),zero.add=c(0.5,0.5))
Runs4CC <- cfa(Mods4, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings", "intercepts"),zero.add=c(0.5,0.5))
lavTestLRT(Runs4AA,Runs4BB,method="satorra.bentler.2010")
lavTestLRT(Runs4AA,Runs4CC,method="satorra.bentler.2010")

#Explained common variance
L<-inspect(Runs4, "coef")$lambda
lsq<-L*L
ECV.Run4<-sum(lsq[,1])/sum(lsq)

##Bifactor Model-3#

Mods5 <- '
ADHD =~ CuSS_S1+CuSS_S2+CuSS_S3+CuSS_S4+CuSS_S5+CuSS_S6+CuSS_S7+CuSS_S8+CuSS_S9+CuSS_S10+CuSS_S11+CuSS_S12+CuSS_S13+CuSS_S14+CuSS_S15+CuSS_S16+CuSS_S17+CuSS_S18
INATT=~CuSS_S3 + CuSS_S1 + CuSS_S5 + CuSS_S7 + CuSS_S9 + CuSS_S11 + CuSS_S13 + CuSS_S15 + CuSS_S17
HYP=~CuSS_S2 + CuSS_S4 + CuSS_S6 + CuSS_S8 + CuSS_S10 + CuSS_S12 
IMP=~CuSS_S14+CuSS_S16+CuSS_S18


ADHD~~0*INATT
ADHD~~0*HYP
ADHD~~0*IMP
INATT~~0*HYP
INATT~~0*IMP
HYP~~0*IMP

'
Runs5 <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", std.lv=TRUE, ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), missing ="pairwise", zero.add=c(0.5,0.5))
summary(Runs5,standardized=TRUE,fit.measures=TRUE, rsquare=TRUE) 

# Reliability

reliability(Runs5) 
HancockMueller(Runs5)

#Explained common variance
M<-inspect(Runs5, "coef")$lambda
lsq<-M*M
ECV.Run5<-sum(lsq[,1])/sum(lsq)


#Measurement invariance: dataset
Runs5A <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset",zero.add=c(0.5,0.5))
Runs5B <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Dataset", group.equal=c("loadings"),zero.add=c(0.5,0.5))
Runs5C <- cfa(Mods5, dat=tempdata2, group="Dataset", group.equal=c("loadings", "intercepts"), missing="pairwise", estimator="ULSMV", zero.add=c(0.5,0.5), ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE,zero.add=c(0.5,0.5))
lavTestLRT(Runs5A,Runs5B, method="satorra.bentler.2010")
lavTestLRT(Runs5A,Runs5C, method="satorra.bentler.2010")

#Measurement invariance: Gender
Runs5AA <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender",zero.add=c(0.5,0.5))
Runs5BB <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings"),zero.add=c(0.5,0.5))
Runs5CC <- cfa(Mods5, dat=tempdata2,estimator="ULSMV", ordered=c("CuSS_S1", "CuSS_S2","CuSS_S3", "CuSS_S4", "CuSS_S5", "CuSS_S6", "CuSS_S7", "CuSS_S8", "CuSS_S9", "CuSS_S10", "CuSS_S11", "CuSS_S12", "CuSS_S13", "CuSS_S14", "CuSS_S15", "CuSS_S16", "CuSS_S17", "CuSS_S18"), std.lv=TRUE, missing="pairwise", group="Gender", group.equal=c("loadings", "intercepts"),zero.add=c(0.5,0.5))
lavTestLRT(Runs5AA,Runs5BB, method="satorra.bentler.2010")
lavTestLRT(Runs5AA,Runs5CC, method="satorra.bentler.2010")
```
