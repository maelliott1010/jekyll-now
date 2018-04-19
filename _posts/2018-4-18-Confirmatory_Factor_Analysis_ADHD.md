---
layout: post
title: Factor Analytic Approach to Understanding Adult ADHD: A Walkthrough of the CFA performed in Park et al. (2018)
---

The expected factor structure of ADHD symptoms has undergone many changes recently. Historically, the creation and validation of this factor structure was based on primarily juvenile male populations. Although it is typically assumed that the structure of symptoms holds across the lifespan, empirical study of the factor structure of ADHD symptoms in adults is quite recent. 

This study analysis systematically examines the one-factor, two-factor, three-factor, and bifactor models of adult ADHD symptoms. 
- First, we evaluate the fit statistics of each of the models. Next, we obtain estimates of model-based reliability (e.g., coefficient omega and coefficient omega hierarchical; Rodriguez, Reise, & Haviland, 2015).
- For the bifactor models, we examine estimates for the unidimensionality of the scale (i.e., explained common variance [ECV]; Rodriguez et al., 2015), and estimates of stability and construct replicability for each factor (e.g., H-factor; Hancock & Mueller, 2001). 
- As exploratory analyses, we examine the invariance of the models across both gender and datasets. 
- Finally, we explore model validity by regressing variables with known associations with ADHD symptoms, such as education, depression, hostility, and both positive and negative parenting on the latent factors as specified in the best-fitting models. 

**In this post, I will walkthrough the R code used for our 2018 paper, which takes a confirmatory factor analysis (CFA) approach to addressing the factor structure of ADHD symptoms in adults.**



