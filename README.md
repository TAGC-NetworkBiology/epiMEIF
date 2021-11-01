# epiMEIF: Detecting high-order epistasis interactions using mixed effect conditional inference forest (epiMEIF)
# Article Information
Title : Epi-MEIF, a flexible and efficient novel method for detection of higher order epistatic interactions from GWAS 
Authors: Saswati Saha 1*, Laurent Perrin 1,2, Christine Brun 1,2, Lionel Spinelli 1*

(1) Aix Marseille Univ, INSERM, TAGC, UMR_S_1090, Turing Centre for Living systems, 13288 Marseille, France (2) CNRS, 13288 Marseille, France 
(*) Authors for correspondence (saswati.saha@uni-amu.fr, lionel.spinelli@univ-amu.fr)

Abstract: Understanding the relationship between genetic variations and variations in complex and quantitative phenotypes remains an ongoing challenge in any population. While Genome-wide association studies (GWAS) have become a vital tool for identifying single-locus associations, we lack methods for identifying higher order epistatic interactions. In this article, we propose a novel method for high-order epistasis detection using mixed effect conditional inference forest (epiMEIF). The epiMEIF model is fitted on a group of potential causal SNPs and the tree structure in the forest facilitates the identification of n-way interactions between the SNPs. Additional testing strategies improve the robustness of the method and thereafter provides a generalized way to identify higher order interactions from any GWAS data. We demonstrate the ability of the method to detect true n-way interactions via extensive simulations in both synthetic cross-sectional and longitudinal dataset. This is further illustrated in an application to reveal epistatic interactions from natural variations of heart period and cardiac aging in Drosophila. Overall, the epiMEIF is a high-performance flexible method for detecting higher-order epistasis interactions that can help us in identifying the genetic architecture linked to complex phenotypes and better understand the cellular and molecular pathways involved in complex GWAS studies. 
Keywords: GWAS, Random Forest, Epistatic interaction. 

# Goal of the github

This github project contains the instructions and material to implement the above method on any cross-sectional and longitudinal dataset. An RMD file (epiMEIF_Illustration) is loaded which retrieves data from the folder data and R codes from the folder codes and shows the step by step implementation of the epiMEIF method. 
In this html report (epiMEIF_Illustration) we illustrate how to run the cforest part of the epiMEIF and obtain higher order interactions from different types of dataset (cross-sectional dataset/longitudinal). We will also illustrate ways to validate the interactions obtained from MEIF using the additional testing strategies- max-t test and anova test and obtain the final interaction network.

# Description of the dataset
The dataset (in folder data) contains a similated cross-sectional and a longitudinal dataset that is simulated based on the cardiac phenotype data of drosopohila population (DGRP). The R code Creating_Age1_Dataset illustrates how the data are simulated.


# Description of the code

The R codes are saved in the folder codes. We have created two primary scripts: Interaction_Score_Age1 (for cross-sectional data), Interaction_Score_Ageing (for longitudinal data). The following section will elaborate on the different functions in Interaction_Score_Age1.

## Interaction_Score_Age1
### Functions
#### 1. cforest_gen

Generalized adaptation of cforest function from partykit R package that allows the user to apply a weighted conditional inference forest on any dataset.

#### Usage
```{r }
cforest_gen(formula, data, weights, subset, offset, cluster, strata, na.action = na.pass, control = ctree_control(teststat = "quad", testtype = "Univ", mincriterion = 0,
saveinfo = FALSE), ytrafo = NULL, scores = NULL, ntree = 500L, perturb = list(replace = FALSE, fraction = 0.632), mtry = ceiling(sqrt(nvar)), applyfun = NULL, cores = NULL, 
trace = FALSE, weight_variable = NULL)
```

#### Parameters
* formula- a symbolic description of the model to be fit.

* data- a data frame containing the variables in the model.

* subset- an optional vector specifying a subset of observations to be used in the fitting process.

* weights- an optional vector of weights to be used in the fitting process. Non-negative integer valued weights are allowed as well as non-negative real weights. Observations are sampled (with or without replacement) according to probabilities weights / sum(weights). The fraction of observations to be sampled (without replacement) is computed based on the sum of the weights if all weights are integer-valued and based on the number of weights greater zero else. Alternatively, weights can be a double matrix defining case weights for all ncol(weights) trees in the forest directly. This requires more storage but gives the user more control.

* offset- an optional vector of offset values.

* cluster- an optional factor indicating independent clusters. Highly experimental, use at your own risk.

* strata- 	an optional factor for stratified sampling.

* na.action- a function which indicates what should happen when the data contain missing value.

* control- a list with control parameters, see ctree_control. The default values correspond to those of the default values used by cforest from the party package. saveinfo = FALSE leads to less memory hungry representations of trees. Note that arguments mtry, cores and applyfun in ctree_control are ignored for cforest, because they are already set.

* ytrafo- an optional named list of functions to be applied to the response variable(s) before testing their association with the explanatory variables. Note that this transformation is only performed once for the root node and does not take weights into account (which means, the forest bootstrap or subsetting is ignored, which is almost certainly not a good idea). Alternatively, ytrafo can be a function of data and weights. In this case, the transformation is computed for every node and the corresponding weights. This feature is experimental and the user interface likely to change.

* scores-	an optional named list of scores to be attached to ordered factors.

* ntree- Number of trees to grow for the forest.

* perturb- a list with arguments replace and fraction determining which type of resampling with replace = TRUE referring to the n-out-of-n bootstrap and replace = FALSE to sample splitting. fraction is the portion of observations to draw without replacement. Honesty (experimental): If fraction has two elements, the first fraction defines the portion of observations to be used for tree induction, the second fraction defines the portion of observations used for parameter estimation. The sum of both fractions can be smaller than one but most not exceed one. Details can be found in Section 2.4 of Wager and Athey (2018).

* mtry- number of input variables randomly sampled as candidates at each node for random forest like algorithms. Bagging, as special case of a random forest without random input variable sampling, can be performed by setting mtry either equal to Inf or manually equal to the number of input variables.

* applyfun- an optional lapply-style function with arguments function(X, FUN, ...). It is used for computing the variable selection criterion. The default is to use the basic lapply function unless the cores argument is specified (see below).

* cores- numeric. If set to an integer the applyfun is set to mclapply with the desired number of cores.

* trace-	a logical indicating if a progress bar shall be printed while the forest grows.

* object- An object as returned by cforest

* newdata- An optional data frame containing test data.

* type- a character string denoting the type of predicted value returned, ignored when argument FUN is given. For "response", the mean of a numeric response, the predicted class for a categorical response or the median survival time for a censored response is returned. For "prob" the matrix of conditional class probabilities (simplify = TRUE) or a list with the conditional class probabilities for each observation (simplify = FALSE) is returned for a categorical response. For numeric and censored responses, a list with the empirical cumulative distribution functions and empirical survivor functions (Kaplan-Meier estimate) is returned when type = "prob". "weights" returns an integer vector of prediction weights. For type = "where", a list of terminal node ids for each of the trees in the forest ist returned.

* OOB- a logical defining out-of-bag predictions (only if newdata = NULL). If the forest was fitted with honesty, this option is ignored.

* FUN- a function to compute summary statistics. Predictions for each node have to be computed based on arguments (y, w) where y is the response and w are case weights.

* simplify- a logical indicating whether the resulting list of predictions should be converted to a suitable vector or matrix (if possible).

* scale- a logical indicating scaling of the nearest neighbor weights by the sum of weights in the corresponding terminal node of each tree. In the simple regression forest, predicting the conditional mean by nearest neighbor weights will be equivalent to (but slower!) the aggregation of means.

* tree- an integer, the number of the tree to extract from the forest.

* \boldsymbol{weight_variable = an optional vector of weights for the predictors to be sampled at the root node. In the cforest function of the partykit package there was no option to assign weight to the predictors sampling scheme. Each predictor was equally probable to be sampled at the root node. We added this feature in our adaptation where the user can assigned weights to all the variables based on some pre-decided criteria.}
