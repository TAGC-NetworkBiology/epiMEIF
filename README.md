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

