# Biostatistics---IRT-functions

Here you may find the description of every R-functions written in the "internship_work" file. The purpose of this work is to facilitate the study of databases of order polytomous items responses and to apply them IRT models. 

# The packages required :

We need several packages to correctly run the functions we are about to describe :
- the "eRm" package : fits Rasch models for dichotomous items ;
- the "sirt" package : for a specific purpose, which is to test local independence ;
- the "ltm" package : fits Rasch models for polytomous items (it is the most useful package we will use) ;
- the "mirt" package : fits IRT models for multidimensionnal latent variable ;

# Likelihood Ratio Test :

In this study, we need to test DIF (Differential Item Functioning). What is DIF ? DIF occurs for an item when people with the same ability but from different groups (groups defined by a covariate, such as gender or age) have unequal probabilities of giving a response. In other words, the conditional probability of success for a DIF-item depends on the value of a covariate. To apply Rasch models, we must have no DIF items. We then need to test that hypothesis when we apply a Rasch model. To this end, we proceed the following way :

 



