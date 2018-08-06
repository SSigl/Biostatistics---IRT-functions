# Biostatistics---IRT-functions

Here you may find the description of every R-functions written in the "internship_work" file. The purpose of this work is to facilitate the study of databases of order polytomous items responses and to apply them IRT models. 

# The packages required :

We need several packages to correctly run the functions we are about to describe :
- the "eRm" package : fits Rasch models for dichotomous items ;
- the "sirt" package : for a specific purpose, which is to test local independence ;
- the "ltm" package : fits Rasch models for polytomous items (it is the most useful package we will use) ;
- the "mirt" package : fits IRT models for multidimensionnal latent variable ;


# Likelihood Ratio Test : "LRtest" (intermediary function)

In this study, we need to test DIF (Differential Item Functioning). What is DIF ? DIF occurs for an item when people with the same ability but from different groups (groups defined by a covariate, such as gender or age) have unequal probabilities of giving a response. In other words, the conditional probability of success for a DIF-item depends on the value of a covariate. To apply Rasch models, we must have no DIF items. We then need to test that hypothesis when we apply a Rasch model. To check if an item is DIF for a specific covariate, we proceed as follows :

- from the DIF-suspected item and the covariate, we create new items : one for each level of the covariate, equal to the original item when the covariate values the studied level, and equal to NA otherwise. 
We take an example : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA ;

- we may then apply two Rasch models : one model to the database with the original item, and one model to the database with the new items we just created (and at the same time we exclude the original item) ;

- then, we may apply a LR-test : if the p-value is under a specified thereshold, then we can reject the null hypothesis : "there is no DIF". 

To this end, we wrote the "LRtest" function, which very basically apply a log-likelihood ratio test. This function takes in argument two models ("model0" and "model1"), which might be for example gpcm objects : the first model must be the one with the original item, and the second model the one with the items created from the covariate. This function has also a "display" argument : "display" must be a boolean, equal to "TRUE" if you wish to display the result in the R-Console. 
The function returns the value of the LR-test, the degrees of freedom and the p-value. 

This function will be used in several other functions and generally won't be used "on its own", that's why we qualify it as an "intermediary function". 



# To plot the expected value : "esperance" (intermediary function)

In the case of polytomous items, we cannot plot the probability of success as a function of the latent variable eta, as we can in the case of dichotomous items : we do not interpret the dichotomous and polytomous items the same way. However, we may plot the expected value of an item for a gpcm-model applied to a database with polytomous items. 

To this end, we wrote the "esperance" function. This functions takes in argument one database ("data"), one item ("item"), one itemlist ("itemlist") and one constraint ("constraint"). The database is the database we extract the information we need from. We restrain this database to an itemlist to apply a gpcm model ("gpcm" stands for "generalized partial credit model"), which might be a 1PL or 2PL, according to the constraint we choose : 
> if the constraint is "rasch", the discrimination parameter is assumed equal for all items and fixed to one ;
> if the constraint is "1PL", the discrimination parameter is assumed equal for all items but is estimated ;
> if the constraint is "gpcm", the discrimination parameters differ between the items and are estimated.
For more details on the gpcm objects, we recommand the reading of the description of the "ltm" R-package.
. The item is the item whom we wish to calculate the expected value. 






















