# Biostatistics---IRT-functions

Here you may find the description of every R-functions written in the "internship_work" file. The purpose of this work is to facilitate the study of databases of order polytomous items responses and to apply them IRT models. 

# The packages required :

We need several packages to correctly run the functions we are about to describe :
- the "eRm" package : fits Rasch models for dichotomous items ;
- the "sirt" package : for a specific purpose, which is to test local independence ;
- the "ltm" package : fits Rasch models for polytomous items (it is the most useful package we will use) ;
- the "mirt" package : fits IRT models for multidimensionnal latent variable ;

# Likelihood Ratio Test :

In this study, we need to test DIF (Differential Item Functioning). What is DIF ? DIF occurs for an item when people with the same ability but from different groups (groups defined by a covariate, such as gender or age) have unequal probabilities of giving a response. In other words, the conditional probability of success for a DIF-item depends on the value of a covariate. To apply Rasch models, we must have no DIF items. We then need to test that hypothesis when we apply a Rasch model. To check if an item is DIF for a specific covariate, we proceed as follows :

- from the DIF-suspected item and the covariate, we create new items : one for each level of the covariate, equal to the original item when the covariate values the studied level, and equal to NA otherwise. We take an example : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women"); otherwise "item1_2" is equal to NA ;
- we may then apply two Rasch models : one model to the database with the original item, and one model to the database with the new items we just created (and at the same time we exclude the original item) ;
- then, we may apply a LR-test : if the p-value is under a specified thereshold, then we can reject the null hypothesis H0 : "there is no DIF". 

To this end, we use the "LRtest" function, which very basically apply a log-likelihood ratio test. The function returns the value of the LRTest, the degrees of freedom and the p-value. 




