# Biostatistics---IRT-functions

Here you may find the description of every R-functions written in the "internship_work" file. The purpose of this work is to facilitate the study of databases of order polytomous items responses and to apply them IRT models. 

# The packages required :

We need several packages to correctly run the functions we are about to describe :
- the "eRm" package : fits Rasch models for dichotomous items ;
- the "sirt" package : for a specific purpose, which is to test local independence ;
- the "ltm" package : fits Rasch models for polytomous items (it is the most useful package we will use) ;
- the "mirt" package : fits IRT models for multidimensionnal latent variable ;


# To lead a likelihood Ratio Test : "LR_test" (intermediary function)

In this study, we need to test DIF (Differential Item Functioning). What is DIF ? DIF occurs for an item when people with the same ability but from different groups (groups defined by a covariate, such as gender or age) have unequal probabilities of giving a response. In other words, the conditional probability of success for a DIF-item depends on the value of a covariate. To apply Rasch models, we must have no DIF items. We then need to test that hypothesis when we apply a Rasch model. To check if an item is DIF for a specific covariate, we proceed as follows :

- from the DIF-suspected item and the covariate, we create new items : one for each level of the covariate, equal to the original item when the covariate values the studied level, and equal to NA otherwise. 
We take an example : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA ;

- we may then apply two Rasch models : one model to the database with the original item, and one model to the database with the new items we just created (and at the same time we exclude the original item) ;

- then, we may apply a LR-test : if the p-value is under a specified thereshold, then we can reject the null hypothesis : "there is no DIF". 

To this end, we wrote the "LR_test" function, which very basically apply a log-likelihood ratio test. This function takes in argument two models ("model0" and "model1"), which might be for example gpcm objects : the first model must be the one with the original item, and the second model the one with the items created from the covariate. This function has also a "display" argument : "display" must be a boolean, equal to "TRUE" if you wish to display the result in the R-Console. 
The function returns a dataframe including the value of the LR-test, the degrees of freedom and the p-value. 

This function will be used in several other functions and generally won't be used "on its own", that's why we qualify it as an "intermediary function". 



# To plot the expected value : "expected_value" (intermediary function)

In the case of polytomous items, we cannot plot the probability of success as a function of the latent variable eta, as we can in the case of dichotomous items : we do not interpret the dichotomous and polytomous items the same way. However, we may plot the expected value of an item for a gpcm-model applied to a database with polytomous items. 

To this end, we wrote the "expected_value" function. This functions takes in argument one database ("data"), one item ("item"), one itemlist ("itemlist") and one constraint ("constraint"). The database is the database we extract the information we need from. We restrain this database to an itemlist to apply a gpcm model ("gpcm" stands for "generalized partial credit model"), which might be a 1PL or 2PL, according to the constraint we choose : 

> if the constraint is "rasch", the discrimination parameter is assumed equal for all items and fixed to one ;

> if the constraint is "1PL", the discrimination parameter is assumed equal for all items but is estimated ;

> if the constraint is "gpcm", the discrimination parameters differ between the items and are estimated.

For more details on the gpcm objects, we recommand the reading of the description of the "ltm" R-package.

The item is the item whom we wish to calculate the expected value. The itemlist corresponds to the subscale the item belongs : when we apply the model, we restrain the studied dataset to this itemlist. The item and the itemlist's type must be "character". The constraint specifies the type of model we wish to apply, as we just detailed ; it may be equal to "rasch", "1PL" or "gpcm". 

Now we may describe the running of the function :

- first, we apply the model to the restrained database ; from this model, we extract the possible values of eta (the latent variable we're interested in) ;
- we then extract (still from the model), for every possible value taken by the polytomous item (for example 1, 2, 3, 4, 5), the probability that the item is equal to each value (for example, the probability that the item would be equal to 1, 2, etc) ; with all these probabilities, we may then calculate the expected value.

The function returns a vector of two columns : the value of eta and, for each one, the corresponding expected value. 

This function is an "intermediary function" : we especially use it in the "dif_poly" function to plot the expected value. 

# To regather some values of a variable when we wish to split this variable according to a differentiating variable : "regather" (intermediary function)


As we previously said (in the description of the LRT function), we have to test DIF on the items of the studied database. For this, we have to split each item according to a covariate, the differentiating variable. From the DIF-suspected item and the covariate, we create new items : one for each level of the covariate, equal to the original item when the covariate values the studied level, and equal to NA otherwise. When we do this, there is a risk that the newly-created items do not take all the same values. We take the same example again : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA.
Now, suppose item1 has 4 possible values : 1, 2, 3, 4. We might have to face the following situation : item1_1 taking all 4 possible values (1, 2, 3, 4) and item1_2 taking only 3 possibles values (2, 3, 4). This might happen if the value "1" is initially rare. 
This is problematic if we wish to plot and compare the expected values of the newly-created items. 

In order to avoid that, the "regather" function gathers the possible values of the studied item so that the items created from a differentiating variable do take the same values. The "regather" function is applied to a database, "data", an item, "item", and a differentiating variable, "diffvar" ; "item" and "diffvar" must be of "character" type ;

- first, we extract the levels of "diffvar" ; we check if "diffvar" includes at least  two levels (if not, "diffvar" has no interest) ;
- for each possible value of "diffvar", we list the values taken by the item when "diffvar" is equal to the studied value ; 
- we then check whether those lists includes elements which are not common to all lists ;
- if this last case is encountered, we then recode the item so that for each value of "diffvar", the item takes the same levels.

The "regather" function returns a dataframe, "data", which is either the original dataframe if all possible values were taken for each level of "diffvar" by "item", either the modified dataframe, with the recoded item. 

This "regather" function is an "intermediary function" : we use it in all functions when we need to split an item according to a differentiating variable.

# To create differenciated items from a covariate : "differenciate" (intermediary function)

To test DIF on an item splited according to a covariate, we need to create the splited items. We take an example : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA. 

The "differenciate" function takes in argument a dataframe, "data", an item, "item", and a differentiating variable, "diffvar" ; "item" and "diffvar" must be of "character" type. 

The function's running is quite intuitive ; one thing we might notice is that we first apply the "regather" function to the dataframe, the item and the differentiating variable. 

The "differenciate" function returns a list of two elements : a dataframe which includes the original dataframe and the newly-differenciated variables, and a vector of names, the names of the newly-differenciated variables. We need the two elements because this function is called in other functions and then it is very useful to have the names of the differenciated variables.

This "differenciate" function is an "intermediary function" : we use it in several function when we need to split an item according to a differenciating variable.

# To test DIF on an item and a covariate : "dif_poly_int" (intermediary function) and "dif_poly"

Now, we wish to test DIF on the items of a order polytomous items database, according to a covariate. For this, we use the "dif_poly_int" function, which takes in argument : 
> "data", the studied dataframe ;
> "item", the item whom we wish to test DIF according to the covariable ; "item" must be a character string ;
> "itemlist", the list of items on which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" belongs to "itemlist" ;
> "diffvar", the covariate from which we differenciate the item "item" ; "diffvar" must be a character string ;
> "constraint", which precise which model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;
> "toPlot", optional, default value set to "FALSE" : we may precise "toPlot = TRUE" to have in addition the plot of the expected values of the differenciated items.

We now describe the running of the function :

- "initializing - create differenciate variables" :
first, we apply the intermediary function "differenciate" to add to "data" the differenciated items created from "item" and "diffvar" ; we then create two models : "fit0" is the model applied to the itemlist (so with the original item "item") and "fit1" is the model applied to the itemlist without the original item, to which we add the newly-differenciated items ("difflist") ; we save the item parameters of "fit1" ; 
- "with or without plot" :
if "toPlot" is equal to TRUE, we calculate the expected value for the differenciated items, calling the "expected_value" function ; we then display the plots in the same graph ; in addition to that, we display the result of the log-likelihood ratio test by calling the "LR_test" function and we return the item parameters of the model with the differenciated items that we saved earler ; if "toPlot" is equal to FALSE, we return a list of two elements : the "LR_test" result and the item parameters we previously saved. 


Now, we generally wish to test DIF on a list of items, not on a single item. For this, we use a second function, "dif_poly". This function has the same arguments as the "dif_poly_int" function, except for "item", which is replaced by "items" : "items" must be a vector of items (a vector of character string). We then apply the dif_poly_int function for each item of the "items" list :

- if "toPlot" is equal to TRUE, the function will display all plots one by one and at the same time the LRT result and the parameters of the model with the differenciated items, for each item of the "items" list ; 
- if "toPlot" is equal to FALSE which is the by default value), the function will return a 


























