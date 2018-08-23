
# Biostatistics---IRT-functions

Here you may find the description of every R-functions written in the "internship_work" file. The purpose of this work is to facilitate the study of databases with polytomous item responses and to analyze them using IRT models. 

# The packages required :

We need several packages to correctly run the functions we are about to describe :
- the "ltm" package : fits IRT models for polytomous item ;
- the "car" package : useful to recode some variables ;
- the "rlist" package : useful to use the "list.append()" function ;
- the "gtools" package : useful to sort some lists in our programm.


# To compute a likelihood Ratio Test : "LR_test" (intermediary function)

In this study, we need to test DIF (Differential Item Functioning). What is DIF ? DIF occurs for an item when people with the same ability but from different groups (groups defined by a covariate, such as gender or age) have unequal probabilities of response. In other words, the conditional probability of success for a DIF-item depends on the value of a covariate. To apply IRT models, we must have no DIF items. We then need to test that hypothesis when we apply a IRT model. To check if an item is DIF for a specific covariate, we proceed as follows :

- from the DIF-suspected item and the covariate, we create new items : one for each level of the covariate, equal to the original item when the covariate values the studied level, and equal to NA otherwise. 
We take an example : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA ; here is an example :

```R
  item1 item2 Sex item1_1 item1_2
1     1     2   1       1      NA
2     1     1   1       1      NA
3     2     2   2      NA       2
4     3     4   1       3      NA
5     2     3   2      NA       2
6     4     2   1       4      NA

```

- we may then apply two IRT models : one model to the database with the original item, and one model to the database with the new items we just created (and at the same time we exclude the original item) ;

- then, we may apply a LR-test : if the p-value is under a specified thereshold, then we can reject the null hypothesis : "there is no DIF". 

To this end, we wrote the "LR_test" function, which very basically apply a log-likelihood ratio test. This function takes in argument two models ("model0" and "model1"), which might be for example gpcm objects : the first model must be the one with the original item, and the second model the one with the items created from the covariate. This function has also a "display" argument : "display" must be a boolean, equal to "TRUE" if you wish to display the result in the R-Console. 
The function returns a dataframe including the value of the LR-test, the degrees of freedom, the value of the chi-square and the p-value. 

This function will be used in several other functions and generally won't be used "on its own", that's why we qualify it as an "intermediary function". 



# To plot the expected value : "expected_value" (intermediary function)

In the case of polytomous items, we cannot plot the probability of success as a function of the latent variable eta, as we can in the case of dichotomous items : we do not interpret the dichotomous and polytomous items the same way. However, we may plot the expected value of an item for a gpcm-model applied to a database with polytomous items. 

To this end, we wrote the "expected_value" function. This functions takes in argument one database ("data"), one item ("item"), one itemlist ("itemlist") and one constraint ("constraint"). The database is the database we extract the information we need from. We restrain this database to an itemlist to apply a gpcm model ("gpcm" stands for "generalized partial credit model"), which might be a 1PL or 2PL, according to the constraint we choose : 

> if the constraint is "rasch", the discrimination parameter is assumed equal for all items and fixed to one ;

> if the constraint is "1PL", the discrimination parameter is assumed equal for all items but is estimated ;

> if the constraint is "gpcm", the discrimination parameters differ between the items and are estimated.

For more details on gpcm objects, we refer to the description of the "ltm" R-package.

The item is the item for which we wish to calculate the expected value. The itemlist corresponds to the subscale the item belongs : when we apply the model, we restrain the studied dataset to this itemlist. The item and the itemlist's type must be "character". The constraint specifies the type of model we wish to apply, as we just detailed ; it may be equal to "rasch", "1PL" or "gpcm". 

Now we may describe the running of the function :

- first, we apply the model to the restrained database and we extract the possible values of eta (the latent variable we're interested in) ;
- we then extract (still from the model), for every possible value of the polytomous item (for example 1, 2, 3, 4, 5), the probability that the item is equal to each value (for example, the probability that the item would be equal to 1, 2, etc) ; with all these probabilities, we may then calculate the expected value.

The function returns a vector of two columns : the value of eta and, for each one, the corresponding expected value of the item. 

This function is an "intermediary function" : we especially use it in the "dif_poly" function to plot the expected value. 

# To collapse some values of a variable when we wish to split this variable according to a differentiating variable : "collapse" (intermediary function)


As we previously said (in the description of the LRT function), we have to test DIF on the items of the studied database. For this, we have to split each item according to a covariate, the differentiating variable. From the DIF-suspected item and the covariate, we create new items : one for each level of the covariate, equal to the original item when the covariate values the studied level, and equal to NA otherwise. When we do this, there is a risk that the newly-created items do not take all the same values. We take the same example again : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA.
Now, suppose item1 has 4 possible values : 1, 2, 3, 4. We might have to face the following situation : item1_1 taking all 4 possible values (1, 2, 3, 4) and item1_2 taking only 3 possibles values (2, 3, 4). This might happen if the value "1" rarely occurs. Here is a simple example : 

```R
  item1 item2 Sex item1_1 item1_2
1     1     2   1       1      NA
2     1     1   1       1      NA
3     1     2   2      NA       1
4     3     4   1       3      NA
5     2     3   2      NA       2
6     4     2   1       4      NA
```
You may see that for "Sex = 1", "item1" takes the values "1, 3, 4" ("item1_1") ; and for "Sex=2", "item1" takes the values "1, 2" ("item1_2") ; we then need to collapse the values of "item1" so that the splitted items "item1_1" and "item1_1" take the same values.

This is problematic if we wish to plot and compare the expected values of the newly-created items. 

In order to avoid that, the "collapse" function groups the possible values of the studied item so that all collapsed items take the same values. The "collapse" function is applied to a database, "data", an item, "item", and a differentiating variable, "diffvar" ; "item" and "diffvar" must be of "character" type ;

- first, we extract the levels of "diffvar" ; we check if "diffvar" includes at least  two levels (if not, "diffvar" has no interest) ;
- for each possible value of "diffvar", we list the values taken by the item when "diffvar" is equal to the studied value ; 
- we then check whether those lists includes elements which are not common to all lists ;
- if this last case is encountered, we then recode the item so that for each value of "diffvar", the item takes the same levels.

The "collapse" function returns a dataframe, "data", which is either the original dataframe if all possible values were taken for each level of "diffvar" by "item", either the modified dataframe, with the recoded item. 

This "collapse" function is an "intermediary function" : we use it in all functions when we need to split an item according to a differentiating variable.

# To create splitted items from a covariate : "split" (intermediary function)

To test DIF on an item with respect to a covariate, we need to create the splitted items. We take an example : suppose we suspect item1 to be a DIF for the covariate "Sex". "Sex" takes only two values : 1 (men) and 2 (women). We then create the variables "item1_1" (for men) and "item1_2" (for women). The variable "item1_1" is equal to the variable "item1" when "Sex" is equal to 1 ("men"), otherwise "item1_1" is equal to NA ; similarly, the variable "item1_2" is equal to the variable "item1" when "Sex" is equal to 2 ("women") ; otherwise "item1_2" is equal to NA. 

The "split" function takes in argument a dataframe, "data", an item, "item", and a differentiating variable, "diffvar" ; "item" and "diffvar" must be of "character" type. 

The function's running is quite intuitive ; one thing to notice is that we first apply the "collapse" function to the dataframe, the item and the differentiating variable. 

The "split" function returns a list of two elements : a dataframe which includes the original dataframe and the newly-splitted variables, and a vector of names, the names of the newly-splitted variables. We need the two elements because this function is called in other functions and then it is very useful to have the names of the splitted variables. Here is an example of the dataframe we could obtain :


```R
  item1 item2 Sex item1_1 item1_2
1     1     2   1       1      NA
2     1     1   1       1      NA
3     2     2   2      NA       2
4     3     4   1       3      NA
5     2     3   2      NA       2
6     4     2   1       4      NA

```

This "split" function is an "intermediary function" : we use it in several function when we need to split an item according to a differenciating variable.

# To test DIF on an item and a covariate : "dif_poly_int" (intermediary function) and "dif_poly"

Now, we wish to test DIF on the items of a order polytomous items database, according to a covariate. For this, we use the "dif_poly_int" function, which takes in argument : 
> "data", the studied dataframe ;

> "item", the item whom we wish to test DIF according to the covariable ; "item" must be a character string ;

> "itemlist", the list of items to which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" belongs to "itemlist" ;

> "diffvar", the covariate from which we split the item "item" ; "diffvar" must be a character string ;

> "constraint", choosing model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "toPlot", optional, default value set to "FALSE" : we may precise "toPlot = TRUE" to have in addition the plot of the expected values of the splitted items.

We now describe the running of the function :

- "initializing - create split variables" :
first, we apply the intermediary function "split" to add to "data" the splitted items created from "item" and "diffvar" ; we then create two models : "fit0" is the model applied to the itemlist (so with the original item "item") and "fit1" is the model applied to the itemlist without the original item, to which we add the newly-splitted items ("difflist") ; we save the item parameters of "fit1" in a dataframe ; 
- "with or without plot" :
if "toPlot" is equal to TRUE, we calculate the expected value for the splitted items, calling the "expected_value" function ; we then display the plots in the same graph ; in addition to that, we display the result of the log-likelihood ratio test by calling the "LR_test" function and we return the item parameters of the model with the splitted items that we saved earler ; if "toPlot" is equal to FALSE, we return a list of two elements : the "LR_test" result and the item parameters we previously saved. 


Now, we generally wish to test DIF for a list of items. For this, we use a second function, "dif_poly". This function has the same arguments as the "dif_poly_int" function, except for "item", which is replaced by "items" : "items" must be a vector of items (a vector of character strings). We then apply the dif_poly_int function for each item of the "items" list :

- if "toPlot" is equal to TRUE, the function will display all plots one by one and, at the same time,for each item of the "items" list, the LRT result and the parameters of the model with the splitted items ; 
- if "toPlot" is equal to FALSE (which is the by default value), the function will return a list of two elements : a dataframe including the results of the LRT test for all items in "items" list and a list of dataframes, one for each item in the "items" list, including the parameters of the model applied on the splitted items. 

The reader should note that it might happen that the dif_poly_int function returns an error because of convergence's problems of the model. This is a problem if we wish to apply the "dif_poly" function on a list of several items because it might interrupt the process : to avoid this, we use "tryCatch" on "dif_poly_int" so that if the "dif_poly_int" returns no result, the loop does pass to the next item in the "items" list. 

This "dif_poly" function may be used for only one item in the "items" argument, so we will not use the "dif_poly_int" function, which is only an intermediary function useful in the "dif_poly" function, and we will only use "dif_poly".


# To check whether the model fits our data : "simulation_1_int" (intermediary function) and "simulation_1"

When we apply a gpcm model on a polytomous dataframe, we wish to perform simulations in order to check model fit. For this, we wrote several functions. First, the most basic ones we use are "simulation_1_int" and "simulation_1". The idea is the following : for an item in a subscale, we plot the mean as a function of the observed total score obtained by the test-takers ; we apply the gpcm-model and we extract the persons parameters and the parameters of the studied item ; we then sample the persons locations (bootstrap) and from the model parameters, we simulate a polytomous responses dataset ; we may then calculate the simulated score and plot it in the same graph as the first plot. We repeat the sampling several times to see, for the studied item, how similar the simulated means are from the observed means. 

The intermediary function "simulation_1_int" takes the arguments :

> the studied dataframe "data" ;

> "item", the item whom we shall trace the mean according to the total score ; let's take an example : we wish to know the mean of item n°3 depending on the total score of the test-takers, which is between 0 and 10 ; then, for each possible value of the score (0, 1, 2, etc), we calculate the mean of the item n°3 on all persons who obtained this possible value of the score ; so, we first restrain the database to the persons who had a total score of 0 and we calculate the mean of the item n°3, then we do the same with the persons who had a total score of 1, etc ; we then obtain for each possible value of the score the corresponding mean of the item, and plot this item as a function of the score ; "item" must be a character string ;

> "itemlist", the list of items on which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" belongs to "itemlist" ;

> "constraint", choosing the model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "B" for "bootstrap", which is the number of simulations we wish to use ;

> "sc_gp" for "score group" : this is an optional parameter, by default set to 1 ; "sc_gp" indicates in how many levels we need to recode the total score, to improve the plot's readibility ; indeed, for example if we study a subscale of 6 items between 1 and 5, the score might goes from 0 (because of missing values) to 30 ; then, we might want to recode this score variable so it has only 3, 4 or 5 levels for instance. 

We now describe the running of the functions :

- "initialiazing -- creation of score -- re-level score" :

first, we restrain the dataset on the subscale defined by "itemlist" and we calculate the corresponding score by summing the result for each item by row ; when "sc_gp" is higher than 1 (which is the by default value), we calculate the quantiles of the score distribution so to split this distribution in a number of "sc_gp" equal parts, and we may then recode the score variable on this cut ; we then have the levels of the score variable ;

- "calculate observed mean" : 

for each level of the score variable, we restrain the dataset to the rows for which the score is equal to the regarded and for this selection, we calculate the mean of the item we are intereted in ;

- "plot observed mean" :

we then plot the mean of the item we previously calculated depending on the score ;

- "model" :

we apply the model on the dataset restrained to the list of items "itemlist" (we must not include the score column in the model) ; we save the items' parameters in "coef" and the persons' parameters in "persons" ; there is a subtlety here : the persons' parameters is a dataframe which includes all possible rows of answers originally in the dataset, and for each row, the corresponding person location and the number of time this row is observed in the dataset ; for example, if we consider the subscale item1, item2, item3, we suppose that 3 people have answered (1,2,1) and 2 people have answered (3,3,4) : then, in the persons' parameters, those anwers will appear only one time, and the "Obs" column precises that the first answer appeared 3 times and that the second answer appeared 2 times ; we later wish to sample the persons' parameters and so, we need to take the number of observations into account ; for this, we calculate for every row its frequency of occurrence (in other words, we weight each row) ;

- "Bootstrap" : now, we may do the simulations strictly speaking ; we repeat the following description a number "B" of times :

> first, we sample the persons' parameters and we only save the persons' location

> "calculate for each person location the probability P(Xiv=x|theta) for x every possible levels of each item" : for each persons' location we previsouly sampled, we calculate, for each level of each item of the subscale, the probability that this item is equal to this level conditionnaly to this persons location ; we apply the link-function of the rasch model ;

> "simulated dataset" : now that we have all the probabilities for each value of all items and for each persons' location, we may simulate a polytomous dataset using a multinomial law ; this way, by using a rasch model, we obtain a simulated dataset which is the same size of the original dataset ;

> "simulated score" : we may now calculate the score of the simulated dataset and, as before, recode this score in "sc_gp" levels, using the same cuts as the ones used for recoding the observed score ; then, we do as before and we calculate the mean of the items for each level of the score ; finally, we might plot this simulated mean in the same graph as the original plot.


The "simulation_1_int" function displays the plot of the observed and simulated means for only one item ; it is an intermediary function. To display the results of this "simulation_1_int" function for several items, we use the "simulation_1" function. 

The "simulation_1" function has the same arguments as the "simulation_1_int" function, except for "item" which is replaced by "items" and must be a vector of the items for which we wish to lead the simulation ; "items" must be a vector of string characters ; note that all items in "items" must also be in "itemlist". The "simulation_1" function successively applys the "simulation_1_int" function for all items in "items". 

This "simulation_1" function may be used for only one item in the "items" argument, so we will not use the "simulation_1_int" function, which is only an intermediary function useful in the "simulation_1" function, and we will only use "simulation_1". 

# To compute a simulation with splitted variables : "simulation_2_int" (intermediary function) and "simulation_2"

Now, we might wish to compute the simulation we previously described not on the original dataset but on the dataset with splitted items. For this, we use two functions : "simulation_2" and "simulation_2_int". The "simulation_2_int" function is only an intermediary function. 

The "simulation_2" function takes the following arguments : 

> "data", the regarded dataset ; "data" must be a data.frame ;

> "item", the item we wish to split ; "item" must be a character string ; 

> "itemlist", the list of items to which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" belongs to "itemlist" ;

> "constraint", choosing the model ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "B" for "bootstrap", which is the number of simulations we wish to lead ;

> "sc_gp" for "score group" : this is an optional parameter, by default set to 1 ; "sc_gp" indicates in how many levels we need to recode the total score, to improve the plot's readibility ; indeed, for example if we study a subscale of 6 items between 1 and 5, the score might goes from 0 (because of missing values) to 30 ; then, we might want to recode this score variable so it has only 3, 4 or 5 levels for instance ;

> "diffvar", the covariate according which we wish to split the item "item" in the dataset "data" ; "diffvar" must be a character string ;

> "unq" for "unique" : this is an optional parameter, by default set to "TRUE" ; "unq" indicates whether the "simulation_2" function is applied only once or in a more general function ; this parameter is useful to select the appropriate par for the plot's display ;

The "simulation_2" function's running is the following : 

- first, we split the regarded item using the "split" function ; we obtain a new dataset with the splitted items and also the names of the newly-splitted items, that we save in "diff_names" ; 

- then, we apply the "simulation_1" function to the dataset restrained to the itemlist and regarding the item "item" (we do not need the splitted items there) ; 

- after this, we apply a second function, "simulation_2_int", to the dataset with the splitted items for the items in "diff_names", i.e. the splitted items.


We now present this second function, "simulation_2_int" ; this function has the following arguments : 

> "data", the regarded dataset, which contains the splitted items we previously created using the "split" function ;

> "item", the item for which we wish to plot the observed and simulated means ; this item must be a splitted item (this is why the "simulation_2_int" function is in "simulation_2" function successively applied to the variables in "diff_names") ; 

> "itemlist", the list of items which are not splitted and that we need to use in the IRT model ; in relation to the "itemlist" parameter of the "simulation_2" function, the itemlist we consider here is the previous itemlist to which we substract the initial item "item" (since we only wish to keep the splitted items, and not the original one) ; we take an example to make things clearer : assume we wish to split "item1" according to a two-levels covariate and we consider the itemlist ("item1, item2, item3") ; the splitted items are "item1_1" and "item1_2" ; we then naturally want to apply the rasch model on ("item1_1, item1_2, item_2, item_3") and not ("item1, item1_1, item1_2, item_2, item_3") ; this is why we do not consider the same itemlist in the two functions ; 

> "diff_names", the names of the splitted variables ;

> "constraint", choosing the model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "B" for "bootstrap", which is the number of simulations we wish to lead ;

> "sc_gp" for "score group" : this is an optional parameter, by default set to 1 ; "sc_gp" indicates in how many levels we need to recode the total score, to improve the plot's readibility ; indeed, for example if we study a subscale of 6 items between 1 and 5, the score might goes from 0 (because of missing values) to 30 ; then, we might want to recode this score variable so it has only 3, 4 or 5 levels for instance. 

Concerning the running of the function, it basically follows the same steps as in the "simulation_1_int" function but there are a few subtleties : 

- there are no changes in the calculation of the score, except that we have to be carefull to the objects we manipulate ("tot_list" and not "itemlist" for example) ;

- for the simulation strictly speaking, there is one major difference : concerning the sampling, we have to sample a dataset for each level of the DIF variable ; initially, we successively restrain the persons parameter dataset for each level of the DIF variable and sample this dataset (bootstrap) ; then, we apply the same steps as in the "simulation_1_int" function for each level of the differenciating variable ; finally, we merge all simulated datasets and obtain one simulated dataset and then plot the simulated mean of the regarded splitted items.


The "simulation_2_int" function displays the plot of the observed and simulated means for one splitted item ; it is an intermediary function. The "simulation_2" function displays the plot of the observed and simulated means for one item and the splitted items obtained from this original item ; all the plots are displayed side by side so that these might be compared. 

If we wish to apply the whole operation to several items, we might use the "simulation_2_tot" function, which will simply successively apply the "simulation_2" function to a vector of items, "items". 


# To better compare the results of simulations : "simulation_3"

The "simulation_3" functions re-use the previous simulation we wrote ; we use this "simulation_3" function if we want to use specific simulations and to display all their results side by side. The "simulation_3" has the same arguments as the "simulation_2" function. The difference between "simulation_2" and "simulation_3" is that "simulation_2" will only lead the simulation for the original item "item" and the splitted items created from "item" and the differenciating variable "diffvar", whereas "simulation_3" will in addition lead the simulation for the items in itemlist which are not the item "item". Then, the "simulation_3" function will display the result of the simulation for all items in itemlist and for the splitted items. 

# To plot simulations and means of splitted variables in the same graph - option to extract plot : "simulation_4_int" (intermediary function) and "simulation_4"

This simulation function may be used to plot in the same graph the observed and simulated means of an item, as well as the observed means of splitted items created from the item and a differenciating variable. These functions do not present any substantial changes in relation to the previous functions we wrote : the "simulation_4_int" function does follow exactly the same steps as the "simulation_2_int" ans "simulation_2" function, but it was necessary to write a new function to obtain the desired result because plot-objects in R are difficult to modify once they are plot, so we needed to directly work on the plot, and so to write another function. 

The "simulation_4_int" function has the same arguments as the "simulation_2" function and an additional optional parameter, "display", which is by default equal to "TRUE". By default, the "simulation_4_int" function returns a single graph with the observed mean of an original item, the simulated mean of this original item, and the observed means of splitted items created from the original item and a differentiating variable ; if display is set to "FALSE", the "simulation_4_int" function returns a dataframe with the coordinates of all the plots (the observed means of the original and split items plus the simulated mean of the original item) ; rather than plotting them, it allows you to export those plots.

If we wish to repeat "simulation_4_int" function on several items instead of only one, we use the "simulation_4" function which calls the "simulation_4_int" function for each item in "items", which is a vector of character strings. By default, the "simulation_4" function plots all graphs side by side for each item in "items" ; if "display" is set to "FALSE", the "simulation_4" function returns a list of dataframes, one for each item in "items", each dataframe containing the coordinates of all the regarded plot. 

This "simulation_4" function may be used for only one item in the "items" argument, so we will not use the "simulation_4_int" function, which is only an intermediary function useful in the "simulation_4" function, and we will only use "simulation_4".



# To apply a simulation to a database with one splitted item : "simulation_5_int" (intermediary function) and "simulation_5"

Assume you face the following case : you try to visually identify DIF and observe that two items of the same subscale seem suspect. You may wish to re-apply the simulation on one of the suspected item by differenciating the other suspected item and applying the model to the subscale with these splitted items. For this, we needed to write two new functions : "simulation_5_int", which is an intermediary function, and "simulation_5".

The "simulation_5" function takes the following arguments :

> "data", the regarded dataset ; "data" must be a data.frame ;

> "item", the item we wish to split and whom we wish to plot the mean for the splitted items ; "item" must be a character string ; 

> "item_dif", the second item we wish to split but whom we do not wish to plot the mean ;

> "itemlist", the list of items to which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" and "item_dif" belong to "itemlist" ;

> "constraint", choosing the model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "B" for "bootstrap", the number of simulations ;

> "sc_gp" for "score group" : this is an optional parameter, by default set to 1 ; "sc_gp" indicates in how many levels we need to recode the total score, to improve the plot's readibility ; indeed, for example if we study a subscale of 6 items between 1 and 5, the score might goes from 0 (because of missing values) to 30 ; then, we might want to recode this score variable so it has only 3, 4 or 5 levels for instance ;

> "diffvar", the covariate according which we wish to split the item "item" and the item "item_dif" in the dataset "data" ; "diffvar" must be a character string ;

> "unq" for "unique" : this is an optional parameter, by default set to "TRUE" ; "unq" indicates whether the "simulation_2" function is applied only once or in a more general function ; this parameter is useful to select the appropriate par for the plot's display.

The "simulation_5" function works like this : 

- first, we split the items "item_dif" and "item" using the "split" function ; we obtain a new dataset with the splitted items and also the names of the newly-splitted items, that we respectively save in "dif_list" and "to_dif_list" ; we also save the levels of the covariate "diff_var" in "level" ;

- then, we apply the "simulation_2_int" function to the dataset restrained to the items that are not splitted, the item "item" and the splitted items of "dif_list", regarding the item "item" ; in other words, we do the simulation for the original item "item", but we have to apply the "simulation_2_int" function because of the other items from "dif_list" which are splitted ; indeed, simulating a dataset with splitted items requires additive steps compared to a dataset with only non-splitted items ;

- after this, we apply a second function, "simulation_5_int", to the dataset with the splitted items for the items in "dif_list", i.e. the splitted items of the item we do not wish to lead the simulation for, and the splitted items for the items in "to_dif_list", i.e. the splitted items of the item we do wish the lead the simulation for.

We now present this second function, "simulation_5_int" ; this function has the following arguments : 

> "data", the regarded dataset, which contains the splitted items we previously created using the "split" function from dif_list and from to_dif_list (which means, the splitted items for the two items we are initially interested in);

> "item", the item for which we wish to plot the observed and simulated means ; this item must be one of the splitted item from "to_dif_list" (this is why the "simulation_5_int" function is in "simulation_5" function successively applied to the variables in "to_dif_list") ; 

> "itemlist", the list of items which are not splitted and that we need to include in the IRT model ; in relation to the "itemlist" parameter of the "simulation_5" function, the itemlist we consider here is the previous itemlist to which we substract the initial items "item" and "item_dif" (since in the IRT model, we only wish to keep the splitted items, and not the original one) ; we take an example to make things clearer : assume we have to split "item1" and "item2" according to a two-levels covariate and we consider the itemlist ("item1, item2, item3") ; the splitted items are ("item1_1, item1_2") and ("item2_1, item2_2") : we wish to lead the simulation for the splitted items from "item2", including splitted items from "item1" ; we then naturally want to apply the rasch model on ("item1_1, item1_2, item2_1, item2_2, item3") and not ("item1, item1_1, item1_2, item2, item2_1, item2_2, item3") ; this is why we do not consider the same itemlist in the two functions ; 

> "dif_list", the names of the splitted variables for "item_dif", that we need to apply the rasch model, and for which we do not wish to lead the simulations for ;

> "to_dif_list", the names of the splitted variables for "item", that we need to apply the rasch model, and for which we do wish to lead the simulation for ; 

> "level", the levels of the covariate "diffvar" which was initially used to split the database ; 

> "constraint", choosing the model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "B" for "bootstrap", which is the number of simulations we wish to lead ;

> "sc_gp" for "score group" : this is an optional parameter, by default set to 1 ; "sc_gp" indicates in how many levels we need to recode the total score, to improve the plot's readibility ; indeed, for example if we study a subscale of 6 items between 1 and 5, the score might goes from 0 (because of missing values) to 30 ; then, we might want to recode this score variable so it has only 3, 4 or 5 levels for instance. 

The function basically follows the same steps as "simulation_1_int" but there are a few subtleties : 

- there are no changes in the calculation of the score, except that we have to be carefull to the objects we manipulate ("tot_list" and not "itemlist" for example) ;

- for the simulation strictly speaking, there is one major difference : concerning the sampling, we have to sample a dataset for each level of the differentiating variable ; initially, we successively have to restrain the persons parameters' dataset for each level of the differenciating variable and sample this dataset (bootstrap) ; then, we apply the same steps as in the "simulation_1_int" function for each level of the differenciating variable ; finally, we may merge all the simulated dataset we obtain to have only one simulated dataset and then plot the simulated mean of the regarded splitted items ; we may note that the fact that there is now two splitted items instead of one does not introduce major changes compared to the "simulation_2_int" function.


The "simulation_5_int" function displays the plot of the observed and simulated means for one splitted item from "to_dif_list" ; the model applied and used for the simulation is the model with all splitted items (from dif_list and to_dif_list) ; "simulation_5_int" is an intermediary function. The "simulation_5" function displays the plot of the observed and simulated means for one item and the splitted items obtained from this original item ; all the plots are displayed side by side so that these might be compared. 


# A global function : "simulation"

To ease the manipulation of those functions, we wrote a global "simulation" function : depending on the arguments you provide, you will obtain a different result.

This "simulation" function may take several arguments : 

> "data", the studied dataset ; "data" must be a dataframe ;

> "item", the item for which we wish to plot the mean and lead the simulation ; "item" must be a character string ;

> "itemlist", the list of items on which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" belongs to "itemlist" ;

> "constraint", which precise which model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "B" for "bootstrap", which is the number of simulations we wish to lead ;

> "sc_gp" for "score group" : this is an optional parameter, by default set to 1 ; "sc_gp" indicates in how many levels we need to recode the total score, to improve the plot's readibility ; indeed, for example if we study a subscale of 6 items between 1 and 5, the score might goes from 0 (because of missing values) to 30 ; then, we might want to recode this score variable so it has only 3, 4 or 5 levels for instance ;

> "diffvar" is an optional parameter ; it precises the covariate according which we wish to split the item "item" in the dataset "data" ; "diffvar" must be a character string ;

> "item_dif" is an optional parameter ; it precises that we wish to split two items and apply a simulation to the splitted item "item" by taking in account another splitted item "item_dif" in the model ; note that if you wish to apply that kind of simulation, "diffvar" must necessarily be provided (or if it is not, the function will simply lead a classic simulation without taking the "item_dif" in account) ;

> "samePlot" is an optional parameter ; it is a boolean, by default set to "FALSE".

Now, we precise which arguments we need to provide for which type of simulation :

- if we want to lead the "classic" simulation, which means apply the "simulation_1" function, we do not need to fill out "diffvar" nor "samePlot" ;

- if we want to split the regarded item, we need to provide "diffvar" ; then, there are three options :

> if "dif_list" is also provided, the simulation lead is applied to the splitted items from "item" once you already splitted the item "item_dif" ; the regarded model is applied to the splitted items obtained from "item" and "item_dif" ; the carried out simulation will then be "simulation_5" ;

> we might want to lead the simulation for every item in itemlist and for the splitted items ; then, we do not need to provide any additionnal parameter ; the carried out simulation will then be "simulation_3" ;

> or, we might want to plot in the same graph the observed mean of an original item, the simulated mean of this original item, and the observed means of splitted items created from the original item and a differentiating variable ; then, we must provide "samePlot = TRUE" ; the carried out simulation will then be "simulation_4" ;

There is no option to lead "simulation_2" because "simulation_3" provides more information. 


# To globally test an item's dif over several datasets : "test_item_int" (intermediary function) and "test_item"

The "dif_poly" function allows us to check DIF for several items and one covariate in one database. Now, suppose we have several databases with the same variables (i.e. items) and for each database several covariate. For one item, we wish to have a global overview of whether DIF is observed in each database and for each covariate. The "test_item_int" and "test_item" functions provide this global overview.

The "test_item_int" function takes the following arguments :

> "dataList", a list containaing all the regarded databases and their respective names ; "dataList" must be a list and their elements need to have a name ; 

> "item", the item for which we wish to test DIF according to the covariate in all databases and for all covariates ; "item" must be a character string ;

> "itemlist", the list of items to which we will restrain the dataframe to apply the gpcm model ; "itemlist" must be a vector of character strings ; note that necessarily, "item" belongs to "itemlist" ;

> "covariates", the list of covariates according to which we wish to test DIF for item "item", for each database in "dataList" ; "covariates" must be a list of same length of "dataList" ; the elements of covariates must have the same names as "dataList" ;

> "constraint", which precise which model we wish to apply ; "constraint" might be equal to "rasch", "1PL" or "gpcm" (c.f. more details in the description of the "expected_value" function) ;

> "eps" stands for "epsilon" ; it indicates the threshold for the log-likelihood ratio test we lead for each covariate in each database : if the p-value is smaller than this threshold then we consider that DIF is indeed observed.

The "test_item_int" function simply applies the "dif_poly" function on the item, for each covariate of each dataset ; it returns for each case the result of the lig-likelihood ratio test.

The "test_item_int" function is an intermediary function : in general, we wish to test dif not only on one item but on several items. For this, we use the "test_item" function : this function has the same arguments as the "test_item_int" function, except for "item", which is replaced by "items" : "items" must be a vector of items (a vector of character strings), all included in the vector "itemlist. We then apply the "test_item_int" function for each item of the "items" vector. Finally, "test_item" returns a list, as long as the number of items in "items" : for each item, the corresponding element of the list is a vector containing for each covariate and each database, the values of the log-likelihood ratio test, the chi-square, the p-value, and a "test" column whether the p-value is smaller than "eps", the provided threshold. Then, the result gives a global overview of the dif-character of the item in all the databases. 

This "test_item" function may be used for only one item in the "items" argument, so we will not use the "test_item_int" function, which is only an intermediary function useful in the "test_item" function, and we will only use "test_item".

# To sample a number of rows according to a variable's levels : "select_random"

Assume the following situation : you have several studies with the same items but not the same number of surveyed people. You may work on the total dataset with all results combined but the differences in cross sections's sizes might have an influence on the results. Then, you may wish to sample, for each study, a same number of people, and then work on the total database with same-sized cross sections. 

The "select_random" function takes the following arguments : 

> "data", the regarded dataset ; "data" must be a dataframe ;

> "refvar", the variable according which we will sample the database : for this, we will create one database for each level of "refvar" and sample the same number of rows in each newly-created database ; finally, we will merge all the sampled databases ;

> "n", the number of lines we wish to sample in each database by level ; "n" must be smaller than the minimal size of the database for each category ; 

The "select_random" function returns a dataframe which merges all sampled databases (one for each level of "refvar").



# To select the right par to display several plots side by side : "select_par" (intermediary function)

This function is used many times in the functions we previously described. It is useful when we have to display plots side by side and to still be able to see it all at the same time.
The argument of "select_par" is "n", the number of displayed plots.

# To remove the rows which only have NA answers : "select_not_na" (intermediary function)

This function which takes a dataset "data" in argument removes in this dataset the lines with only "NA" answers, and returns the modified dataset.
















