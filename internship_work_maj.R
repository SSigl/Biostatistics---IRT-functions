#===================================================================================#
#                                                                                   #
#              INTERNSHIP - MASTER I/II - BIOSTATISTICS - Personnal work            #
#                       Latent variables model - IRT - Rasch                        #
#                                     S. SIGL                                       #
#                                                                                   #
#===================================================================================#

#===================================================================================#
# the environnment

#%%%%%%%%%%%%%%

#   packages
#rm(list = ls())

#install.packages("eRm",repos="http://cran.rstudio.com/")
#library(eRm)

#install.packages("sirt", repos="http://cran.rstudio.com/") #pour tester dépendance locale
#library(sirt)

#install.packages("ltm", repos="http://cran.rstudio.com/")
#library(ltm)

#install.packages("mirt", repos="http://cran.rstudio.com/")
#library(mirt)


#===================================================================================#
# first database : dichotomous items

# function to lead a likelihood ratio test
LRtest = function(model0,model1,display){
  L0 <- logLik(model0)
  L1 <- logLik(model1)
  nb0 <- attr(L0, "df")
  nb1 <- attr(L1, "df")
  df <- nb1 - nb0
  LRT = -2*(L0[1]-L1[1])
  pvalue = round(1 - pchisq(LRT,df),digits=5)
  result = data.frame("df"=c(df),"pvalue"=c(pvalue))
  if(display==TRUE){
  text0 = "LR test :"
  text1 = paste("Log-likelihood ratio test",LRT,sep=":")
  text2 = paste("Degrees of freedom",df,sep=":")
  text3 = paste("p-value",pvalue,sep=":")
  LRtext = paste(text0,text1,text2,text3,sep="\n")
  cat("\n")
  cat(LRtext)
  cat("\n")}
  else{return(result)}
}

# NB : this function may only be used for gpcm objects
splitFct = function(data,item,itemlist,diffvar,constraint,toPlot){
  if(class(item)!="character"){
  item = colnames(data)[item]}
  # levels of the diff variable
  level = as.numeric(levels(as.factor(data[,c(diffvar)])))
  # number of levels
  n = length(level) 
  # number of columns of the complete database
  m = dim(data)[2] 
  # empty vector which will countain the names of the new variables for the suspected item by level
  var_names = c()
  
  # we create the names of the variables and the variables themselves
  for(i in 1:n){
    var_names[i] = paste(item,i,sep="_")
    data$var <- data[,c(item)]
    data$var[data[,c(diffvar)] != level[i]] <- NA 
    colnames(data)[colnames(data)=="var"] <- var_names[i]
  }
  
  # we run the "original" model, so we can compare likelihood ratio
  it0 <- data[,c(itemlist)]
  # we remove the suspected item from the itemlist (it is easier to proceed that way)
  k_1 = which(colnames(data)==item)
  k_1 = which(itemlist==k_1)
  itemlist = itemlist[-k_1]
  # we add the new created items : (m+1):(m+n)
  it1 <- data[,c(itemlist,(m+1):(m+n))]
  # models
  fit0 <- gpcm(it0,constraint = constraint)
  fit1 <- gpcm(it1,constraint=constraint)
  
  if(toPlot == TRUE){
  # we consider the variables whom the name starts with the name of the suspected item
  # i.e. all new created variables
  l = grep(paste(item,"",sep="_"), colnames(it1))
  
  # we plot 
  v1 <- plot(fit1, items= l, plot=FALSE,xlim=c(-4,4))
  eta = v1$z
  plot(eta,v1$pr[[var_names[1]]][,2],col=c(1),type="l",xlab="eta",ylab="P(y=1|eta)",ylim=c(0,1))
  for(i in 2:n){
    lines(eta,v1$pr[[var_names[i]]][,2],col=c(i))
  }
  legend("topleft",legend=c(var_names),col=c(1:n),lty=rep(1,n))
  titre = paste("DIF-test for the item",item,sep=" ")
  title(sub=titre)
  
  # we also want the likelihood ratio
  cat(paste("For the item",item,sep=" "))
  LRtest(fit0,fit1,display=TRUE)}else{
    return(LRtest(fit0,fit1,display=FALSE))
    
  }
}

# now we use this first function which only worked for one item, and we apply it for a list of suspected items
splitTot = function(data,items,itemlist,diffvar,constraint){
  len = length(items)
  par(mfrow=selectPar(len))
  for(i in 1:len){
    splitFct(data,item = items[i],itemlist,diffvar,constraint)
    if(i<len){
    invisible(readline(prompt="Press [enter] to continue"))}
  }
}

#splitTot(AMTS,items = c(1,4,7),itemlist=c(1:10),diffvar="agegrp",constraint="rasch")
#par(mfrow=c(1,1))

#===================================================================================#
# polytomous items

# returns a dataframe with the eta and expected value corresponding to a gpcm-fit object
esperance = function(data,item,itemlist,constraint){
  # first we want all the values taken by the polytomous item
  factIt <- factor(data[,item])
  values = levels(factIt)
  n = length(values)
  # model
  it <- data[,c(itemlist)]
  fit <- gpcm(it,constraint=constraint)
  v <- plot(fit, plot=FALSE)
  eta = v$z
  m = length(eta)
  # we create a dataframe which will countain all probability by value of the item
  proba = c()
  Probas= data.frame(var=rep(0,m))
  Probas$eta <- eta
  for(i in 1:n){
    proba[i] = paste("p",i,sep="_")
    Probas$var <- v$pr[[item]][,i]
    colnames(Probas)[colnames(Probas)=="var"] <- proba[i]}
  # now we want to calculate the expected value ('esperance' in french, so "esp")
  Probas$esp <- Probas[,paste("p",1,sep="_")]*as.numeric(values[1])
  for(j in 2:n){
    Probas$esp <- Probas$esp + Probas[,paste("p",j,sep="_")]*as.numeric(values[j])}
  # now we may plot it !
  #plot(eta,Probas$esp,col=c(1),type="l",xlab="eta",ylab="Espérance")
  return(Probas)
}

# this function takes in parameter a database, an item and a diffvar
# we split the item regarding the diffvar 
# it may happen that the item does not take every value for each level of the diffvar
# then we want to gather some values of the item to re-scale it
regroupe = function(data,item,diffvar){
  # levels of the diff variable (we need to "factorize" first)
  fact <- factor(data[,c(diffvar)])
  level = levels(fact) 
  n=length(level)
  # we create a list which will contains, for each level of the diffvar, the values
  # taken by the item we're interested in
  list_levels = list()
  # we extract all values taken by the interesting item for each level
  for(i in 1:n){
    data_sub = subset(data,data[,diffvar]==level[i])
    assign(paste("level",i,sep="_"),levels(factor(data_sub[,item])))
    list_levels = append(list_levels,list(get(paste("level",i,sep="_"))))
    #assign(paste("n",i,sep="_"),length(levels(factor(data_sub[,item]))))
  }
  # we now want to check if the vectors of the values of the item
  # for each level contain the same values
  # for this we use the setdiff function which extract the values which are 
  # in one vector and not in another
  # we check all lists of values 2 by 2
  same_level = c()
  for(i in 1:(n-1)){ # rq : we necessarily have n>1 or the diffvar has no interest
    l1 = get(paste("level",i,sep="_"))
    l2 = get(paste("level",i+1,sep="_"))
    l = intersect(l1,l2)
    bool = c(length(setdiff(l1,l)) > 0 | length(setdiff(l2,l) >0))
    same_level = c(same_level,bool)
  }
  # we have now a vector which says if the two succedding lists of items do not contain the same values
  # if these is one TRUE, then we have to recode
  if(TRUE %in% same_level){
    # we extract the values common to all the vector
    # this vector cannot be empty in our study
    # (if by chance it was, it may imply DIF ?)
    common = as.numeric(Reduce(intersect,list_levels))
    p = length(common)
    new_level = c(1:p) # we will later recode the final variable with this new values (to re-scale)
    # this is the recoded variable
    data$var <- data[,item]
    # we recode by  categorie
    data$var[data$var <= common[1]] <- common[1]
    for(j in 1:(p-1)){
      data$var[data$var > common[j] & data$var <= common[j+1]] <- common[j+1]
    }
    data$var[data$var > common[p]] <- common[p]
    # now we recode with the new values (from 1 to p)
    for(j in 1:p){
      data$var[data$var == common[j]] <- new_level[j]
    }
    # we delete the initial item and rename the new variable "var" we just created
    k = which(colnames(data)==item)
    data <- data[,-k]
    colnames(data)[colnames(data)=="var"] <- item
    }
  return(data)
     }



# same function as previously for split but for polytomous item
difPoly = function(data,item,itemlist,diffvar,constraint,toPlot){
  data = subset(data,is.na(data[,diffvar])==FALSE)
  data = regroupe(data,item,diffvar) 
  if(class(item)!="character"){
    item = colnames(data)[item]}
  # levels of the diff variable (we need to "factorize" first)
  fact <- factor(data[,c(diffvar)])
  level = levels(fact) 
  # number of levels
  n = length(level) 
  # number of columns of the complete database
  m = dim(data)[2] 
  # empty vector which will countain the names of the new variables for the suspected item by level
  var_names = c()
  
  # we create the names of the variables and the variables themselves
  for(i in 1:n){
    var_names[i] = paste(item,level[i],sep="_")
    data$var <- data[,c(item)]
    data$var[data[,c(diffvar)] != level[i]] <- NA 
    colnames(data)[colnames(data)=="var"] <- var_names[i]
  }
  
  # in order to perform a log-likelihood test : 
  it0 <- data[,c(itemlist)]  
  # we remove the suspected item from the itemlist (it is easier to proceed that way)
  if(class(itemlist)!="character"){
  k_1 = which(colnames(data)==item)
  k_1 = which(itemlist==k_1)
  itemlist = itemlist[-k_1]}else{
    k_1 = which(itemlist==item)
    itemlist = itemlist[-k_1]
  }
  # we add the new created items : (m+1):(m+n)
  if(class(itemlist)!="character"){
  it1 <- data[,c(itemlist,(m+1):(m+n))]}else{
    it1 <- data[,c(itemlist,colnames(data)[(m+1):(m+n)])]
  }  
  # we create the model (even if we do it in the esperance function, in order to make a LR test)
  fit0 = gpcm(it0,constraint=constraint)
  fit1 = gpcm(it1,constraint=constraint)
  
  # two options : either we want the p-value and the plot, either we only want the p-value
  if(toPlot==TRUE){
  # we "re-initialize" itemlist to apply the "esperance" function
  len = dim(it1)[2]
  itemlist <- c(1:len)
  # we consider the variables whom the name starts with the name of the suspected item
  # i.e. all new created variables
  l = grep(paste(item,"",sep="_"), colnames(it1))
  # we create several database "Probas" containing the expected values for every new created item
  for(j in l){
    k = n-(len - j)
    assign(paste("dataPr",k,sep="_"),esperance(data=it1,item=c(j),itemlist=itemlist,constraint))
  }
  y_sup = max(as.integer(levels(as.factor(data[,item]))))
  plot(dataPr_1[,"eta"],dataPr_1[,"esp"],col=c(1),type="l",xlab="eta",ylab="Expected value",ylim=c(0,y_sup))
  for(j in 2:n){
    data = get(paste("dataPr",j,sep="_"))
    lines(data[,"eta"],data[,"esp"],col=c(j))
  }
  legend("bottomright",legend=c(var_names),col=c(1:n),lty=rep(1,n))
  titre = paste("DIF-test for the item",item,sep=" ")
  title(sub=titre)
  
  # we also want the likelihood ratio
  cat(paste("For the item",item,sep=" "))
  LRtest(fit0,fit1,display=TRUE)}
  if(toPlot==FALSE){
    result = LRtest(fit0,fit1,display=FALSE)
    rownames(result)=item
    return(result) # if we only want the p-value
  }
}

# same as before, we now want to directly apply the function on a list of suspected items
difPolyTot = function(data,items,itemlist,diffvar,constraint,toPlot){
  len = length(items)
  if(toPlot==TRUE){
  par(mfrow=selectPar(len))
  for(i in 1:len){
    difPoly(data,item = items[i],itemlist,diffvar,constraint,toPlot)
    if(i<len){
      invisible(readline(prompt="Press [enter] to continue"))}
  }
  }else{
    result = data.frame("df"=c(0),"pvalue"=c(0))
    for(i in 1:len){
      result[i,] = difPoly(data,item = items[i],itemlist,diffvar,constraint,toPlot)
    }
    rownames(result)=items
    return(result)
    }
}

# if we want to export each plot automatically
difPolyTotAutoPlot = function(data,items,itemlist,diffvar,constraint,path){
  len = length(items)
  for(i in 1:len){
    mypath <- file.path(path,paste(paste(items[i], diffvar, sep = "_"),".jpg",sep=""))
    jpeg(file=mypath)
    difPoly(data,item = items[i],itemlist,diffvar,constraint)
    dev.off()
    #if(i<len){
    #  invisible(readline(prompt="Press [enter] to continue"))}
  }
}



#===================================================================================#
# simulation - confidence intervals

simulation_1 = function(data,item,itemlist,constraint,B,scoreGrp){
  if(class(item)!="character"){
    item = colnames(data)[item]}
  if(class(itemlist)!="character"){
  itemlist = colnames(data)[itemlist]}
  # first, we create the score column (for the itemlist)
  data$score <- apply(subset(data,select=itemlist),1,sum,na.rm = TRUE)
  # scoreGrp is an option to "gather" some score categories so we have
  # enough people in each category 
  # scoreGrp specifies how many categories we wish to create apart from the score
  if(scoreGrp >1){
  data$score <- cut(data$score,breaks=scoreGrp,labels=FALSE)
  }
  R = as.integer(levels(factor(data$score)))
  # for the item we're interested in, we are going to create a data.frame with,
  # for everey corresponding score, the sum for all persons of the value of the item
  # when the score is equal to the corresponding score
  tab = data.frame(level_R=R,S = rep(0,length(R)))
  # we go over every possible value of R
  for(i in 1:length(R)){
    # we considerer the restrained database of the initial data, for every score value equal to the value of R we consider
    data_r = subset(data,data$score == as.numeric(R[i]) & is.na(data[,item])==FALSE,select=item)
    tab[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
  }
  plot(tab$level_R,tab$S,type="l",xlab="score",main=paste("Score for the item",item,sep=" "),col="blue")
  # we note we do not obtain a monotonic function for each item
  
  # we calculate the gpcm-fit object corresponding to our data
  fit = gpcm(data[,c(itemlist)],constraint=constraint)
  coef = coef(fit)
  persons = data.frame(factor.scores(fit)$score.dat)
  # rq : "persons" and "data" are of different size because some lines are duplicated
  # there is a probleme if we consider that the persons who are repeated
  # should normally have a higher probability of being selected
  # so we have to replicate the lines if obs > 1
  N = dim(persons)[1]
  p = N
  for(i in 1:N){
    k = persons[i,"Obs"]
    if(k>1){
      for(j in 1:(k-1)){
        persons[p+1,] = persons[i,]
        p = p+1
        }
        }
  }
  N = dim(persons)[1]
  q = length(itemlist)

  
  # now we do B boostrap simulations using the gpcm fit object
  for(b in 1:B){
    # boostrap : we sample some rows
    data_B = persons[sample(nrow(persons),size=N,replace=TRUE),]
    data_B = subset(persons,select=c("z1"))
    # in order to avoid any index problem :
    rownames(data_B) <- NULL
    
    # we calculate, for each person location, the probability P(Xiv=x|theta) for x =0,1,2,... (possible levels of each item)
    for(k in 1:q){
      item_B = itemlist[k]
      values = as.integer(levels(factor(data[,item_B])))
      n = length(values)
      # parameters :
      parameters = as.data.frame(coef[item_B,])
      # we add the corresponding parameter for the lowest value's category (which is set to 0 - convention)
      parameters["Catgr.0",] <-0
      # we simply sort the dataframe by the name of the row (it may avoid technical problems in what follows)
      parameters <- data.frame(parameters[ sort(row.names(parameters)), ],row.names = sort(row.names(parameters)))
      colnames(parameters)[colnames(parameters)==colnames(parameters)] <- "value"
      # empty vector which will countain the names of the calculated probabilities for each possible value of the item
      var_names = c()
      
      #before, we calculate useful entities
      # discriminatino parameter
      alpha = parameters["Dscrmn",]
      # required parameters :
      eta = c=("eta"=c(0))
      for(j in 1:(n-1)){
        beta_j = parameters[paste("Catgr",j,sep="."),]
        eta=c(eta,eta[j]-beta_j)
        }
      quotient = 0
      for(j in 1:n){
        quotient = quotient + exp(alpha*(values[j]*data_B$z1+eta[j]))
      }
      quotient <- 1/quotient
      
      # we now create the names of the variables and the variables themselves
      for(j in 1:n){
        var_names[j] = paste("proba",item_B,values[j],sep="_")
        
        # we select the good parameters
        eta_j = eta[j]
        h=values[j]
        
        # we calculate the corresponding probability
        data_B$var <- exp(alpha*(h*data_B$z1+eta_j))*quotient
        colnames(data_B)[colnames(data_B)=="var"] <- var_names[j]
      }
    }
    
    # now we simulate an item-dataset, by using a multinomial law
    simul_data = data[1,itemlist]
    for(i in 1:N){
      for(k in 1:q){
        item_B = itemlist[k]
        values = as.integer(levels(factor(data[,item_B])))
        n = length(values)
        var_names = c()
        for(j in 1:n){
          var_names[j] = paste("proba",item_B,values[j],sep="_")
        }
        prob = data_B[i,var_names]
        simul_data[i,k] <- values[which(as.data.frame(rmultinom(n=1,size=1,prob))==1)]
      }
    }
    simul_data$score <- apply(simul_data,1,sum)  
    if(scoreGrp >1){
      simul_data$score <- cut(simul_data$score,breaks=scoreGrp,labels=FALSE)
    }
    
    # now we may as previously compute the mean for each level of the score
    R = as.integer(levels(factor(simul_data$score)))
    # for the item we're interested in, we are going to create a data.frame with,
    # for everey corresponding score, the sum for all persons of the value of the item
    # when the score is equal to the corresponding score
    tab_b = data.frame(level_R=R,S = rep(0,length(R)))
    # we go over every value of R
    for(i in 1:length(R)){
      # we considerer the restrained database of the initial data, for every score value equal to the value of R we consider
      data_r = subset(simul_data,simul_data$score == as.numeric(R[i]),select=item)
      tab_b[i,][2] <- sum(data_r[,item])/dim(data_r)[1]
    }
    # we may finally plot the result
    lines(tab_b$level_R,tab_b$S,lty = 3,col="red")
  }
  lines(tab$level_R,tab$S,col="blue",lwd=3)
}

# now, we want to add a split option according to a differenciation variable diffvar
# for that, we need a very simple differenciation function 

diff_fct = function(data,item,diffvar){
  data = subset(data,is.na(data[,diffvar])==FALSE)
  data = regroupe(data,item,diffvar) 
  if(class(item)!="character"){
    item = colnames(data)[item]}
  # levels of the diff variable (we need to "factorize" first)
  fact <- factor(data[,c(diffvar)])
  level = levels(fact) 
  # number of levels
  n = length(level) 
  # number of columns of the complete database
  m = dim(data)[2] 
  # empty vector which will countain the names of the new variables for the item we wish to differenciate by level
  var_names = c()
  
  # we create the names of the variables and the variables themselves
  for(i in 1:n){
    var_names[i] = paste(item,level[i],sep="_")
    data$var <- data[,c(item)]
    data$var[data[,c(diffvar)] != level[i]] <- NA 
    colnames(data)[colnames(data)=="var"] <- var_names[i]
  }
  
  result = list("data"=data,"var_names"=var_names)
  return(result)
}

# same as before, we now want to directly apply the function on a list of items
simulation_1_Tot = function(data,items,itemlist,constraint,B,scoreGrp)
  {
  len = length(items)
  par(mfrow=selectPar(len))
  for(i in 1:len){
    simulation_1(data=data,item=items[i],itemlist=itemlist,constraint=constraint,B=B,scoreGrp=scoreGrp)
  }
}

# now we want a diff option
simulation_1_diff = function(data,item,itemlist,constraint,B,scoreGrp,diffvar){
  result = diff_fct(data,item,diffvar)
  data = result$data
  var_names = result$var_names
  len = length(var_names)
  par(mfrow=selectPar(len+1))
  # first, we do the simulation for the splited item 
  simulation_1(data=data,item,itemlist=itemlist,constraint=constraint,B=B,scoreGrp=scoreGrp)
  # now, we change itemlist
  k = which(itemlist == item)
  itemlist = itemlist[-k]
  itemlist = c(itemlist, var_names)
  for(i in 1:len){
    simulation_1(data=data,item=var_names[i],itemlist=itemlist,constraint=constraint,B=B,scoreGrp=scoreGrp)
  }
}

simulation_1_diff_Tot = function(data,items,itemlist,constraint,B,scoreGrp,diffvar){
  len = length(items)
  for(i in 1:len){
    simulation_1_diff(data,items[i],itemlist,constraint,B,scoreGrp,diffvar)
    if(i<len){
      invisible(readline(prompt="Press [enter] to continue"))}
  }
  }
  
  


#===================================================================================#
# some useful functions 

# to select columns in a dataframe whom the name starts with expr
selectItems = function(data,expr){
  l = grep(paste(expr,"",sep=""), colnames(data))
  return(l)
}

# to delete the columns with only NA values in a dataset
selectNotNa = function(data){
  # to delete the items with only NA answers (not relevant - means there were not used for the study)
  data <- data[,colSums(is.na(data)) != nrow(data)]
  return(data)
}

selectPar = function(n){
  if(n <=3){
    return(c(1,n))}
  if(n == 4 ){
    return(c(2,2))
  }
  if(n >4 & n <= 6){
    return(c(2,3))  
  }
  if(n >6 & n <=9){
    return(c(3,3))  
  }
  if(n > 9){
    return(c(1,3))  
  }
}


#===================================================================================#
# to test how reliable an item is

# we write a function which will test, for every dataset in "dataList"
# and for every covariate in "covariate"
# if when we differenciate regarding the covariate, there is a DIF prb
# (we report each p-value)
# rq : here, we necessarily have "item" and "itemlist" from "character" class
# because there are several datasets
testItem = function(dataList,item,itemlist,covariates,constraint,eps){
  n = length(dataList)
  names = names(dataList)
  covar = as.character(covariates[,names[[1]]][1])
  study = names[[1]]
  resultItem = data.frame("covar"=c(covar),"study" = c(study),"pvalue"=c(0))
  for(i in 1:n){
    study = names[[i]]
    covar = levels(covariates[,study])
    m= length(covar)
    for(j in 1:m){
    diffvar = covar[j]
    #options(warn=-1) 
    #options(show.error.messages = -1)
    #result = try(difPoly(dataList[[study]],item=item,itemlist=itemlist,diffvar=diffvar,constraint=constraint,toPlot=FALSE),silent=TRUE)
    result = tryCatch(difPoly(dataList[[study]],item=item,itemlist=itemlist,diffvar=diffvar,constraint=constraint,toPlot=FALSE),error=function(e) NA)
    #if(class(result)!="try-error"){
    #data_int = data.frame("covar"=c(diffvar),"study"=c(study),"pvalue"=result[,"pvalue"])
    #}else{
    #data_int = data.frame("covar"=c(diffvar),"study"=c(study),"pvalue"="ERROR")  
    #}
    if(is.na(result)){
      data_int = data.frame("covar"=c(diffvar),"study"=c(study),"pvalue"="ERROR")    
    }else{
      data_int = data.frame("covar"=c(diffvar),"study"=c(study),"pvalue"=result[,"pvalue"])  
    }
    resultItem = rbind(resultItem,data_int)
  }}
  resultItem = resultItem[-1,]
  rownames(resultItem) <- NULL
  resultItem$test <- (resultItem$pvalue < eps)  
  resultItem$test[resultItem$pvalue=="ERROR"] <- "ERROR"
  return(resultItem)
}


testItem_Tot = function(dataList,items,itemlist,covariates,constraint,eps){
  result_Tot = list()
  len = length(items)
  for(i in 1:len){
    result = testItem(dataList,items[i],itemlist,covariates,constraint,eps)
    result_Tot = list.append(result_Tot,result)
  }
  names(result_Tot) <- items
  return(result_Tot)
}











