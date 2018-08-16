#===================================================================================#
#                                                                                   #
#              INTERNSHIP - MASTER I/II - BIOSTATISTICS - Personal work            #
#                       Latent variables model - IRT - Rasch                        #
#                                     S. SIGL                                       #
#                                                                                   #
#===================================================================================#

#===================================================================================#
# the environnment

#%%%%%%%%%%%%%%
#   packages
rm(list = ls())

#install.packages("eRm",repos="http://cran.rstudio.com/")
library(eRm)

#install.packages("sirt", repos="http://cran.rstudio.com/") #pour tester dépendance locale
library(sirt)

#install.packages("ltm", repos="http://cran.rstudio.com/")
library(ltm)

#install.packages("mirt", repos="http://cran.rstudio.com/")
library(mirt)

#install.packages("car") # useful to recode
library(car)

#install.packages("rlist") # to use the list.append() function
library(rlist)

#install.packages("gtools") # useful to sort some list
library(gtools)
#===================================================================================#

#===================================================================================#
# function to lead a likelihood ratio test
LR_test = function(model0,model1,display){
  L0 <- logLik(model0)
  L1 <- logLik(model1)
  nb0 <- attr(L0, "df")
  nb1 <- attr(L1, "df")
  df <- nb1 - nb0
  LRT = -2*(L0[1]-L1[1])
  pvalue = round(1 - pchisq(LRT,df),digits=3)
  LRT = round(LRT,digits=3)
  ChiSQ = round(pchisq(LRT,df),digits=3)
  result = data.frame("LRT" =c(LRT), "df"=c(df), "chisq"=ChiSQ, "pvalue"=c(pvalue))
  if(display==TRUE){
    text0 = "LR test :"
    text1 = paste("Log-likelihood ratio test",LRT,sep=":")
    text2 = paste("Degrees of freedom",df,sep=":")
    text3 = paste("Chi square",ChiSQ,sep=":")
    text4 = paste("p-value",pvalue,sep=":")
    LRtext = paste(text0,text1,text2,text3,text4,sep="\n")
    cat("\n")
    cat(LRtext)
    cat("\n")}
  else{return(result)}
  }
#===================================================================================#



#===================================================================================#
expected_value = function(data,item,itemlist,constraint){
  # levels of the polytomous item
  values = levels(as.factor(data[,item]))
  n = length(values)
  # application of the model
  fit <- gpcm(data[,c(itemlist)],constraint=constraint)
  v <- plot(fit, plot=FALSE)
  # we calculate the expected value
  df_exp = data.frame("eta"=v$z)
  df_exp$exp <- 0
  for(i in 1:n){
    df_exp$exp <- df_exp$exp + v$pr[[item]][,i]*as.numeric(values[i])
  }
  return(df_exp)
}
#===================================================================================#



#===================================================================================#
regather = function(data,item,diffvar){
  # levels of the diffvar 
  level = levels(as.factor(data[,c(diffvar)]))
  n = length(level)
  stopifnot(n>1)
  # for each level of "level", all values taken by the item 
  list_levels = list()
  for(i in 1:n){
    list_levels =append(list_levels,list(levels(as.factor(data[which(data[,diffvar]==level[i]),][,item]))))
  }
  same_level = c()
  # levels met for each level of the covariate
  common = as.numeric(Reduce(intersect,list_levels))
  stopifnot(is.null(common)==FALSE)
  for(i in 1:n){ 
    same_level = c(same_level,c(length(setdiff(list_levels[[i]],common))>0))
  }
  if(TRUE %in% same_level){
    p = length(common)
    data$var <- data[,item] # recoded variable
    data$var[data$var <= common[1]] <- common[1] # recoding
    for(j in 1:(p-1)){
      data$var[data$var > common[j] & data$var <= common[j+1]] <- common[j+1]
    }
    data$var[data$var > common[p]] <- common[p]
    for(j in 1:p){
      data$var[data$var == common[j]] <- j # recoding the new values from 1 to p
    }
    k = which(colnames(data)==item)
    data <- data[,-k]
    colnames(data)[colnames(datpa)=="var"] <- item
  }
  return(data)
}
#===================================================================================#


#===================================================================================#
differenciate = function(data,item,diffvar){
  # initializing
  data = subset(data,is.na(data[,diffvar])==FALSE)
  data = regather(data,item,diffvar) 
  level = levels(as.factor(data[,c(diffvar)])) 
  n = length(level) 
  diff_names = paste(item,level,sep="_")
  # create differenciated variables
  for(i in 1:n){
    data$var <- data[,c(item)]
    data$var[data[,c(diffvar)] != level[i]] <- NA 
    colnames(data)[colnames(data)=="var"] <- diff_names[i]
  }
  result = list("data"=data,"diff_names"=diff_names)
  return(result)
}
#===================================================================================#



#===================================================================================#
dif_poly_int = function(data,item,itemlist,diffvar,constraint,toPlot=FALSE){
  # initializing - creating differenciate variables
  result = differenciate(data,item,diffvar)
  data = result$data
  var_names = result$diff_names
  n=length(var_names)
  # model with original item
  fit0 = gpcm(data[,c(itemlist)],constraint=constraint)
  # model with differenciated items
  difflist = c(itemlist[-which(itemlist==item)],var_names)
  fit1 = gpcm(data[,c(difflist)],constraint=constraint)
  parameters = as.data.frame(fit1$coefficients)
  
  # with or without plot
  if(toPlot==TRUE){
    # expected value for differenciated items
    for(i in 1:n){
      assign(paste("data_pr",i,sep="_"),expected_value(data=data,item=var_names[i],itemlist=difflist,constraint))
    }
    y_sup = max(as.integer(levels(as.factor(data[,item]))))
    y_inf = min(as.integer(levels(as.factor(data[,item]))))
    # plot 
    plot(data_pr_1[,"eta"],data_pr_1[,"exp"],col=c(1),type="l",xlab="eta",ylab="Expected value",ylim=c(y_inf,y_sup))
    for(i in 2:n){
      data = get(paste("data_pr",i,sep="_"))
      lines(data[,"eta"],data[,"exp"],col=c(i))
    }
    legend("bottomright",legend=c(var_names),col=c(1:n),lty=rep(1,n),cex=0.6)
    title(paste("DIF-test for the item",item,sep=" "))
    # LRT
    cat(paste("For the item",item,sep=" "))
    LR_test(fit0,fit1,display=TRUE)
    print(parameters)
  }else{
    LRT = LR_test(fit0,fit1,display=FALSE)
    rownames(LRT)=item
    result=list("LRT"=LRT,"parameters"=parameters)
    return(result) 
  }
}

dif_poly = function(data,items,itemlist,diffvar,constraint,toPlot=FALSE){
  len = length(items)
  if(toPlot){
    par(mfrow=select_par(len))
    options(warn=-1) 
    for(i in 1:len){
      tryCatch(dif_poly_int(data,item=items[i],itemlist=itemlist,diffvar=diffvar,constraint=constraint,toPlot=TRUE),error=function(e) NA)
      if(i<len){
        invisible(readline(prompt="Press [enter] to continue"))}
    }
    options(warn=1) 
  }else{
    result_LRT = data.frame("LRT"=c(0),"df"=c(0),"chisq"=c(0),"pvalue"=c(0))
    result_parameters = list()
    options(warn=-1) 
    for(i in 1:len){
      result_dif = tryCatch(dif_poly_int(data,item=items[i],itemlist=itemlist,diffvar=diffvar,constraint=constraint,toPlot=FALSE),error=function(e) NA)   
      if(is.na(result_dif)==FALSE){
        result_LRT[i,] = result_dif[["LRT"]]
        result_parameters = append(result_parameters,list(result_dif[["parameters"]]))
      }
    }
    options(warn=-1) 
    names(result_parameters) = items
    rownames(result_LRT)=items
    result = list("LRT"=result_LRT,"parameters"=result_parameters)
    return(result)
  }
}
#===================================================================================#



#===================================================================================#
simulation_1_int = function(data,item,itemlist,constraint,B,sc_gp=1){
  # initialiazing -- creation of score -- re-level score
  data <- data[,c(itemlist)] 
  data$score <- apply(data,1,sum,na.rm = TRUE)
  if(sc_gp >1){
    quantiles = as.numeric(quantile(data$score,probs=0:sc_gp/sc_gp))
    data$score <- quantcut(data$score,q=sc_gp,na.rm=TRUE,labels=1:sc_gp)
  }
  R = as.integer(levels(factor(data$score)))
  tab = data.frame(level_R=R,S = rep(0,length(R)))
  # calculate real sucess rate
  for(i in 1:length(R)){
    data_r = subset(data,data$score == as.numeric(R[i]) & is.na(data[,item])==FALSE,select=item)
    tab[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
  }
  # plot real sucess rate
  y_inf = min(data[,itemlist],na.rm=TRUE)
  y_sup = max(data[,itemlist],na.rm=TRUE)
  if(sc_gp>1){x_lab="regathered score"}else{x_lab = "score"}
  plot(tab$level_R,tab$S,type="l",xlab=x_lab,ylab="sucess rate",main=paste("Sucess rate for item",item,sep=" "),col="blue",ylim=c(y_inf,y_sup))
  # model
  fit = gpcm(data[,c(itemlist)],constraint=constraint)
  coef = coef(fit)
  persons = data.frame(factor.scores(fit)$score.dat)
  n_pers = sum(persons[,"Obs"])
  n_item = length(itemlist)
  persons[,"Obs"] <- persons[,"Obs"]/n_pers
  
  # Bootstrap
  for(b in 1:B){
    data_B = persons[sample(nrow(persons),size=n_pers,replace=TRUE, prob=persons[,"Obs"]),]
    data_B = subset(data_B,select=c("z1"))
    rownames(data_B) <- NULL
    # calculate for each person location the probability P(Xiv=x|theta) for x every possible levels of each item
    for(k in 1:n_item){
      item_B = itemlist[k]
      values = as.integer(levels(factor(data[,item_B])))
      n_values = length(values)
      # extract and calculate useful entities : discrimination, eta, beta, quotient
      parameters = as.data.frame(coef[item_B,])
      parameters["Catgr.0",] <-0
      parameters <- data.frame(parameters[ sort(row.names(parameters)), ],row.names = sort(row.names(parameters)))
      colnames(parameters)[colnames(parameters)==colnames(parameters)] <- "value"
      alpha = parameters["Dscrmn",]
      eta = c=("eta"=c(0))
      for(j in 1:(n_values-1)){
        beta_j = parameters[paste("Catgr",j,sep="."),]
        eta=c(eta,eta[j]-beta_j)
      }
      quotient = 0
      for(j in 1:n_values){
        quotient = quotient + exp(alpha*(values[j]*data_B$z1+eta[j]))
      }
      quotient <- 1/quotient
      var_names = paste("proba",item_B,values,sep="_")
      for(j in 1:n_values){
        data_B$var <- exp(alpha*(values[j]*data_B$z1+eta[j]))*quotient
        colnames(data_B)[colnames(data_B)=="var"] <- var_names[j]
      }
    }
    # simulated dataset
    simul_data = data[1,itemlist]
    for(i in 1:n_pers){
      for(k in 1:n_item){
        item_B = itemlist[k]
        values = as.integer(levels(factor(data[,item_B])))
        var_names = paste("proba",item_B,values,sep="_")
        simul_data[i,k] <- values[which(as.data.frame(rmultinom(n=1,size=1,data_B[i,var_names]))==1)]
      }
    }
    # simulated score
    simul_data$score <- apply(simul_data,1,sum)
    if(sc_gp >1){
      simul_data$score <- cut(simul_data$score,breaks=quantiles,labels=FALSE)
    }
    R = as.integer(levels(factor(simul_data$score)))
    tab_b = data.frame(level_R=R,S = rep(0,length(R)))
    for(i in 1:length(R)){
      data_r = subset(simul_data,simul_data$score == as.numeric(R[i]),select=item)
      tab_b[i,][2] <- sum(data_r[,item])/dim(data_r)[1]
    }
    lines(tab_b$level_R,tab_b$S,lty = 3,col="red")
  }
  lines(tab$level_R,tab$S,col="blue",lwd=2)
  legend("bottomright",legend=c("real success rate","sim success rate"),col=c("blue","red"),lty=c(1,3),cex=0.6)
}

simulation_1 = function(data,items,itemlist,constraint,B,sc_gp=1){  
  len = length(items)
  par(mfrow=select_par(len))
  for(i in 1:len){
    simulation_1_int(data=data,item=items[i],itemlist=itemlist,constraint=constraint,B=B,sc_gp=sc_gp)
  }
}
#===================================================================================#



#===================================================================================#
simulation_2_int = function(data,item,itemlist,diff_names,constraint,B,sc_gp=1){
  # initializing - creation of score
  tot_list <- c(itemlist,diff_names)  
  data = subset(data,select=c(tot_list))
  data$score <- apply(data,1,sum,na.rm = TRUE)
  if(sc_gp >1){
    quantiles = as.numeric(quantile(data$score,probs=0:sc_gp/sc_gp))
    data$score <- quantcut(data$score,q=sc_gp,na.rm=TRUE,labels=1:sc_gp)
  }
  R = as.integer(levels(factor(data$score)))
  tab = data.frame(level_R=R,S = rep(0,length(R)))
  # calculate real sucess rate
  for(i in 1:length(R)){
    data_r = subset(data,data$score == as.numeric(R[i]) & is.na(data[,item])==FALSE,select=item)
    tab[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
  }
  # plot real sucess rate
  y_inf = min(data[,tot_list],na.rm=TRUE)
  y_sup = max(data[, tot_list],na.rm=TRUE)
  if(sc_gp>1){x_lab="regathered score"}else{x_lab = "score"}
  plot(tab$level_R,tab$S,type="l",xlab=x_lab,ylab="sucess rate",main=paste("Sucess rate for item",item,sep=" "),col="blue",ylim=c(y_inf,y_sup))
  # model
  fit = gpcm(data[,c(tot_list)],constraint=constraint)
  coef = coef(fit)
  persons = data.frame(factor.scores(fit)$score.dat)
  n_item = length(itemlist) + 1
  
  # Bootstrap
  for(b in 1:B){
    diff_n = length(diff_names)
    data_names = paste("data_B",diff_names,sep="_") # sampled databases
    simul_names = paste("simul",diff_names,sep="_") # simulated databases
    for(i in 1:diff_n){
      data_B = subset(persons,is.na(persons[,diff_names[i]])==FALSE)
      n_pers = sum(data_B[,"Obs"])
      data_B[,"Obs"] <- data_B[,"Obs"]/n_pers
      data_B = data_B[sample(nrow(data_B),size=n_pers,replace=TRUE,prob=data_B[,"Obs"]),]
      data_B = subset(data_B,select=c("z1"))
      rownames(data_B) <- NULL	
      assign(data_names[i], data_B) 
    }
    # calculate for each person location the probability P(Xiv=x|theta) for x every possible levels of each item
    for(l in 1:diff_n){
      for(k in 1:n_item){
        tot_list <- c(itemlist,diff_names[l])
        item_B = tot_list[k]
        values = as.integer(levels(factor(data[,item_B])))
        n_values = length(values)
        # extract and calculate useful entities : discrimination, eta, beta, quotient
        parameters = as.data.frame(coef[item_B,])
        parameters["Catgr.0",] <-0
        parameters <- data.frame(parameters[ sort(row.names(parameters)), ],row.names = sort(row.names(parameters)))
        colnames(parameters)[colnames(parameters)==colnames(parameters)] <- "value"
        alpha = parameters["Dscrmn",]
        eta = c=("eta"=c(0))
        for(j in 1:(n_values-1)){
          beta_j = parameters[paste("Catgr",j,sep="."),]
          eta=c(eta,eta[j]-beta_j)
        }
        quotient = 0
        for(j in 1:n_values){
          quotient = quotient + exp(alpha*(values[j]*get(data_names[l])[,"z1"]+eta[j]))
        }
        quotient <- 1/quotient
        var_names = paste("proba",item_B,values,sep="_")
        data_B = get(data_names[l]) # "transition" database
        for(j in 1:n_values){
          data_B$var <- exp(alpha*(values[j]*data_B$z1+eta[j]))*quotient
          colnames(data_B)[colnames(data_B)=="var"] <- var_names[j]
        }
        assign(data_names[l],data_B)
      }
    }
    # simulated datasets
    for(l in 1:diff_n){
      tot_list <- c(itemlist,diff_names[l])
      data_B = get(data_names[l])
      simul_data = data[1,tot_list]
      n_pers= dim(data_B)[1]
      for(i in 1:n_pers){
        for(k in 1:n_item){
          item_B = tot_list[k]
          values = as.integer(levels(factor(data[,item_B])))
          var_names = paste("proba",item_B,values,sep="_")
          simul_data[i,k] <- values[which(as.data.frame(rmultinom(n=1,size=1,data_B[i,var_names]))==1)]
        }
      }
      for(m in 1:diff_n){
        if(m !=l){
          simul_data$var <- NA
          colnames(simul_data)[colnames(simul_data)=="var"] <- diff_names[m]
          assign(simul_names[l],simul_data)
        }
      }
    }
    
    # final simulated dataset
    simul_data = get(simul_names[1])
    for(l in 2:diff_n){
      simul_data <- rbind(simul_data,get(simul_names[l]))
    }
    simul_data$score <- apply(simul_data,1,sum,na.rm=TRUE)
    if(sc_gp >1){
      simul_data$score <- cut(simul_data$score,breaks=quantiles,labels=FALSE)
    }
    R = as.integer(levels(factor(simul_data$score)))
    tab_b = data.frame(level_R=R,S = rep(0,length(R)))
    for(i in 1:length(R)){
      data_r = subset(simul_data,simul_data$score == as.numeric(R[i]) & is.na(simul_data[,item])==FALSE,select=item)
      tab_b[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
    }
    lines(tab_b$level_R,tab_b$S,lty = 3,col="red")
  }
  lines(tab$level_R,tab$S,col="blue",lwd=2)
  legend("bottomright",legend=c("real success rate","sim success rate"),col=c("blue","red"),lty=c(1,3),cex=0.6)
}


simulation_2 = function(data,item,itemlist,constraint,B,sc_gp=1,diffvar,unq = TRUE){
  result = differenciate(data,item,diffvar)
  data = result$data
  diff_names = result$diff_names
  len = length(diff_names)
  if(unq){
  par(mfrow=select_par(len+1))}
  # simulation for the original item
  simulation_1(data=data,item,itemlist=itemlist,constraint=constraint,B=B,sc_gp=sc_gp)
  # simulation for the differenciated items
  k = which(itemlist == item)
  itemlist = itemlist[-k]
  for(i in 1:len){
    simulation_2_int(data=data,item=diff_names[i],itemlist=itemlist,diff_names=diff_names,constraint=constraint,B=B,sc_gp=sc_gp)
  }
}

simulation_2_tot = function(data,items,itemlist,constraint,B,sc_gp=1,diffvar){
  len_1 = length(items)
  len_2 = length(levels(as.factor(data[,diffvar])))
  par(mfrow=select_par(len_2+1))
  for(i in 1:len_1){
    simulation_2(data,items[i],itemlist,constraint,B,sc_gp,diffvar)
    if(i<len_1){
      invisible(readline(prompt="Press [enter] to continue"))}
  }
}

#===================================================================================#



#===================================================================================#
simulation_3 = function(data,item,itemlist,constraint,B,sc_gp=1,diffvar){
  n = length(itemlist)
  len = length(levels(as.factor(data[,diffvar]))) + n
  par(mfrow=select_par(len))
  simulation_2(data,item,itemlist,constraint,B,sc_gp,diffvar,unq=FALSE)
  for(i in 1:n){
    if(itemlist[i]!=item){
      simulation_1_int(data=data,item=itemlist[i],itemlist,constraint=constraint,B=B,sc_gp=sc_gp)
    }
  }
}
#===================================================================================#




#===================================================================================#
simulation_4_int = function(data,item,itemlist,constraint,B,sc_gp=1,diffvar){
  # initializing - creation of score
  result = differenciate(data,item,diffvar)
  data = result$data
  diff_names = result$diff_names
  len = length(diff_names)
  all_items = c(item,diff_names) # all studied items 
  data = subset(data,select=c(itemlist,diff_names))
  difflist <- c(itemlist,diff_names) 
  k = which(difflist == item)
  difflist = difflist[-k] # without original item
  
  data$score <- apply(subset(data,select=itemlist),1,sum,na.rm = TRUE)
  if(sc_gp >1){
    quantiles = as.numeric(quantile(data$score,probs=0:sc_gp/sc_gp))
    data$score <- quantcut(data$score,q=sc_gp,na.rm=TRUE,labels=1:sc_gp)
  }
  R = as.integer(levels(factor(data$score)))
  tab = data.frame(level_R=R,S = rep(0,length(R)))
  for(i in 1:length(R)){
    data_r = subset(data,data$score == as.numeric(R[i]) & is.na(data[,item])==FALSE,select=item)
    tab[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
  }
  
  # plot real score
  y_inf = min(data[, itemlist],na.rm=TRUE)
  y_sup = max(data[, itemlist],na.rm=TRUE)
  if(sc_gp>1){x_lab="regathered score"}else{x_lab = "score"}
  plot(tab$level_R,tab$S,type="l",xlab=x_lab,ylab="sucess rate",main=paste("Simulation and differenciation for the",item,sep=" "),col="blue",ylim=c(y_inf,y_sup))
  
  # model
  fit = gpcm(data[,c(itemlist)],constraint=constraint)
  coef = coef(fit)
  persons = data.frame(factor.scores(fit)$score.dat)
  n_pers = sum(persons[,"Obs"])
  n_item = length(itemlist)
  persons[,"Obs"] <- persons[,"Obs"]/n_pers
  
  
  # Bootstrap
  for(b in 1:B){
    data_B = persons[sample(nrow(persons),size=n_pers,replace=TRUE, prob=persons[,"Obs"]),]
    data_B = subset(data_B,select=c("z1"))
    rownames(data_B) <- NULL
    
    # calculate for each person location the probability P(Xiv=x|theta) for x every possible levels of each item
    for(k in 1:n_item){
      item_B = itemlist[k]
      values = as.integer(levels(factor(data[,item_B])))
      n_values = length(values)
      # extract and calculate useful entities : discrimination, eta, beta, quotient
      parameters = as.data.frame(coef[item_B,])
      parameters["Catgr.0",] <-0
      parameters <- data.frame(parameters[ sort(row.names(parameters)), ],row.names = sort(row.names(parameters)))
      colnames(parameters)[colnames(parameters)==colnames(parameters)] <- "value"
      alpha = parameters["Dscrmn",]
      eta = c=("eta"=c(0))
      for(j in 1:(n_values-1)){
        beta_j = parameters[paste("Catgr",j,sep="."),]
        eta=c(eta,eta[j]-beta_j)
      }
      quotient = 0
      for(j in 1:n_values){
        quotient = quotient + exp(alpha*(values[j]*data_B$z1+eta[j]))
      }
      quotient <- 1/quotient
      var_names = paste("proba",item_B,values,sep="_")
      for(j in 1:n_values){
        data_B$var <- exp(alpha*(values[j]*data_B$z1+eta[j]))*quotient
        colnames(data_B)[colnames(data_B)=="var"] <- var_names[j]
      }
    }
    # simulated dataset
    simul_data = data[1,itemlist]
    for(i in 1:n_pers){
      for(k in 1:n_item){
        item_B = itemlist[k]
        values = as.integer(levels(factor(data[,item_B])))
        var_names = paste("proba",item_B,values,sep="_")
        simul_data[i,k] <- values[which(as.data.frame(rmultinom(n=1,size=1,data_B[i,var_names]))==1)]
      }
    }
    # simulated score
    simul_data$score <- apply(simul_data,1,sum)
    if(sc_gp >1){
      simul_data$score <- cut(simul_data$score,breaks=quantiles,labels=FALSE)
    }
    R = as.integer(levels(factor(simul_data$score)))
    tab_b = data.frame(level_R=R,S = rep(0,length(R)))
    for(i in 1:length(R)){
      data_r = subset(simul_data,simul_data$score == as.numeric(R[i]),select=item)
      tab_b[i,][2] <- sum(data_r[,item])/dim(data_r)[1]
    }
    lines(tab_b$level_R,tab_b$S,lty = 3,col="red")
  }
  lines(tab$level_R,tab$S,col="blue",lwd=2)
  
  # sucess rate for the differenciated items
  R = as.integer(levels(factor(data$score)))
  # list of colors to plot the differenciated items
  if(len<10){
    list_colors = c("green3","yellow3","pink","orange","purple","cyan","magenta","gray","aquamarine","coral")
  }else{
    list_colors = sample(colors(),3*len)
    for(name in c("blue","red")){
      if(name %in% list_colors){
        k = as.integer(which(list_colors==name))
        list_colors = list_colors[-k]}
    }
  }
  for(l in 1:len){
    item = diff_names[l]
    tab_b = data.frame(level_R=R,S = rep(0,length(R)))
    for(i in 1:length(R)){
      data_r = subset(data,data$score == as.numeric(R[i]) & is.na(data[,item])==FALSE,select=item)
      tab_b[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
    }
    lines(tab_b$level_R,tab_b$S,lty = 1,col=list_colors[l])
  }
  legend = paste(c("real rate for item","sim rate for item"),all_items[1],sep=" ")
  legend = c(legend,paste("real rate for item",diff_names,sep=" "))
  list_colors = c("blue","red",list_colors[1:len])
  legend("bottomright",legend=legend,col=list_colors,lty=c(1,3,rep(1,len)),cex=0.5)
}

simulation_4 = function(data,items,itemlist,constraint,B,sc_gp=1,diffvar){
  len = length(items)
  par(mfrow=select_par(len))
  for(i in 1:len){
    simulation_4_int(data,items[i],itemlist,constraint,B,sc_gp,diffvar)
    if(i<len & len > 9){
      invisible(readline(prompt="Press [enter] to continue"))}
  }
}
#===================================================================================#




#===================================================================================#
simulation_5_int = function(data,item,itemlist,dif_list,to_dif_list,level,constraint,B,sc_gp){
  # initializing - creation of score
  tot_list <- c(itemlist,dif_list,to_dif_list)  
  data = subset(data,select=c(tot_list))
  data$score <- apply(data,1,sum,na.rm = TRUE)
  if(sc_gp >1){
    quantiles = as.numeric(quantile(data$score,probs=0:sc_gp/sc_gp))
    data$score <- quantcut(data$score,q=sc_gp,na.rm=TRUE,labels=1:sc_gp)
  }
  R = as.integer(levels(factor(data$score)))
  tab = data.frame(level_R=R,S = rep(0,length(R)))
  # calculate real sucess rate
  for(i in 1:length(R)){
    data_r = subset(data,data$score == as.numeric(R[i]) & is.na(data[,item])==FALSE,select=item)
    tab[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
  }
  # plot real sucess rate
  y_inf = min(data[,tot_list],na.rm=TRUE)
  y_sup = max(data[, tot_list],na.rm=TRUE)
  if(sc_gp>1){x_lab="regathered score"}else{x_lab = "score"}
  plot(tab$level_R,tab$S,type="l",xlab=x_lab,ylab="sucess rate",main=paste("Sucess rate for item",item,sep=" "),col="blue",ylim=c(y_inf,y_sup))
  # model
  fit = gpcm(data[,c(tot_list)],constraint=constraint)
  coef = coef(fit)
  persons = data.frame(factor.scores(fit)$score.dat)
  n_item = length(itemlist) + 2
  n_dif = length(to_dif_list)
  
  # Bootstrap
  for(b in 1:B){
    data_names = paste("data_B",level,sep="_") # sampled databases
    simul_names = paste("simul",level,sep="_") # simulated databases    
    for(i in 1:n_dif){
      data_B = subset(persons,is.na(persons[,dif_list[i]])==FALSE)
      n_pers = sum(data_B[,"Obs"])
      data_B[,"Obs"] <- data_B[,"Obs"]/n_pers
      data_B = data_B[sample(nrow(data_B),size=n_pers,replace=TRUE,prob=data_B[,"Obs"]),]
      data_B = subset(data_B,select=c("z1"))
      rownames(data_B) <- NULL	
      assign(data_names[i], data_B) 
    }
    # calculate for each person location the probability P(Xiv=x|theta) for x every possible levels of each item
    for(l in 1:n_dif){
      for(k in 1:n_item){
        tot_list <- c(itemlist,dif_list[l],to_dif_list[l])
        item_B = tot_list[k]
        values = as.integer(levels(factor(data[,item_B])))
        n_values = length(values)
        # extract and calculate useful entities : discrimination, eta, beta, quotient
        parameters = as.data.frame(coef[item_B,])
        parameters["Catgr.0",] <-0
        parameters <- data.frame(parameters[ sort(row.names(parameters)), ],row.names = sort(row.names(parameters)))
        colnames(parameters)[colnames(parameters)==colnames(parameters)] <- "value"
        alpha = parameters["Dscrmn",]
        eta = c=("eta"=c(0))
        for(j in 1:(n_values-1)){
          beta_j = parameters[paste("Catgr",j,sep="."),]
          eta=c(eta,eta[j]-beta_j)
        }
        quotient = 0
        for(j in 1:n_values){
          quotient = quotient + exp(alpha*(values[j]*get(data_names[l])[,"z1"]+eta[j]))
        }
        quotient <- 1/quotient
        var_names = paste("proba",item_B,values,sep="_")
        data_B = get(data_names[l]) # "transition" database
        for(j in 1:n_values){
          data_B$var <- exp(alpha*(values[j]*data_B$z1+eta[j]))*quotient
          colnames(data_B)[colnames(data_B)=="var"] <- var_names[j]
        }
        assign(data_names[l],data_B)
      }
    }
    
    # simulated datasets
    for(l in 1:n_dif){
      tot_list <- c(itemlist,dif_list[l],to_dif_list[l])
      data_B = get(data_names[l])
      simul_data = data[1,tot_list]
      n_pers= dim(data_B)[1]
      for(i in 1:n_pers){
        for(k in 1:n_item){
          item_B = tot_list[k]
          values = as.integer(levels(factor(data[,item_B])))
          var_names = paste("proba",item_B,values,sep="_")
          simul_data[i,k] <- values[which(as.data.frame(rmultinom(n=1,size=1,data_B[i,var_names]))==1)]
        }
      } 
      for(m in 1:n_dif){
        if(m !=l){
          simul_data$var <- NA
          colnames(simul_data)[colnames(simul_data)=="var"] <- to_dif_list[m]
          simul_data$var <- NA
          colnames(simul_data)[colnames(simul_data)=="var"] <- dif_list[m]
          assign(simul_names[l],simul_data)
        }
      }
    }
    
    # final simulated dataset
    simul_data = get(simul_names[1])
    for(l in 2:n_dif){
      simul_data <- rbind(simul_data,get(simul_names[l]))
    }
    simul_data$score <- apply(simul_data,1,sum,na.rm=TRUE)
    if(sc_gp >1){
      simul_data$score <- cut(simul_data$score,breaks=quantiles,labels=FALSE)
    }
    R = as.integer(levels(factor(simul_data$score)))
    tab_b = data.frame(level_R=R,S = rep(0,length(R)))
    for(i in 1:length(R)){
      data_r = subset(simul_data,simul_data$score == as.numeric(R[i]) & is.na(simul_data[,item])==FALSE,select=item)
      tab_b[i,][2] <- sum(data_r[,item],na.rm=TRUE)/dim(data_r)[1]
    }
    lines(tab_b$level_R,tab_b$S,lty = 3,col="red")
  }  
  lines(tab$level_R,tab$S,col="blue",lwd=2)
  legend("bottomright",legend=c("real success rate","sim success rate"),col=c("blue","red"),lty=c(1,3),cex=0.5)
}

simulation_5 = function(data,item,item_dif,itemlist,constraint,B,sc_gp=1,diffvar,unq = TRUE){
  # first differenciation
  result = differenciate(data,item_dif,diffvar)
  data = result$data
  dif_list=  result$diff_names
  itemlist = c(itemlist[-which(itemlist==item_dif)])
  
  # second differenciation
  result = differenciate(data,item,diffvar)
  data=result$data
  to_dif_list = result$diff_names
  len = length(to_dif_list)
  
  # levels of the covariate
  level = levels(as.factor(data[,diffvar]))
  
  if(unq){
    par(mfrow=select_par(len+1))}
  # simulation for the original item
  simulation_2_int(data=data,item=item,itemlist=itemlist,diff_names=dif_list,constraint=constraint, B=B,sc_gp=sc_gp)
  # simulation for the differenciated items
  itemlist = itemlist[-which(itemlist == item)]
  for(i in 1:len){
    simulation_5_int(data=data,item=to_dif_list[i],itemlist=itemlist,dif_list=dif_list,to_dif_list=to_dif_list,level=level,constraint=constraint,B=B,sc_gp=sc_gp)
  }
}




#===================================================================================#




#===================================================================================#
simulation = function(data,items,itemlist,constraint,B,sc_gp=1,diffvar=NULL,item_dif=NULL,samePlot=FALSE){
  if(is.null(diffvar)==FALSE){
    if(is.null(item_dif)==FALSE){
      simulation_5(data,item,item_dif,itemlist,constraint,B,sc_gp,diffvar)
    }else{
    if(samePlot==TRUE){
      simulation_4(data,items,itemlist,constraint,B,sc_gp,diffvar)
    }else{
      simulation_3(data,items,itemlist,constraint,B,sc_gp,diffvar)}}
  }else{
    simulation_1(data,items,itemlist,constraint,B,sc_gp)  
  }
}
#===================================================================================#


#===================================================================================#

test_item_int = function(dataList,item,itemlist,covariates,constraint,eps){
  n = length(dataList)
  names = names(dataList)
  covar = as.character(covariates[,names[[1]]][1])
  study = names[[1]]
  result_item = data.frame("covar"=c(covar),"study" = c(study),"LRT"=c(0),"chisq"=c(0),"pvalue"=c(0))
  for(i in 1:n){
    study = names[[i]]
    covar = levels(covariates[,study])
    m= length(covar)
    for(j in 1:m){
      diffvar = covar[j]
      result = tryCatch(dif_poly(dataList[[study]],item=item,itemlist=itemlist,diffvar=diffvar,constraint=constraint,toPlot=FALSE),error=function(e) NA)
      if(is.na(result)){
        data_int = data.frame("covar"=c(diffvar),"study"=c(study),"LRT"="ERROR","chisq"="ERROR","pvalue"="ERROR")    
      }else{
        result = result$LRT
        data_int = data.frame("covar"=c(diffvar),"study"=c(study),"LRT"=result[,"LRT"],"chisq"=result[,"chisq"],"pvalue"=result[,"pvalue"])  
      }
      result_item = rbind(result_item,data_int)
    }}
  result_item = result_item[-1,]
  rownames(result_item) <- NULL
  result_item$test <- (result_item$pvalue < eps)  
  result_item$test[result_item$pvalue=="ERROR"] <- "ERROR"
  return(result_item)
}


test_item = function(dataList,items,itemlist,covariates,constraint,eps){
  result_Tot = list()
  len = length(items)
  for(i in 1:len){
    result = test_item_int(dataList,items[i],itemlist,covariates,constraint,eps)
    result_Tot = list.append(result_Tot,result)
  }
  names(result_Tot) <- items
  return(result_Tot)
}
#===================================================================================#


#===================================================================================#
select_aleat = function(data,refvar,n){
  values = levels(as.factor(data[,refvar]))  
  p=length(values)
  var_names = paste(refvar,1:p,sep="_")
  for(i in 1:p){
    data_t = subset(data,data[,refvar]==values[i])
    data_t = data_t[sample(nrow(data_t),size=n,replace=FALSE),]
    rownames(data_t) <- NULL
    assign(var_names[i],data_t)
  }
  data_t = get(var_names[1])
  for(i in 2:p){
    data_t=rbind(data_t,get(var_names[i]))
  }
  return(data_t)  
}
#===================================================================================#

#===================================================================================#
select_items = function(data,expr){
  l = grep(paste(expr,"",sep=""), colnames(data))
  return(l)
}

select_not_na = function(data){
  data <- data[,colSums(is.na(data)) != nrow(data)]
  return(data)
}

select_par = function(n){
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


#===================================================================================#
# importation of the dataset
#setwd("your_path")
#data <- read.csv("heiq.csv")

# creation and transformation of the databases
study_1 <- subset(data,data$study==1)
study_2 <- subset(data,data$study==2)
study_3 <- subset(data,data$study==3)
study_4 <- subset(data,data$study==4)
study_5 <- subset(data,data$study==5)


study_1 = select_not_na(study_1)
study_2 = select_not_na(study_2)
study_3 = select_not_na(study_3)
study_4 = select_not_na(study_4)
study_5 = select_not_na(study_5)

datalist = list("study_1" = study_1,"study_2" =study_2,"study_3" =study_3,"study_4" =study_4,"study_5" =study_5)
covariables = data.frame("study_1"=c("Sex","over64"),"study_2"=c("Sex","over48"),"study_3"=c("Sex","over72"),"study_4"=c("Sex","over76"),"study_5"=c("Sex","over51"))

# sub_scales
healthBehaviour =paste("Heiq",c(1,9,13,19),sep="")
engagementLife = paste("Heiq",c(2,5,8,10,15),sep="")
emotionWellbeing = paste("Heiq",c(4,7,12,14,18,21),sep="")
selfMonitoring = paste("Heiq",c(3,6,11,16,17,20),sep="")
beingConstructive = paste("Heiq",c(27,34,36,39,40),sep="")
skillAcquisition = paste("Heiq",c(23,25,26,30),sep="")
socialIntegration = paste("Heiq",c(22,28,31,35,37),sep="")
healthNavigation = paste("Heiq",c(24,29,32,33,38),sep="")





