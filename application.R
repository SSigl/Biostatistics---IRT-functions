#===================================================================================#
#                                                                                   #
#          INTERNSHIP - MASTER I/II - BIOSTATISTICS - Work on the true database    #
#                       Latent variables model - IRT - Rasch                        #
#                                   S. SIGALLA                                      #
#                                                                                   #
#===================================================================================#

#%%%%%%%%%%%%%%
# we import the functions we created in another file
source("your_path/github_code.R")

#%%%%%%%%%%%%%%
# importation of the dataset
setwd("your_path")
data <- read.csv("heiq.csv")

# creation and transformation of the databases
study_1 <- subset(data,data$study==1)
study_2 <- subset(data,data$study==2)
study_3 <- subset(data,data$study==3)
study_4 <- subset(data,data$study==4)
study_5 <- subset(data,data$study==5)


study_1 = selectNotNa(study_1)
study_2 = selectNotNa(study_2)
study_3 = selectNotNa(study_3)
study_4 = selectNotNa(study_4)
study_5 = selectNotNa(study_5)

dataList = list("study_1" = study_1,"study_2" =study_2,"study_3" =study_3,"study_4" =study_4,"study_5" =study_5)
covariates = data.frame("study_1"=c("Sex","over64"),"study_2"=c("Sex","over48"),"study_3"=c("Sex","over72"),"study_4"=c("Sex","over76"),"study_5"=c("Sex","over51"))

# sub_scales
healthBehaviour =paste("Heiq",c(1,9,13,19),sep="")
engagementLife = paste("Heiq",c(2,5,8,10,15),sep="")
emotionWellbeing = paste("Heiq",c(4,7,12,14,18,21),sep="")
selfMonitoring = paste("Heiq",c(3,6,11,16,17,20),sep="")
beingConstructive = paste("Heiq",c(27,34,36,39,40),sep="")
skillAcquisition = paste("Heiq",c(23,25,26,30),sep="")
socialIntegration = paste("Heiq",c(22,28,31,35,37),sep="")
healthNavigation = paste("Heiq",c(24,29,32,33,38),sep="")


#%%%%%%%%%%%%%%
# we apply the functions

difPolyTot(data=data,items=selfMonitoring,itemlist=selfMonitoring,
           diffvar="study",constraint="gpcm",toPlot=TRUE)

simulation_1_Tot(data=data,items=healthBehaviour,itemlist=healthBehaviour,
                 constraint="gpcm",B=30,scoreGrp=3)

simulation_1_diff_Tot(data=data,items=healthBehaviour,itemlist=healthBehaviour,
                      constraint="gpcm",B=10,scoreGrp=3,diffvar="Sex")



