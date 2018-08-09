#===================================================================================#
#                                                                                   #
#          INTERNSHIP - MASTER I/II - BIOSTATISTICS - Work on the true database     #
#                       Latent variables model - IRT - Rasch                        #
#                                      S. SIGL                                      #
#                                                                                   #
#===================================================================================#

#%%%%%%%%%%%%%%
# we import the functions we created in another file
source("your_path/github_code.R")

#%%%%%%%%%%%%%%
# we apply the functions

difPolyTot(data=data,items=selfMonitoring,itemlist=selfMonitoring,
           diffvar="study",constraint="gpcm",toPlot=TRUE)

simulation(data=data,items=healthBehaviour,itemlist=healthBehaviour,
           constraint="gpcm",B=10,scoreGrp=5)

simulation(data=data,items="Heiq1",itemlist=healthBehaviour,
                      constraint="gpcm",B=10,scoreGrp=1,diffvar="Sex")

simulation(data=data,items="Heiq1",itemlist=healthBehaviour,constraint="gpcm",B=10,scoreGrp=3,diffvar="Sex",samePlot=TRUE)


testItem_Tot(dataList=datalist,item="healthBehaviour",itemlist=healthBehaviour,
             covariates=covariables,constraint="gpcm",eps=0.05)


