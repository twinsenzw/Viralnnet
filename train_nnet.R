library(nnet)

#load data
data <- read.csv("350k_train.csv",header=FALSE)
Tlabel <- class.ind(data$V1)

#optional: prepare training set from "seen species"
datatrain <- sample(1:351436,250000)
datatest <- setdiff(1:250000,datatrain)

#Fitting nnet
dataNN=nnet(data[datatrain,-1],Tlabel[datatrain,],size=25,softmax=TRUE,maxit=1000,decay=0.01,MaxNWts=15000)

#Predicting on testset
table(predict(dataNN, data[datatest,-1], type="class"),data[datatest,]$V1)

#save nnet object
save(dataNN,file="NNmodel_350k.rda")