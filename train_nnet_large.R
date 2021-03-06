library(nnet)
data <- read.csv("all_kmer_prep.csv",header=FALSE)
Tlabel <- class.ind(data$V1)
datatrain <- sample(1:4579531,3200000)
datatest <- setdiff(1:4579531,datatrain)
dataNN=nnet(data[datatrain,-1],Tlabel[datatrain,],size=60,softmax=TRUE,maxit=2000,decay=0.01,MaxNWts=15000)
table(predict(dataNN, data[datatest,-1], type="class"),data[datatest,]$V1)
save(dataNN,file="dataNN_large.rda")
