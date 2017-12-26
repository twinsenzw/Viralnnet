library(nnet)

#load data and nnet model
load(file="NNmodel_350k.rda")
data<-read.csv("kmer_prep.csv",header=FALSE)

#Predicting on testset
pred <- predict(dataNN, data[,-1], type="class")

#write predicted output
write(pred,file="predoutput.txt")
