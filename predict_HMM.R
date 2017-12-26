library("HMM")
library ("lattice")

#data for all host sequences
data<-read.csv("HMM_sorted.csv",header=FALSE)

#hidden states and observed states (symbols)
hstates <- c("human","viral")
ostates <- c("human","viral")

#parse prior parameters for each fasta entry
get_prior <- function(x)
{
datarow <- x
datarow <- datarow[datarow!=""]

#host seq ID
seqID <- datarow[1]

#transmission probability
transp <- matrix(nrow=2,byrow=TRUE,c(datarow[2:5]))
class(transp) <- "numeric"

#emission probability
emitp <- matrix(nrow=2,byrow=TRUE,c(datarow[6:9]))
class(emitp) <- "numeric"

#symbol sequence
obs_state <- datarow[11:length(datarow)]

#intialize hmm
hmm <- initHMM(hstates,ostates,NULL,transp,emitp)

#compute most probable
obs_state.pos <- viterbi(hmm, obs_state)

#ouput to file
line<- c(seqID,obs_state.pos)
write(line, file="hmmout.txt",ncolumns=length(line),append=TRUE)
}

#apply to matrix
apply(data, 1, get_prior)
