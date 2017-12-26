t.test2 <- function(x,m0=0,equal.variance=FALSE)
{
    n1 <- x[2]
    class(n1) <- "numeric"
    n2 <- x[3]
    class(n2) <- "numeric"
    m1 <- x[4]
    class(m1) <- "numeric"
    m2 <- x[5]
    class(m2) <- "numeric"
    s1 <- x[6]
    class(s1) <- "numeric"
    s2 <- x[7]
    class(s2) <- "numeric"
    
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    p=2*pt(-abs(t),df) 
    line <- c(x[1],p,x[8])
    write(line, file="viral_t-test.out",ncolumns=length(line),append=TRUE)
}

#read
data <- read.csv("bed_cov.out")
apply(data,1,t.test2)