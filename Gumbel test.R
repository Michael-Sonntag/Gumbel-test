# Gumbel test

# Test with constant volatility

ctest<-function(dInc,p,quantil){
m<-length(which(dInc!=0))
a<- sqrt(2*log(2*m));b<-a-((log(log(2*m))+log(4*pi))/(2*a))

K<-mean(abs(rnorm(1000000))^p)
x1<-mean(abs(dInc)^p)
d1<-(x1/K)^(-1/p)

tStatistik<-a*(max(abs(dInc)*d1)-b)
return(tStatistik>=-log(-log(quantil)))

}

# Test with time varying volatility

test<-function(dInc,p,alpha,quantil){    
m<-length(which(dInc!=0))
n<-length(dInc)
a<- sqrt(2*log(2*m));b<-a-((log(log(2*m))+log(4*pi))/(2*a))
                                                                  
hn<-floor(n^(2*alpha/(2*alpha+1)))         
temp<-floor(n/hn)                                                                             
l<- hn*temp

K<-mean(abs(rnorm(1000000))^p)

A1<-matrix(abs(dInc[1:l]),nrow=temp,ncol=hn,byrow=T)
x1<-rowSums(abs(A1)^p)
d1<-(x1/(hn*K))^(-1/p)
tStatistik<-a*(max(A1*d1)-b)
return(tStatistik>=-log(-log(quantil)))

 }

# Sequential test and jump positions

seqTestPos<-function(dInc,p,alpha,quantil){

K<-mean(abs(rnorm(2000000))^p)


stop<-test(dInc,p,alpha,quantil)
while(stop){


m<-length(which(dInc!=0))
n<-length(dInc)
a<- sqrt(2*log(2*m))                     
b<-a-((log(log(2*m))+log(4*pi))/(2*a))

hn<-floor(n^(2*alpha/(2*alpha+1)))         
temp<-floor(n/hn)                                                                             
l<- hn*temp

A1<-matrix(abs(dInc[1:l]),nrow=temp,ncol=hn,byrow=T)
x1<-rowSums(abs(A1)^p)
d1<-(x1/(hn*K))^(-1/p)
obs<-(a*A1*d1)-a*b
temp1<-c(t(obs))


dInc[which(temp1==max(temp1))]<-0
stop<-test(dInc,p,alpha,quantil)


         }
return(dInc)
   }








