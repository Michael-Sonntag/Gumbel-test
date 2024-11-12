# Power for fixed sample size

# n: sample size
# H: Hurst exponent
# m: Monte-Carlo simulations
# J: Jump size
# pos1: Jump position
# k: Number of increases of jump size

h<-function(n,H,m,J,pos1,k){

res1<-matrix(nrow=m,ncol=k)

H2 <- 2*H
matcov <- matrix(0, n, n)
for(i in (1:n)) {
r <- 0.5*(abs(i)^H2 + abs(1:n)^H2 - abs(i-(1:n))^H2)
r <- r/(n^H2)
matcov[i, 1:n] <- r
matcov[1:n,i] <- matcov[i,1:n]
}
L <- chol(matcov)


p<-0.9;quantil<-0.95
m1<-n-1
a<- sqrt(2*log(2*m1));b<-a-((log(log(2*m1))+log(4*pi))/(2*a))
K<-mean(abs(rnorm(1000000))^p)



for(j in 1:m){

Z <- rnorm(n)
fBm <- t(L) %*% Z
fBm <- c(0, fBm)


for(i in 1:k){

y1<-numeric(n+1)
y1[(pos1:(n+1))]<-J*(i-1)
fBmJumps<-fBm+y1




dInc<-diff(fBmJumps,lag=1,differences=2)

x1<-mean(abs(dInc)^p)
d1<-(x1/K)^(-1/p)
tStatistik<-a*(max(abs(dInc)*d1)-b)
res1[j,i]<-(tStatistik>=-log(-log(quantil)))
      }

    }
return(res1)

 }

n<-500;H<-0.3


# Example

a1<-h(500,0.3,10000,0.01,200,150)


plot(seq(0,1.49,0.01),colMeans(a1),ylab="Empirical Power",xlab="Jump Size")
LA<-n^(-(H-0.01))*sqrt(2*log(n))
abline(v=LA, col="blue")
a<-sqrt(2*log(n));b<-a-((log(log(2*n))+log(4*pi))/(2*a))
ROT<- -log(-log(0.95))/(n^(H)*a)+(2*b)/(n^(H))
abline(v=ROT,col="red")
legend("topleft",legend=c("Rule of Thumb","Local Alternative"),col=c("red","blue"),lty=1)

