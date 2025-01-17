# R-Code to reproduce Figure 4 

library(ggplot2)

# n: number of observations
# H: Hurst exponent
# m: number of Monte-Carlo simulations
# J: Jump size
# pos1: position of jump

h1<-function(n,H,m,J,pos1){

res<-matrix(nrow=m,ncol=6)

# Covariance matrix fBm

H2 <- 2*H
matcov <- matrix(0, n, n)
for(i in (1:n)) {
r <- 0.5*(abs(i)^H2 + abs(1:n)^H2 - abs(i-(1:n))^H2)
r <- r/(n^H2)
matcov[i, 1:n] <- r
matcov[1:n,i] <- matcov[i,1:n]
}
L <- chol(matcov)


for(j in 1:m){


Z <- rnorm(n)
fBm <- t(L) %*% Z
fBm <- c(0, fBm)


for(i in 1:6){

y1<-numeric(n+1)
y1[(pos1:(n+1))]<-J*(i-1)
fBmJumps<-fBm+y1



C<-1/(2*log(2))
num<-sum((fBmJumps[3:n]-fBmJumps[1:(n-2)])^2)
dem<-sum((fBmJumps[2:n]-fBmJumps[1:(n-1)])^2)
res[j,i]<-(C*log(num/dem))

       }

    }
return(res)

 }

# Example

a1<-h1(2000,0.2,5000,0.5,200)
a2<-h1(2000,0.4,5000,0.5,200)
a3<-h1(2000,0.6,5000,0.5,200)


# Display of graphic

EstimatedValues<-c(a1[,1],a1[,2],a1[,3],a1[,4],a1[,5],a1[,6],a2[,1],a2[,2],a2[,3],a2[,4],a2[,5],a2[,6],a3[,1],a3[,2],a3[,3],a3[,4],a3[,5],a3[,6])
SimulatedData<-rep(c("H=0.2","H=0.4","H=0.6"),each=30000);JumpSize<-rep(c("0","0.5","1","1.5","2","2.5"),each=5000) 
sim<-data.frame(SimulatedData,JumpSize,EstimatedValues)
ggplot(sim,aes(x=SimulatedData,y=EstimatedValues,fill=JumpSize))+
geom_boxplot()+
xlab("")

# Randomly generated jump time

# n: number of observations
# H: Hurst exponent
# m: number of Monte-Carlo simulations
# J: Jump size


h2<-function(n,H,m,J){

res<-matrix(nrow=m,ncol=6)

# Covariance matrix fBm

H2 <- 2*H
matcov <- matrix(0, n, n)
for(i in (1:n)) {
r <- 0.5*(abs(i)^H2 + abs(1:n)^H2 - abs(i-(1:n))^H2)
r <- r/(n^H2)
matcov[i, 1:n] <- r
matcov[1:n,i] <- matcov[i,1:n]
}
L <- chol(matcov)

pos<-sample(x=n-1,size=m, replace=T)


for(j in 1:m){


Z <- rnorm(n)
fBm <- t(L) %*% Z
fBm <- c(0, fBm)


for(i in 1:6){

y1<-numeric(n+1)
y1[((pos[j]+1):(n+1))]<-J*(i-1)
fBmJumps<-fBm+y1



C<-1/(2*log(2))
num<-sum((fBmJumps[3:n]-fBmJumps[1:(n-2)])^2)
dem<-sum((fBmJumps[2:n]-fBmJumps[1:(n-1)])^2)
res[j,i]<-(C*log(num/dem))

       }

    }
return(res)

 }

# Example

a1<-h2(2000,0.2,2000,0.5)
a2<-h2(2000,0.4,2000,0.5)
a3<-h2(2000,0.6,2000,0.5)



# Display of graphic

EstimatedValues<-c(a1[,1],a1[,2],a1[,3],a1[,4],a1[,5],a1[,6],a2[,1],a2[,2],a2[,3],a2[,4],a2[,5],a2[,6],a3[,1],a3[,2],a3[,3],a3[,4],a3[,5],a3[,6])
SimulatedData<-rep(c("H=0.2","H=0.4","H=0.6"),each=12000);JumpSize<-rep(c("0","0.5","1","1.5","2","2.5"),each=2000) 
sim<-data.frame(SimulatedData,JumpSize,EstimatedValues)
ggplot(sim,aes(x=SimulatedData,y=EstimatedValues,fill=JumpSize))+
geom_boxplot()+
xlab("")


