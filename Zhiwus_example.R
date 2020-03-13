#slide 4
#QTNs on CHR 1-5, leave 6-10 empty
myGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
source("http://zzlab.net/StaGen/2020/R/G2P.R")
source("http://zzlab.net/StaGen/2020/R/GWASbyCor.R")
X=myGD[,-1]
index1to5=myGM[,2]<6
X1to5 = X[,index1to5]
set.seed(99164)
mySim=G2P(X= X1to5,h2=.75,alpha=1,NQTN=10,distribution="norm")
p= GWASbyCor(X=X,y=mySim$y)

#slide 5
#False positives
color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot(t(-log10(p))~seq(1:m),col=color.vector[myGM[,2]])
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")

#slide 6
#QQ plot
p.obs=p[!index1to5]
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),-log10(p.obs[order.obs]))
abline(a = 0, b = 1, col = "red")

#slide 7
#Phenotypes by genotypes
order.obs=order(p.obs)
X6to10=X[,!index1to5]
Xtop=X6to10[,order.obs[1]]

boxplot(mySim$y~Xtop)

#slide 8
#correlation
PCA=prcomp(X)
plot(mySim$y,PCA$x[,2])
cor(mySim$y,PCA$x[,2])

#slide 9
#linear regression
set.seed(99164)
s=sample(length(mySim$y),10)
plot(mySim$y[s],PCA$x[s,2])
cor(mySim$y[s],PCA$x[s,2])

#slide 11
#example from ten individuals
cbind(mySim$y[s],1, PCA$x[s,2],Xtop[s])

#slide 15
#action in R
y=mySim$y
X=cbind(1, PCA$x[,2],Xtop)
LHS=t(X)%*%X
C=solve(LHS)
RHS=t(X)%*%y
b=C%*%RHS
yb=X%*%b
e=y-yb
n=length(y)
ve=sum(e^2)/(n-1)
vt=C*ve
t=b/sqrt(diag(vt))
p=2*(1-pt(abs(t),n-2))

#slide 16
#phenotypes by genotypes
LM=cbind(b, t, sqrt(diag(vt)), p)
rownames(LM)=cbind("Mean", "PC2","Xtop")
colnames(LM)=cbind("b", "t", "SD","p")
LM

#slide 17
#loop through genome
G=myGD[,-1]
n=nrow(G)
m=ncol(G)
P=matrix(NA,1,m)
for (i in 1:m){
  x=G[,i]
  if(max(x)==min(x)){
    p=1}else{
      X=cbind(1, PCA$x[,2],x)
      LHS=t(X)%*%X
      C=solve(LHS)
      RHS=t(X)%*%y
      b=C%*%RHS
      yb=X%*%b
      e=y-yb
      n=length(y)
      ve=sum(e^2)/(n-1)
      vt=C*ve
      t=b/sqrt(diag(vt))
      p=2*(1-pt(abs(t),n-2))
    } #end of testing variation
  P[i]=p[length(p)]
} #end of looping for markers

#slide 18
#QQ plot
p.obs=P[!index1to5]
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),
     -log10(p.obs[order.obs]), ylim=c(0,7))
abline(a = 0, b = 1, col = "red")

#slide 19
#using three PCs
G=myGD[,-1]
n=nrow(G)
m=ncol(G)
P=matrix(NA,1,m)
for (i in 1:m){
  x=G[,i]
  if(max(x)==min(x)){
    p=1}else{
      X=cbind(1, PCA$x[,1:3],x)
      LHS=t(X)%*%X
      C=solve(LHS)
      RHS=t(X)%*%y
      b=C%*%RHS
      yb=X%*%b
      e=y-yb
      n=length(y)
      ve=sum(e^2)/(n-1)
      vt=C*ve
      t=b/sqrt(diag(vt))
      p=2*(1-pt(abs(t),n-2))
    } #end of testing variation
  P[i]=p[length(p)]
} #end of looping for markers

#slide 20
#QQ plot
p.obs=P[!index1to5]
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),
     -log10(p.obs[order.obs]), ylim=c(0,7))
abline(a = 0, b = 1, col = "red")

#slide 21
#QQ plot
p.obs=P
m2=length(p.obs)
p.uni=runif(m2,0,1)
order.obs=order(p.obs)
order.uni=order(p.uni)

plot(-log10(p.uni[order.uni]),
     -log10(p.obs[order.obs]), )
abline(a = 0, b = 1, col = "red")

#slide 22
color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
plot(t(-log10(P))~seq(1:m),col=color.vector[myGM[,2]])
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")
