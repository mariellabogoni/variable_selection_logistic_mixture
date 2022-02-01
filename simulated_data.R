#########
# Generating the data
#########

beta_true<-cbind(c(1,-1,0,1,0), c(-1,0,1,0,1),c(-0.5,0,-0.5,0,-0.5)) 
k<-ncol(beta_true)
w_true<-rep(1/k,k)
n=200
N<-50
param<-nrow(beta_true)
X<-rep(1,n)
for (t in 1:(param-1)) X<-cbind(X,rnorm(n,0,1))

S_true<-NULL
y<-NULL
p<-NULL
for (i in 1:n) {
  S_true[i]<-rDiscreta(w_true)
  p[i]<-logit_inverse(X[i,]%*%beta_true[,S_true[i]])
  y[i]<-rbinom(1,N,p[i])
}
remove(p)