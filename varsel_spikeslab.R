library(Matrix)
library(gtools)
library(mvtnorm)
library(compiler)
library(BayesLogit)
library(MASS)
library(label.switching)
library(coda)
library(bayestestR)
#
variable_selection<-function(n_k,beta_k,y_k,N_k,X_k,ind_k,prior_mean_k,cov_matrix_k){
  prob<-NULL
  prob_ind1<-NULL
  prob_ind0<-NULL
  p<-length(prior_mean_k)
  prob[1]<-1
  ind_k[1]<-1
  #  
  # prior inclusion probability
  p_prior<-1/2     
  #polya-gamma
  w_polya<-rpg.devroye(n_k,N_k,matrix(X_k[,which(ind_k==1)],nrow=n_k)%*%matrix(beta_k[which(ind_k==1)],nrow=sum(ind_k)))
  z<-(y_k-N_k/2)/w_polya
  
    
  for (t in 2:p){
    ind_k[t]<-1
    matrix_aux<-solve(cov_matrix_k[which(ind_k==1),which(ind_k==1)])
    post_var<-solve(matrix_aux+t(matrix(X_k[,which(ind_k==1)],nrow=n_k))%*%diag(w_polya, nrow=n_k)%*%matrix(X_k[,which(ind_k==1)],nrow=n_k))
    A<-(matrix_aux%*%prior_mean_k[which(ind_k==1)] + t(matrix(X_k[,which(ind_k==1)],nrow=n_k))%*%diag(w_polya, nrow=n_k)%*%matrix(z,ncol=1))
    post_mean<-post_var%*%A
    
    log_det_postvar<-log(det(post_var))
    prob_ind1[t]<-log(p_prior)+(0.5)*(log_det_postvar-sum(ind_k)*log(cov_matrix_k[1,1]))-(0.5)*(t(prior_mean_k[which(ind_k==1)])%*%matrix_aux%*%prior_mean_k[which(ind_k==1)]-t(post_mean)%*%A) 
    #    
    ind_k[t]<-0
    matrix_aux<-solve(cov_matrix_k[which(ind_k==1),which(ind_k==1)])
    post_var<-solve(matrix_aux+t(matrix(X_k[,which(ind_k==1)],nrow=n_k))%*%diag(w_polya, nrow=n_k)%*%matrix(X_k[,which(ind_k==1)],nrow=n_k))
    A<-(matrix_aux%*%prior_mean_k[which(ind_k==1)]+t(matrix(X_k[,which(ind_k==1)],nrow=n_k))%*%diag(w_polya, nrow=n_k)%*%matrix(z,ncol=1))
    post_mean<-post_var%*%A
    
    log_det_postvar<-log(det(post_var))
    prob_ind0[t]<-log(1-p_prior)+(0.5)*(log_det_postvar-sum(ind_k)*log(cov_matrix_k[1,1]))-(0.5)*(t(prior_mean_k[which(ind_k==1)])%*%matrix_aux%*%prior_mean_k[which(ind_k==1)]-t(post_mean)%*%A)
       
    #posterior inclusion probability
    prob[t]<-exp(-log(1 + exp(prob_ind0[t]-prob_ind1[t])))  
    ind_k[t]<-rbinom(1,1,prob = prob[t])
  }
  #  
  matrix_aux<-solve(cov_matrix_k[which(ind_k==1),which(ind_k==1)])
  post_var<-solve(matrix_aux+t(matrix(X_k[,which(ind_k==1)],nrow=n_k))%*%diag(w_polya, nrow=n_k)%*%matrix(X_k[,which(ind_k==1)],nrow=n_k))
  A<-(matrix_aux%*%prior_mean_k[which(ind_k==1)]+t(matrix(X_k[,which(ind_k==1)],nrow=n_k))%*%diag(w_polya, nrow=n_k)%*%matrix(z,ncol=1))
  post_mean<-post_var%*%A
  
  # updating the betas with ind==1 and setting the others as zero
  beta_k<-matrix(0,nrow=nrow(beta_k),ncol=1)
  beta_k[which(ind_k==1)]<-rmvnorm(1,mean=post_mean,sigma=post_var, method=c("chol"), pre0.9_9994 = TRUE)
  beta_k[which(ind_k==0)]<-0
  return(list(ind_k,beta_k))
}
#
update_slatent<-function(w_iter,betas,y,X,N,ind){
  s<-NULL
  aux<-NULL
  for(j in 1:length(y)){ 
    for (i in 1:length(w_iter)){
      aux[i]<-log(w_iter[i])+dbinom(y[j],N,logit_inverse(X[j,ind[,i]==1]%*%betas[ind[,i]==1,i]),log=TRUE)
    }
    aux<-exp(aux-max(aux))
    s[j]<-rDiscreta(aux/sum(aux))
    if (is.na(s[j])==TRUE) s[j]<-rDiscreta(rep(1/length(w_iter),length(w_iter)))
  }
  return(s)
}
#
rDiscreta<-function(p){
  u<-runif(1)
  P<-cumsum(p)
  val<-sum(P<u)+1
  val
}
#
logit_inverse <- function(x) {
  1/(1+exp(-x))
}
#
steps<-55000
burn<-5000
saltos<-10
caminho<-"~/Documents/"

#########
# Generating the data
#########

source("simulated_data.R", local = FALSE)

#########
# Starting the method
#########

tpr<-NULL;fpr<-NULL
S_final<-NULL
beta_est<-NULL
w_est<-NULL
indicator_est<-NULL
w<-NULL
S<-NULL
g_size<-NULL
beta<-matrix(0,nrow=param,ncol=k)
indicator<-matrix(0,nrow=param,ncol=k)
prior_covmatrix<-list()
#hyperparameters of beta and w prior
alfa<- c(rep(1,k))
prior_mean<-0
prior_var<-10

# Initializing the components
w<-rdirichlet(1, alpha = alfa)
S<-sample(1:k, n, replace = TRUE, prob = w[1,])

for ( f in 1:k) {
  g_size[f]<-length(y[S==f])
}

#Dealing with empty components
if(length(which(g_size==0))!=0){  
  c<-which(g_size==0) #empty component indexes
  for (t in 1:length(c)) {
    comp<-which.max(g_size)                                    #choosing the biggest component to take an observation from
    pos<-sample(x=c(which(S==comp)),size=1, replace = FALSE)   #randomly choosing the observation in the component 
    S[pos]<-c[t]
    cat("\n", c(t,c[t],S[pos]))
  }
  for ( f in 1:k){ g_size[f]<-length(y[S==f]) }
}

#Inicializing Betas through the prior distribution
for (f in 1:k) {
  prior_covmatrix[[f]]<-diag(rep(prior_var,param))
  beta[,f]<-rmvnorm(1, mean =rep(prior_mean,param), sigma = prior_covmatrix[[f]], method=c("chol"), pre0.9_9994 = TRUE)
  indicator[,f]<-rep(0,param)
  
}


enableJIT(3)
for (i in 2:steps)   #main loop
{
  
  for (h in 1:k) {
    out<-variable_selection(g_size[h],matrix(beta[,h],ncol=1),y[S==h],N,
                            matrix(X[S==h,],nrow=g_size[h]),indicator[,h],rep(prior_mean,param),
                            prior_covmatrix[[h]])
    indicator[,h]<-out[[1]]
    beta[,h]<-out[[2]]
    
  }
  
  # Updating w and S
  w<-rdirichlet(1,g_size+alfa)   
  if(k==1){S<-rep(1,n)
  }else{ S<-update_slatent(w,beta,y,X,N,indicator)}
  
  if(i>burn & i%%saltos==0)
  {        
    cat('\n',beta,file=paste(caminho,"Betas",".txt",sep=""),append=T)
    cat('\n',indicator,file=paste(caminho,"Indicators",".txt",sep=""),append=T)
    cat('\n',w,file=paste(caminho,"w",".txt",sep=""),append=T)
    cat('\n',S,file=paste(caminho,"S",".txt",sep=""),append=T)
    logvero<-0
    for (j in 1:length(y)){
      logvero<- logvero + log(w[S[j]])+ dbinom(y[j],N,logit_inverse(X[j,indicator[,S[j]]==1]%*%beta[indicator[,S[j]]==1,S[j]]),log=TRUE)
    }
    cat('',logvero,file=paste(caminho,"log_like",".txt",sep=""),append=T)
  }
  
  for ( g in 1:k) {
    g_size[g]<-length(y[S==g])
  }
  
  #Dealing with empty components
  if(length(which(g_size==0))!=0){  
    c<-which(g_size==0) #empty component indexes
    for (t in 1:length(c)) {
      comp<-which.max(g_size)                                    #choosing the biggest component to take an observation from
      pos<-sample(x=c(which(S==comp)),size=1, replace = FALSE)   #randomly choosing the observation in the component 
      S[pos]<-c[t]
    }
    for ( f in 1:k){ g_size[f]<-length(y[S==f]) }
  }
  
}

######
#  Analyzing the results
######

size<-((steps-burn)/saltos) 
betas<-scan(file=paste(caminho,"Betas",".txt",sep=""))
betas<-matrix(betas,nrow=size,ncol=param*k,byrow=TRUE)
indicator<-scan(file=paste(caminho,"Indicators",".txt",sep=""))
indicator<-matrix(indicator,nrow=size,ncol=param*k,byrow=TRUE)
loglike<-scan(file=paste(caminho,"log_like",".txt",sep=""))
w<-scan(file=paste(caminho,"w",".txt",sep=""))
w<-matrix(w,nrow=size,ncol=k,byrow=TRUE)
S<-scan(file=paste(caminho,"S",".txt",sep=""))
S<-matrix(S, nrow = size, ncol = length(y), byrow = T)

############################
#  label switching         
###########################  
if(k!=1){
  maxindex<-which.max(logvero)
  J<-2*param + 1
  mcmc<-array(data = NA, dim = c(size,k,J))
  j=1
  for (a in 1:k) {
    mcmc[,a,]<-cbind(betas[,j:(a*param)],indicator[,j:(a*param)],w[,a])
    j<-j+param
    
  }
  L<-label.switching(method = "ECR", zpivot = S_true, z = S, K = k)
  
  mcmc_final<-permute.mcmc(mcmc,L$permutations$ECR)   #new mcmc sample
  w<-mcmc_final$output[,,J]
  
  #Uptading Betas
  j=1
  for (b in 1:k) { 
    
    betas[,j:(b*param)]<-mcmc_final$output[,b,1:param]  
    j<-j+param
  }
  j=1
  for (b in 1:k) { 
    indicator[,j:(b*param)]<-mcmc_final$output[,b,(param+1):(J-1)]  
    j<-j+param
  }
  #Uptading S
  for (i in 1:size){
    for(j in 1:length(y))
    {    aux<-S[i,j]
    S[i,j]<-L$permutations$ECR[i, aux]}
  }
}

# checking convergence 
log_v<-mcmc(loglike)  
gewekee<-geweke.diag(log_v)[[1]]   

# Estimates      
w_est<-apply(w, 2, mean)
beta_est<-apply(betas,2,mean)
beta_est<-matrix(beta_est, nrow = param)
indicator_est<-apply(indicator, 2, mean)
indicator_est<-round(matrix(indicator_est, nrow = param),2)

# Classifing the observations
for (i in 1:length(y)){
  prob<-NULL
  for (j in 1:k) prob[j]<-sum(S[,i]==j)/nrow(S)
  S_final[i]<-rDiscreta(prob)
}

#95% highest probability credible interval
for (h in 1:k) {
    ciw<-ci(w[,h], method = "HDI")
    cat('\n',c(ciw$CI_low,ciw$CI_high ),file=paste(caminho,"CI_weights",sep=""),append=T)
}
for (h in 1:(param*k)) {
    cib<-ci(betas[,h], method = "HDI")
    cat('\n', c(cib$CI_low,cib$CI_high ),file=paste(caminho,"CI_Betas",sep=""),append=T)
}

#calculating TPR and FPR
for (h in 1:k) {
    tpr[h]<-length(intersect(which(indicator_est[,h]>=0.5),which(beta_true[,h]!=0)))/length(which(beta_true[,h]!=0))
    fpr[h]<-length(intersect(which(indicator_est[,h]>=0.5),which(beta_true[,h]==0)))/length(which(beta_true[,h]==0))
  }

#calculating TCO
tco = length(which(S_final==S_true))/n 
