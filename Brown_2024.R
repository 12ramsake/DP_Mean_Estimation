library(dplyr)

# https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
# is_inv = function(m){class(try(solve(m),silent=T))=="matrix"}

is_inv = function(X) !inherits(try(solve(X), silent = TRUE), "try-error")
#y is m x d
LargestGoodSubset=function(y, lambda){
  
  m=nrow(y)
  d=ncol(y)
  S=1:m
  
  new_S=S>0
  # out=S
  # 
  # while(length(out)!=0 && length(S)>0){
  #   mat=matrix(0,ncol=d,nrow=d)
  #   for(i in S){
  #     mat=mat+matrix(y[i,],ncol=1)%*%matrix(y[i,],ncol=d)/m
  #   }
  #   inv_mat=expm::sqrtm(solve(mat))
  # 
  #   out=c()
  #   for(i in S){
  #     if(sum((inv_mat%*%matrix(y[i,],ncol=1))^2)>lambda)
  #       out=c(out,i)
  #   }
  #   if(length(out)>0)
  #     S=S[!(S %in% out)]
  # }
    
  while((sum(new_S)>0) & (length(S)!=0)){
    # print(sum(new_S))
    if(length(S)>1){
      YYt=matrix(rowMeans(apply(y[S,],1,function(x){x%*%t(x)})),ncol=d)
    }
    else{
      YYt=matrix(y[S,],ncol=1)%*%matrix(y[S,],ncol=d)
    }
    if(is_inv(YYt)){
      inverted=expm::sqrtm(solve(YYt))
      new_S=apply(matrix(y[S,],nrow=length(S),ncol=d),1,function(x){sum((inverted%*%matrix(x,ncol=1))^2)>lambda})
      S=S[!new_S]
    }
    else
      break
  }
  
  return(S)
}





StableCovariance=function(data, lambda_0, k){
  
  n=nrow(data)
  d=ncol(data)
  m=floor(n/2)
  
  # /* pair and rescale */
  y=sapply(1:m,function(i){(data[i,]-data[i+m,])/sqrt(2)})%>%t()
  
  
  S_ell=sapply(0:(2*k), function(i){LargestGoodSubset(y,lambda_0*exp(i/k))})
  

  score=min(c(k,unlist(lapply(1:(k+1),function(x){m-length(S_ell[[x]])+x-1}))))
  
  # w=c()
  # for(i in 1:m){
  #   count=0
  #   for(j in (k+1):(2*k)){
  #     count=count+(i %in% S_ell[[j+1]])/(k*m)
  #   }
  #   w=c(w, count)
  # }
  # 
  w=sapply(1:m,function(x){ lapply(S_ell[(k+2):(2*k+1)],function(z){x %in% z})%>%unlist()%>%sum() })/(k*m)
  
  # summ=0
  # for(i in 1:m){
  #   summ=summ+w[i]*matrix(y[i,],ncol=1)%*%matrix(y[i,],ncol=d)
  # }
  
  
  y_w=apply(y,2,'*',sqrt(w))
  Sigma=matrix(rowSums(apply(y_w,1,function(x){x%*%t(x)})),ncol=d)
  
  return(list("Sigma"=Sigma,"score"=score))
} 





LargestCore=function(data, Sigma, lambda, tau,R){
  
  n=nrow(data)
  inv=solve(Sigma)
  test=function(i){
    step_1=apply(-data[R,],1,'+',data[i,])%>%t()
    step_2=diag(step_1%*%solve(Sigma)%*%t(step_1))
    return(R[step_2<=lambda])
  }

  N_i=sapply(1:nrow(data),test)
  step_3=lapply(N_i,function(x){length(x)>tau})%>%unlist()
  ind=1:nrow(data)
  return(ind[step_3])
  # 
  # n=nrow(data)
  # N=c()
  # N_is=list()
  # for(i in 1:n){
  #   N_i_count=0
  #   N_i=c()
  #   for(j in R){
  #     add=mahalanobis(data[i,] , data[j,], cov=inv, inverted = TRUE)^2<=lambda
  #     N_i_count=N_i_count+add
  #     if(add)
  #       N_i=c(N_i,j)
  #   }
  #   N_is=append(N_is,list(N_i))
  #   N=c(N,N_i_count)
  # }
  # ind=1:n
  # 
  # return(ind[N>=tau])
  
}

StableMean=function(data, Sigma, lambda_0, k,R){
  
  n=nrow(data)
  d=ncol(data)
  
  S_ell=sapply(0:(2*k),function(i){LargestCore(data, Sigma, lambda_0*exp(i/k), length(R)-i , R)})
  
  
  
  scores=sapply(0:k,function(i){n-length(S_ell[[i+1]])+i})
  score=min(c(k,min(scores)))
  c_i= sapply(1:n,function(i){sum(sapply((k+1):(2*k),function(ell){i %in% S_ell[[ell+1]]}))})
  Z=sum(c_i)
  w_i=c_i/Z
  mu_hat=colSums(apply(data,2,'*',w_i))
  return(list("muhat"=mu_hat,"score"=score))
}


PTR=function(score,eps,delta){
  tau=2*log((1-delta)/delta)/eps+4
  pass_prob=1-delta*exp(eps*(score-2)/2)
  if(score==0)
    return(TRUE)
  else if(score>=tau)
    return(FALSE)
  else
    return(rbinom(1,1,pass_prob)==1)
}

GB_ME=function(data,eps,delta,lambda=5){
  
  n=nrow(data)
  d=ncol(data)
  
  requirement=n >=192*exp(2)*lambda*log(6/delta)/eps+160*exp(2)*lambda 
  if(!requirement){
    print("Requirement n >=192*exp(2)*lambda*log(6/delta)/eps+160*exp(2)*lambda not met, decrease lambda or adjust privacy budget ")
    print(paste0("n",n))
    print(paste0("RHS",192*exp(2)*lambda*log(6/delta)/eps+160*exp(2)*lambda ))
    return(NULL)
  }
  k=ceiling((6*log(6/delta))/eps)+4
  M=6*k+ceiling(18*log(16*n/delta))
  c_sq=720*exp(2)*lambda*log(12/delta)/(eps^2*n^2)
  R=sample(1:n,M)
  
  np_ce=StableCovariance(data, lambda, k)
  np_ce$Sigma
  # np_ce$Sigma=lqmm::make.positive.definite(np_ce$Sigma)
  np_me=StableMean(data,np_ce$Sigma, lambda, k, R)
  overall_score=max(c(np_ce$score,np_me$score))
  test_result=PTR(overall_score,eps/3,delta/6)
  if(test_result){
    return(MASS::mvrnorm(1,mu=np_me$mu_hat,Sigma=c_sq*np_ce$Sigma))
  }
  else{
    print("FAILED")
    return(NA)
  }
  
}

GB_ME_np=function(data,eps,delta,lambda=5){
  
  n=nrow(data)
  d=ncol(data)
  
  requirement=n >=192*exp(2)*lambda*log(6/delta)/eps+160*exp(2)*lambda 
  if(!requirement){
    print("Requirement n >=192*exp(2)*lambda*log(6/delta)/eps+160*exp(2)*lambda not met, decrease lambda or adjust privacy budget ")
    print(paste0("n",n))
    print(paste0("RHS",192*exp(2)*lambda*log(6/delta)/eps+160*exp(2)*lambda ))
    return(NULL)
  }
  k=ceiling((6*log(6/delta))/eps)+4
  M=6*k+ceiling(18*log(16*n/delta))
  c_sq=720*exp(2)*lambda*log(12/delta)/(eps^2*n^2)
  R=sample(1:n,M)
  
  np_ce=StableCovariance(data, lambda, k)
  np_ce$Sigma
  # np_ce$Sigma=lqmm::make.positive.definite(np_ce$Sigma)
  np_me=StableMean(data,np_ce$Sigma, lambda, k, R)
  overall_score=max(c(np_ce$score,np_me$score))
  return(np_me$muhat)
  
}






set.seed(14)
dd=2
n=5000
data=replicate(dd,rnorm(n,5))
# A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
# data=data%*%A
# A%*%t(A)
# norm(A%*%t(A))
# t(A)%*%A

colMeans(data)
# overall_dp_mean(data)


GB_ME(data,5,0.1,1)



req=function(lambda){
  return(192*exp(2)*lambda*log(6/0.1)/5+160*exp(2)*lambda )
}

curve(req(x),1,10,xlab="lambda",ylab='n')


GB_ME_np(data,5,0.1,2)
