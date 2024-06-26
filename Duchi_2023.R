MH_d=function(x,Sigma){
  return(sqrt(matrix(x,ncol=length(x))%*%solve(Sigma)%*%matrix(x,ncol=1)))
}

COVSAFE=function(data, B, mm, z,w){
  n=nrow(data)
  d=ncol(data)
  m=floor(n/2)
  
  # /* pair and rescale */
  y=sapply(1:m,function(i){(data[i,]-data[i+m,])/sqrt(2)})%>%t()
  Sigma_0=matrix(rowMeans(apply(y,1,function(x){x%*%t(x)})),ncol=d)
  con=F
  t=0
  Rt=c()
  Rs=list(Rt)
  Sigmas=list(Sigma_0)
  
  while(!con){
    t=t+1
    truncated=0
    new_loop=(1:m)[!((1:m) %in% Rt)]
    for(i in new_loop){

      LHS=2*log(MH_d(y[i,],Sigma_0))+z[i]+z[m+1]
      if(LHS>log(B)){
        Rt=c(Rt,i)
        truncated=truncated+1
      }
    }
    if(truncated==0){
      con=T
      Tt=t
    }
    if(length(Rt)>0)
      # I am thinking the paper is wrong about this... 
      Sigma_0=matrix(rowMeans(apply(y[-Rt,],1,function(x){x%*%t(x)})),ncol=d)
    Rs=append(Rs,list(Rt))
    Sigmas=append(Sigmas,list(Sigma_0))
  }
  Gamma=list(Rs,Sigmas,Tt)
  if(length(Rt)>(mm+w)){
    print("Potentially unstable")
    return(list(NA,Gamma))
  }
  return(list(Sigma_0,Gamma))
}

set.seed(14)
dd=2
n=5000
data=replicate(dd,rnorm(n,5))
A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
data=data%*%A
A%*%t(A)
norm(A%*%t(A))
t(A)%*%A
tv=(t(A)%*%A)[1,1]

COVSAFE(data, 10, 100, rep(0,n/2+1),0)[[1]]


one_run=function(){
  dd=2
  n=5000
  data=replicate(dd,rnorm(n,5))
  A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
  data=data%*%A
  A%*%t(A)
  norm(A%*%t(A))
  t(A)%*%A
  
  COVSAFE(data, 10, 100, rep(0,n/2+1),0)[[1]][1,1]-tv
}
# res=replicate(100,one_run())
# mean(res)

logdiam=function(X,A){
  n=nrow(X)
  pairs=combn(n,2)
  return(log(max(apply(X[pairs[1,],]-X[pairs[2,],],1,function(x){MH_d(x,A)}))))
}

topK=function(x,k,xi1,xi2){
  y1=x+xi1
  y2=x+xi2
  # R=seq_along(y1)[rank(-y1) <= k]
  x_j=rep(NA,length(x))
  x_j[rank(-y1) <= k]=y2[rank(-y1) <= k]
  return(x_j)
}

MEANSAFE=function(data, A, B,b,k, partm, z,z_prime,w,zN){
  n=nrow(data)
  d=ncol(data)
  m=floor(n/b)
  D=c()
  for(j in 1:m){
    D=c(D,logdiam(data[ partm[[j]],],A))
  }
  Dtilde=topK(D,k,z,z_prime)
  R=c()
  t=0
  for(j in 1:m){
    if(!is.na(Dtilde[j]) && Dtilde[j]>log(sqrt(B)/4)){
      R=c(R,partm[[j]])
      t=t+1
    }
  }
  if(length(R)>0)
    muhat=colMeans(data[-R,])
  else
    muhat=colMeans(data)
  if(t>2*k/3+w)
    mutilde=NA
  else
    mutilde=muhat+expm::sqrtm(A)%*%matrix(zN,ncol=1)
  Gamma=list(D,Dtilde,R,t,muhat)
  return(list(mutilde,Gamma))
}

get_partition=function(n,ngrps){
  perm=sample(1:n,n)
  group_factor = gl(n/ngrps,ngrps,n)
  split_vec = split(perm, group_factor)
  return(split_vec)
}


set.seed(14)
dd=2
n=500
data=replicate(dd,rnorm(n,5))
A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
data=data%*%A
A%*%t(A)
norm(A%*%t(A))
t(A)%*%A

b=20
# (data, A, B,b,k, partm, z,z_prime,w,zN)
colMeans(data)
d=ncol(data)
MEANSAFE(data, t(A)%*%A, 1000000000,20,50, get_partition(n,b), rep(0,n/b),rep(0,n/b),0,rep(0,d))[[1]]

PRIVMEAN=function(data,B,eps,delta){
  n=nrow(data)
  d=ncol(data)
  m=16*log(1/delta)/eps
  m_max=m+16*log((1+exp(eps/4))/delta)/eps
  sigma_z=32*sqrt(exp(1))*B*(m_max+1)/(n*eps)
  sigma_w=16/eps
  Z_cov=jmuOutlier::rlaplace(n/2+1,0,sigma_z)
  W_cov=jmuOutlier::rlaplace(1,0,sigma_w)
  covv=COVSAFE(data,B,m,Z_cov,W_cov)
  if(any(is.na(covv[[1]])))
    return(NA)
  b=1+log(6*n^2/delta ,2)
  k=24 *log(3/delta)/eps-3
  sigma_top=8*k*B*sqrt(exp(1))/((n*eps)*(1-B*sqrt(exp(1))/n))
  #
  t_N_1=exp(3*sigma_top*log(12*n/(b*delta))-log(n*eps))
  sigma_N=(20*b*sqrt(B))*t_N_1
  sigma_w=8/eps
  S=get_partition(n,b)
  z_top=jmuOutlier::rlaplace(floor(n/b),0,sigma_top)
  z_top_p=jmuOutlier::rlaplace(floor(n/b),0,sigma_top)
  w=jmuOutlier::rlaplace(1,0,sigma_w)
  ZN=rnorm(d,0,sigma_N)
  mu_tilde=MEANSAFE(data, covv[[1]], B,b,k, S, z_top, z_top_p,w,ZN)[[1]]
  return(mu_tilde)
  
}

PRIVMEAN(data,B=100000,eps=5,delta=0.1)

ADAMEAN=function(data,eps,delta){
  t=2
  while(T){
    mu=PRIVMEAN(data,B=2^(t-1),eps=eps/t^2,delta/t^2)[[1]]
    t=t+1
    print(t)
    if(!is.na(mu))
      return(mu)
  }
}

ADAMEAN(data,eps=5,delta=0.01)
