

## COIN PRESS MEAN ESTIMATION
norm_2=function(x){sqrt(sum(x^2))}

MVM=function(data,C,r,rho,beta){
  d=ncol(data)
  n=nrow(data)
  gamma_1=sqrt(d+2*sqrt(d*log(n/beta))+2*log(n/beta))
  gamma_2=sqrt(d+2*sqrt(d*log(1/beta))+2*log(1/beta))
  
  #project data
  new_r=r+gamma_1
  centered=t(apply(data,1,'-',C))
  norms=apply( centered ,1,norm_2)
  data_proj=data
  outside_ball=norms>new_r
  if(any(outside_ball)){
    # print(outside_ball)
    if(sum(outside_ball)!=1)
      data_proj[outside_ball,]=apply(centered[outside_ball,]*new_r/norms[outside_ball],1,'+',C)
    else
      data_proj[outside_ball,]=centered[outside_ball,]*new_r/norms[outside_ball]+C
  }
  
  
  Delta=2*new_r/n
  
  Z=colMeans(data)+rnorm(d,0,sd=sqrt(Delta^2/(2*rho)))
  r_prime=gamma_2*sqrt(1/n+2*new_r^2/(n^2*rho))
  return(c(Z,r_prime))
}



# The sum of the rho vector is rho-CDP
# We just need to give nice initial values for C and r
# Beta is the probability that we want to be accurate within... super weird... 
MM_est=function(data,C,r,t,rho_vec,beta=0.01){
  
  d=ncol(data)
  n=nrow(data)
  
  for(i in 1:(t-1)){
    new_Cr=MVM(data,C,r,rho_vec[i],beta/(4*(t-1)))
    C=new_Cr[1:d]
    r=new_Cr[d+1]
  }
  mean_est=MVM(data,C,r,rho_vec[t],beta/4)[1:d]
  return(mean_est)
}

dd=10
data=replicate(dd,rnorm(100,4))

MM_est(data,rep(0,dd),100,5,rep(1/5,5))
colMeans(data)



### Covariance 
# d=2
# n=10
# beta=0.1
# data=matrix(1:20,ncol=2,byrow = T)
# A=diag(d)/sqrt(100)
# data
# W=data%*%A
# gamma=sqrt(d+2*sqrt(d*log(n/beta))+2*log(n/beta))
# eta = 0.5*(2*(sqrt(d/n)) + (sqrt(d/n))**2)
# W
# 
# norms=apply(W ,1,norm_2)
# # norms
# data_proj=W
# outside_ball=norms>gamma
# data_proj[outside_ball,]=data_proj[outside_ball,]*gamma/norms[outside_ball]
# data_proj
# 
# 
# 
# Z=t(data_proj)%*%data_proj/n
# 
# eta = 0.5*(2*(sqrt(d/n)) + (sqrt(d/n))**2)
# nu=0
# 
# U=Z+(eta+nu)*diag(nrow=d)
# 
# # print(solve(U))
# 
# U=lqmm::make.positive.definite(U)
# sqrt_m=expm::sqrtm(solve(U))
# # print(sqrt_m)
# A_prime=sqrt_m%*%A
# A_prime

MVC=function(data,A,rho,beta){
  
  d=ncol(data)
  n=nrow(data)

  
  W=data%*%t(A)
  # W=data%*%A
  gamma=sqrt(d+2*sqrt(d*log(1/beta))+2*log(1/beta))
  
  #project data
  norms=apply(W ,1,norm_2)
  data_proj=W
  outside_ball=norms>gamma
  if(any(outside_ball))
    data_proj[outside_ball,]=data_proj[outside_ball,]*gamma/norms[outside_ball]
  
  Delta=sqrt(2)*gamma^2/n
  
  noise=matrix(0,d,d)
  num_g=d^2-(d^2-d)/2
  noise[lower.tri(noise,T)]=rnorm(num_g,0,sd=sqrt(Delta^2/(2*rho^2)))
  noise[upper.tri(noise)] = t(noise)[upper.tri(noise)]
  Z=t(data_proj)%*%data_proj/n+noise
  
  # Using what is in their code....
  # eta=2*sqrt(d/n)+2*sqrt(2*log(beta/2)/n)+(sqrt(d/n)+sqrt(2*log(beta/2)/n))^2
  # nu_1=6*(1+sqrt(log(d))*((log(d)/d)^(1/3)))/sqrt(log(1+(log(d)/d)^(1/3)))
  # nu=(gamma^2/(n*sqrt(rho)))*(  2*sqrt(d)+2*d^(1/6)*(log(d))^(1/3)+ nu_1+2*sqrt(2*log(1/beta)) )
  # 
  eta = 0.5*(2*(sqrt(d/n)) + (sqrt(d/n))**2)
  nu=0
  
  Z=lqmm::make.positive.definite(Z)
  U=Z+(eta+nu)*diag(nrow=d)
  
  # print(solve(U))
  
  # U=lqmm::make.positive.definite(U)
  U_inv=solve(U)
  # SVD=svd(U_inv)
  sqrt_m=expm::sqrtm(U_inv)
  
  # print(sqrt_m)
  A_prime=sqrt_m%*%A
  # print(A_prime)
  return(list("A"=A_prime,"Est"=Z))

}


# The sum of the rho vector is rho-CDP
# We just need to give nice initial values for C and r
# Beta is the probability that we want to be accurate within... super weird... 
MVCRec=function(data,K,t,rho_vec,beta=0.1){

  
  d=ncol(data)
  n=nrow(data)
  new_A=list()
  new_A$A=diag(d)/sqrt(K)
  A_0=new_A$A
  
  for(i in 1:(t-1)){
    new_A=MVC(data,A_0,rho_vec[i],beta/(4*(t-1)))
    # print(new_A$A)
    A_0=new_A$A
  }
  
  cov_est=MVC(data,new_A$A,rho_vec[i],beta/(4*(t-1)))$Est
  A_inv=solve(new_A$A)
  cov_est=A_inv%*%cov_est%*%A_inv
  lqmm::make.positive.definite(cov_est)
  return(cov_est)
}

set.seed(14)
dd=2
n=1000
data=replicate(dd,rnorm(n,0))
MVCRec(data,10,3,10*rep(1/5,3),0.1)
cov(data)


set.seed(14)
dd=20
n=1000
data=replicate(dd,rnorm(n,0))
A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
data=data%*%A
A%*%t(A)
norm(A%*%t(A))
t(A)%*%A
MVCRec(data,1000,5,5*rep(1,5),0.1)
cov(data)


norm_2(MVCRec(data,1000,5,5*rep(1,5),0.1)-cov(data))

overall_dp_mean=function(data,C=0,r=100,t_1=2,t_2=5,
                         rho_vec_1=5*rep(1/2,2),rho_vec_2=5*rep(1/3,5),K=NULL,beta=0.1){
  
  # difference
  if(is.null(K))
    K=50*sqrt(ncol(data))
  
  diffed=data[seq(1,nrow(data),2),]-data[seq(2,nrow(data),2),]
  
  
  cov_est=MVCRec(diffed,K,t_2,rho_vec_2,beta)/2
  # cov_est
  # SVD=svd(cov_est)
  # D=diag(sqrt(SVD$d))
  # sqrt_mat=SVD$u%*%D%*%t(SVD$v)
  # solve(sqrt_mat)
  # expm::sqrtm(solve(cov_est))
  
  
  
  adj=expm::sqrtm(solve(cov_est))
  
  whitened=data%*%adj
  mean_est=MM_est(whitened,C,r,t_1,rho_vec_1,beta)
  mean_est=mean_est%*%expm::sqrtm(cov_est)
  return(mean_est)
}

# 
# lazy_dp_mean=function(data,C=0,r=100,t_1=2,rho_vec_1=2*rep(1/2,2),beta=0.1){
#   
#   mean_est=MM_est(data,C,r,t_1,rho_vec_1,beta)
#   return(mean_est)
# }



set.seed(14)
dd=2
n=10000
data=replicate(dd,rnorm(n,5))
A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
data=data%*%A
A%*%t(A)
norm(A%*%t(A))
t(A)%*%A

colMeans(data)
overall_dp_mean(data)

cov(data)


set.seed(14)
dd=20
n=10000
data=replicate(dd,rnorm(n,5))
A=matrix(rnorm(dd^2),ncol=dd,nrow=dd)
data=data%*%A
A%*%t(A)
norm(A%*%t(A))
t(A)%*%A

colMeans(data)
overall_dp_mean(data)

