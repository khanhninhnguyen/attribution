
                #### Used functions

#### Function that simulates a series with or without outliers

SimulatedSeries <- function(n,P,prob.outliers,size.outliers){
  
  cluster.true=rep(1,n)
  
  if (P ==1){
    Y=rnorm(n,0,1)
  } 
  if (P==2) {
   outliers=sample(c(size.outliers,0), n, replace=TRUE,prob=c(prob.outliers,1-prob.outliers))
   pos.outliers=which(outliers!=0)
   Y=rnorm(n,0,1)+outliers
   
   cluster.true[pos.outliers]=2
  }
  if (P==3){
   outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
   pos.outliers=which(outliers!=0 )
   Y=rnorm(n,0,1)+outliers
   
   cluster.true[outliers>0]=2
   cluster.true[outliers<0]=3
 }

  invisible(list(Y =Y,cluster.true=cluster.true))
}

####
EM.init_imp <- function(Y,P,option.init="CAH"){
  
  # initialisation avec CAH
  if (option.init=="CAH"){
    dist.Y=dist(Y)
    Clust.cah<- hclust(dist.Y^2, method = "ward.D")
    cluster.init <- cutree(Clust.cah, k =P)
  }
  
  # initialisation with k.means
  if (option.init=="k.means"){
    res.kmeans <- kmeans(Y,centers =P,nstart = 20) 
    cluster.init <- res.kmeans$cluster
  }

  prop = matrix(0,nrow=1,ncol=P)
  mu = matrix(0,nrow=1,ncol=P)
  s2 = matrix(1,nrow=1,ncol=P)
  
  Y.cluster=data.frame(Y=Y,cluster.init=cluster.init)
  res.stat <- Y.cluster %>% group_by(cluster.init) %>% summarise(mu=mean(Y),prop=n())
  mu=res.stat$mu
  prop=res.stat$prop/n
 
  Id.cluster1=which.min( abs(mu))
  #Id.cluster1
  mu[Id.cluster1]=0
  phi=list(prop=prop,mu=mu,s2=s2)
  
  invisible(list(phi = phi,Id.cluster1=Id.cluster1))
  
}


####
EM.init <- function(Y,P,option.init="CAH"){
  
  # initialisation avec CAH
  if (option.init=="CAH"){
    dist.Y=dist(Y)
    Clust.cah<- hclust(dist.Y^2, method = "ward.D")
    cluster.init <- cutree(Clust.cah, k =P)
  }
  
  # initialisation with k.means
  if (option.init=="k.means"){
    res.kmeans <- kmeans(Y,centers =P,nstart = 20) 
    cluster.init <- res.kmeans$cluster
  }
  
  prop = matrix(0,nrow=1,ncol=P)
  mu = matrix(0,nrow=1,ncol=P)
  s2 = matrix(1,nrow=1,ncol=P)
  
  Y.cluster=data.frame(Y=Y,cluster.init=cluster.init)
  res.stat <- Y.cluster %>% group_by(cluster.init) %>% summarise(mu=mean(Y),prop=n(),s2=var(Y))
  mu=res.stat$mu
  prop=res.stat$prop/n
  #s2=res.stat$s2
  
  b    = order(mu)
  mu    = sort(mu)
  s2    = s2[b]
  prop = prop[b]
  
  phi=list(prop=prop,mu=mu,s2=s2)
  
  invisible(list(phi = phi))
  
}




####
EM.algo_imp <- function(Y,phi,P,Id.cluster1){
  
  ### Algo
  delta = 1
  empty = 0
  dv    = 0
  
  iter  = 0
  eps   = 1e-6
  tau   = matrix(1,nrow = n,ncol = P)
  np    = apply(tau,2,sum)
  
  while ((delta>=1e-4) &  (min(np)>eps) & (iter<=5000)){
    iter=iter+1
    phi_temp=phi
    log_dens=log_densite(Y,phi,P)
    
    # E step
    Estepout   = E_Step(phi,log_dens)
    tau        = Estepout$tau
    lvinc      = Estepout$lvinc
    
    # M step
    phi	= M_Step_imp(Y,tau,Id.cluster1)

    # 
    np = apply(tau,2,sum)
    rg         = which(unlist(phi)!=0)
    delta      = max(abs(unlist(phi_temp)[rg]-unlist(phi)[rg])/unlist(phi)[rg])
  } 
  
  
  if (min(np)<eps){
    empty = 1
    lvinc = -Inf
  }
  
  if (iter>5000){
    dv    = 2
    lvinc = -Inf
  }
  rm(delta,log_dens)
  invisible(list(phi = phi,tau = tau,lvinc = lvinc,empty = empty,dv = dv))
  
}


####
EM.algo <- function(Y,phi,P){
  
  ### Algo
  delta = 1
  empty = 0
  dv    = 0
  
  iter  = 0
  eps   = 1e-6
  tau   = matrix(1,nrow = n,ncol = P)
  np    = apply(tau,2,sum)
  
  while ((delta>=1e-4) &  (min(np)>eps) & (iter<=5000)){
    iter=iter+1
    phi_temp=phi
    log_dens=log_densite(Y,phi,P)
    
    # E step
    Estepout   = E_Step(phi,log_dens)
    tau        = Estepout$tau
    lvinc      = Estepout$lvinc
    
    # M step
    phi	= M_Step(Y,tau)
    
    # 
    np = apply(tau,2,sum)
    rg         = which(unlist(phi)!=0)
    delta      = max(abs(unlist(phi_temp)[rg]-unlist(phi)[rg])/unlist(phi)[rg])
  } 
  
  
  if (min(np)<eps){
    empty = 1
    lvinc = -Inf
  }
  
  if (iter>5000){
    dv    = 2
    lvinc = -Inf
  }
  rm(delta,log_dens)
  invisible(list(phi = phi,tau = tau,lvinc = lvinc,empty = empty,dv = dv))
  
}



#####
log_densite=function(Y,phi,P){
  
  mu  = phi$mu
  s2  = phi$s2
  log_dens = matrix(-Inf,ncol=P,nrow=1)
  
  log_dens=sapply(1:P, function(p){
    -0.5*log(2*pi*s2[p])- 0.5*((Y  - mu[p])^2)/s2[p]
  } )
  invisible(log_dens)
} 

####
M_Step_imp=function(Y,tau,Id.cluster1){
  P=ncol(tau)
  n=nrow(tau)
  prop = matrix(0,nrow=1,ncol=P)
  mu = matrix(0,nrow=1,ncol=P)
  s2 = matrix(1,nrow=1,ncol=P)

  prop = apply(tau,2,sum)/n
  
  seq.P=1:P
  for (p in (seq.P[-Id.cluster1])){
    mu[p]  = (Y %*% tau[,p])/sum(tau[,p])
    #s2[p] = (((Y-mu[p])^2) %*%tau[,p])/sum(tau[,p])  
  }   
  
  phi=list(prop=prop,mu=mu,s2=s2)
  invisible(phi)
}

####
M_Step=function(Y,tau){
  P=ncol(tau)
  n=nrow(tau)
  prop = matrix(0,nrow=1,ncol=P)
  mu = matrix(0,nrow=1,ncol=P)
  s2 = matrix(1,nrow=1,ncol=P)
  s=c()
  prop = apply(tau,2,sum)/n
  
  seq.P=1:P
  for (p in seq.P){
    mu[p]  = (Y %*% tau[,p])/sum(tau[,p])
    #s[p]= ((Y-mu[p])^2) %*% tau[,p]  
    #s2[p] = (((Y-mu[p])^2) %*%tau[,p])/sum(tau[,p])  
  }   

  #s2 = rep(sum(s)/n,P)

  
  b    = order(mu)
  mu    = sort(mu)
  s2    = s2[b]
  prop = prop[b]

  phi=list(prop=prop,mu=mu,s2=s2)
  invisible(phi)
}


####
E_Step=function(phi,log_dens){
  
  n = nrow(log_dens)
  P = ncol(log_dens)
  
  tau     = matrix((log( phi$prop)),n,P,byrow=TRUE)+log_dens
  tau_max = apply(tau,1,max)
  tau     = exp(tau-matrix(tau_max,n,P))
  lvinc   = sum(log( apply(tau,1,sum)) + tau_max) 
  tau     = sweep(tau,1,STATS = apply(tau,1,sum), FUN="/")
  
  invisible(list(tau=tau,lvinc=lvinc)) 
}


