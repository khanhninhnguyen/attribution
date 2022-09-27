
                #### Used functions

#### Function that simulates a series with or without outliers

SimulatedSeries <- function(n,P,prob.outliers,size.outliers, rho, theta){
  
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
  if (P==4) {
    outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0)
    Y=rnorm(n,0,1)
    Y[pos.outliers] <- Y[pos.outliers] + rnorm(length(Y[pos.outliers]), mean = 0, sd = 2)
    
    cluster.true[pos.outliers]=2
  }
  if (P==5) {
    outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0)
    Y=rnorm(n,0,1)
    Y[pos.outliers] <-  outliers[which(outliers!=0)]
    
    cluster.true[pos.outliers]=2
  }
  if (P==6) {
    outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0)
    Y=rsn(n = 1000, xi = 0, omega = 1, alpha=0.6)
    Y[pos.outliers] <- Y[pos.outliers] + rnorm(length(Y[pos.outliers]), mean = size.outliers, sd = 1)
    
    cluster.true[pos.outliers]=2
  }
  if (P==7){
    outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0 )
    if(theta == 0 & rho == 0){
      Y=rnorm(n,0,1)
    }else{
      ratio = (1 + 2*theta*rho + theta**2)/(1 - rho**2)
      Y=arima.sim(model = list(ar = rho, ma = theta), n = n0, sd = sqrt(1/ratio) )
    }
    Y[pos.outliers] <- Y[pos.outliers] + rnorm(length(Y[pos.outliers]), mean = 0, sd = 2)
    cluster.true[outliers>0]=2
    cluster.true[outliers<0]=3
  } 
  if (P==8) { # replace tail of normal by tail of the outlier distributio, (0, sigma >1)
    # per.out = sum(pnorm(c(-size.outliers:-1000), mean = 0, sd = 2))*2
    # n1 = round(n*prob.outliers/per.out)
    # Y=rnorm(n,0,1)
    # out.lier = rnorm(n1, mean = 0, sd = 2)
    # n2 = length(which(abs(out.lier)>size.outliers))
    # outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(n2/(2*n),n2/(2*n),1-n2/n))
    Y = rnorm(n,0,1)
    # Y <- Y[which(abs(Y)<=size.outliers)] 

    outliers=sample(c(-size.outliers,size.outliers,0), n, replace=TRUE,prob=c(prob.outliers/2,prob.outliers/2,1-prob.outliers))
    pos.outliers=which(outliers!=0)
    n0 = length(pos.outliers)
    n1 = round(n0/2)
    n2 = n0-n1
    
    # add tails from the crosspoints
    c.p = size.outliers*sqrt(2*log(size.outliers, base = exp(1)))/sqrt(size.outliers**2-1) 
    print(c.p)
    tail.right = rnormt( n1, range = c(c.p, Inf), mu = 0) 
    tail.left = rnormt( n2, range = c(-Inf, -c.p), mu = 0) 
    outlier.reordered <- sample( c(tail.left, tail.right), replace = FALSE)
    # list.replaced = out.lier[which(abs(out.lier)>3)]
    # list.replaced <- list.replaced[order(abs(list.replaced))]
    # Y[pos.outliers] <- list.replaced[1:length(pos.outliers)]
    Y[pos.outliers] <- outlier.reordered
    cluster.true[pos.outliers] = 2
  }

  invisible(list(Y =Y,cluster.true=cluster.true))
}

#### generate random series in tails of normal distribution
rnormt <- function(n, range, mu, s = 1) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
  
}

####
EM.init_imp <- function(Y,P,option.init="CAH"){
  n = length(Y)
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
  n = length(Y)
  
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
  n = length(Y)
  
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
  n = length(Y)
  
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



####
EM.algo_SameMean_PropVariance <- function(Y,phi,P){
  n = length(Y)
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
    phi	= M_Step_SameMean_PropVariance(Y,tau,phi)
    
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
M_Step_SameMean_PropVariance=function(Y,tau,phi){
  P=ncol(tau)
  n=nrow(tau)
  
  a=phi$s2[2]/phi$s2[1]
  
  prop = matrix(0,nrow=1,ncol=P)
  mu = matrix(0,nrow=1,ncol=P)
  s2 = matrix(1,nrow=1,ncol=P)
  
  sum.tau=(tau[,1]+tau[,2]/a)
  
  #1. Proportions
  prop = apply(tau,2,sum)/n
  #2. Means
  mu.same=c()
  mu.same=(Y %*% sum.tau)/sum(sum.tau)
  mu=rep(mu.same,2)
  #3. Variances
  var.same =c()
  var.same =(((Y-mu.same[1])^2) %*% sum.tau)/n
  #4. a?
  a  = c()
  a=(((Y-mu.same[1])^2) %*% tau[,2])/sum(tau[,2]*var.same[1])
  a=max(c(1,a))
  
  s2[1]=var.same
  s2[2]=var.same*a
  
  b    = order(s2)
  s2    = sort(s2)
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

####

GMM_imp = function(P, Y, thres){
  Out.EM.init_imp=EM.init_imp(Y,P=P,option.init="CAH")
  Out.EM_imp =EM.algo_imp(Y,Out.EM.init_imp$phi,P=P,Out.EM.init_imp$Id.cluster1)
  tau_imp = Out.EM_imp$tau 
  cluster_imp0 = apply(tau_imp, 1, which.max)
  main.g = as.numeric(names(sort(table(cluster_imp0),decreasing=TRUE)[1]))
  cluster_imp = sapply(tau_imp[,main.g], function(x) ifelse(x>thres, 1, 2) ) # change here to modify how to choose 
  outliers = which(cluster_imp!=1)
  return(list(outliers = outliers, cluster = cluster_imp))
}


####
GMM_sameMean_uequalvar = function(P, Y, thres){
  dist.Y=dist(Y)
  Clust.cah <- hclust(dist.Y^2, method = "ward.D")  
  classif.outliers <- Mclust(Y,G=P,modelName="V",initialization=Clust.cah)
  
  if (is.null(classif.outliers) == TRUE){
    outliers <- c()
    cluster_SameMean_PropVariance <- rep(1, length(Y))
  } else{
    mu=classif.outliers$parameters$mean
    s2=classif.outliers$parameters$variance$sigmasq
    prop=classif.outliers$parameters$pro
    
    b    = order(s2)
    s2   = sort(s2)
    mu   = mu[b]
    prop = prop[b]
    
    phi=list()
    phi=list(prop=prop,mu=mu,s2=s2)
    # 
    Out.EM_SameMean_PropVariance = EM.algo_SameMean_PropVariance(Y,phi,P)
    phi_SameMean_PropVariance = Out.EM_SameMean_PropVariance$phi
    tau_SameMean_PropVariance =Out.EM_SameMean_PropVariance$tau
    empty_SameMean_PropVariance = Out.EM_SameMean_PropVariance$empty
    dv_SameMean_PropVariance = Out.EM_SameMean_PropVariance$dv
    lvinc_SameMean_PropVariance = Out.EM_SameMean_PropVariance$lvinc
    cluster_imp0 = apply(tau_SameMean_PropVariance, 1, which.max)
    main.g = as.numeric(names(sort(table(cluster_imp0),decreasing=TRUE)[1]))
    cluster_SameMean_PropVariance = sapply(tau_SameMean_PropVariance[,main.g], function(x) ifelse(x>thres, 1, 2) )
    
    outliers = which(cluster_SameMean_PropVariance != 1)
  }
  return(list(outliers = outliers, cluster = cluster_SameMean_PropVariance))
}

testOS <- function(x, thres){
  outliers <- c(rep(1, length(x)))
  prob.z=dnorm(x,0,1)
  f=2
  OS=log(prob.z)^(2*f)
  ScOS=(OS-min(OS))/(max(OS)-min(OS))*10
  rg=which(ScOS>thres)
  outliers[rg] <- 2
  return(list(outliers = rg, cluster = outliers))
} 

