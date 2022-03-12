#===================================================#
#                  load package                     #
#===================================================#
list.of.packages = c("philentropy","caret","blockmatrix")
installed_pkgs = installed.packages()[,"Package"]
new.packages = list.of.packages[!(list.of.packages %in% installed_pkgs)]
if(!length(new.packages) == 0){install.packages(new.packages)}
library("philentropy")
library("caret")
library("blockmatrix")
#===================================================#
#                  function                         #
#===================================================#
is.singular.matrix = function(x,tol = 1e-8){
  if( !(nrow(x) == ncol(x)) ) stop("argument x is not a square matrix")
  if( !is.numeric(x) ) stop("argument x is not a numeric matrix")
  return( abs(det(x)) < tol )
}

sigma.spatial = function(period,D,par,q){
  
  if( length(par) == 3) par = c(1,par)
  c = par[1]
  tau = par[2]
  kappa = par[3]
  H = par[4]
  design = matrix(data = NA,nrow = period,ncol = period)
  for(i in c(1:period)){
    design[i,c(i:period)] = c(0:(period-i))
    design[c(i:period),i] = c(0:(period-i))
  }
  k = design %x% matrix(data = 1,nrow = nrow(D)/period,ncol = ncol(D)/period)
  rho = (1/(c*(k**2)+1))*(kappa/(kappa+tau))*exp(-(D/H)^(q))
  diag(rho) = 1
  sigma = rho*kappa
  return(sigma)
  
}

zeropois = function(n,B1,B2,SIGMA1,SIGMA2){
  
  ## Generate X  
  x0 = rep(x = 1,times = n*period)
  x1 = rnorm(n = n*period,mean = 0,sd = 1)
  x2 = rnorm(n = n*period,mean = 0,sd = 1)
  x3 = rnorm(n = n*period,mean = 0,sd = 1)
  x4 = rnorm(n = n*period,mean = 0,sd = 1)
  x5 = runif(n = n*period,min = -1,max = 1)
  x6 = rbinom(n = n*period,size = 1,prob = 0.5)
  X1 = cbind(x0,x1,x2,x3)  #alpha的covariate
  X2 = cbind(x0,x4,x5,x6)  #lambda的covariate
  
  ## alpha、lambda
  alpha = exp(X1 %*% B1)/(1 + exp(X1 %*% B1))
  lambda = exp(X2 %*% B2)
  
  ## 生成判斷結構零之機率
  XX1 = rnorm(n*period)
  XX1.cor = XX1 %*% chol(SIGMA1)
  Y1 = c(qbinom(p = pnorm(XX1.cor),size = 1, prob = alpha))
  
  ## 生成一般卜瓦松變數
  XX2 = rnorm(n*period)
  XX2.cor = XX2 %*% chol(SIGMA2)
  Y2 = c(qpois(p = pnorm(XX2.cor),lambda = lambda))
  
  ## 零膨脹之卜瓦松隨機變數
  Y = ifelse(test = {Y1 == 1},yes = {0},no = Y2)
  
  return(list(X1 = X1,X2 = X2,Y = Y,alpha = alpha,lambda = lambda))
}

newton = function(b,Y,X1,X2,R1,R2,stop.num,D.break = "T"){
  z = 1
  beta00 = matrix(data = c(b),nrow = length(c(b)),ncol = 1)
  repeat{
    
    #===============================================================
    beta1hat = matrix(data = c(beta00[ c(1:length(B1)), ]),ncol = 1)
    beta2hat = matrix(data = c(beta00[-c(1:length(B1)), ]),ncol = 1)
    alphahat = exp(X1 %*% beta1hat)/(1 + exp(X1 %*% beta1hat))
    lambdahat = exp(X2 %*% beta2hat)
    Index = ifelse(test = {Y == 0},yes = {1},no = {0})
    #==================================================================
    mu1 = alphahat + (1-alphahat)*exp(-lambdahat)
    mu2 = (1-alphahat)*lambdahat
    D11 = -t(X1) %*% diag(c((1-exp(-lambdahat))*alphahat*(1-alphahat)))
    D12 =  t(X1) %*% diag(c(alphahat*(1-alphahat)*lambdahat))
    D21 =  t(X2) %*% diag(c((1-alphahat)*exp(-lambdahat)*lambdahat))
    D22 = -t(X2) %*% diag(c((1-alphahat)*lambdahat))
    D = rbind(cbind(D11,D12),cbind(D21,D22))
    #================================================================     
    V1 = diag(x = c(mu1*(1-mu1)))
    V2 = diag(x = c(mu2*(1+lambdahat)-(mu2)^2))
    V1_structural = sqrt(V1) %*% R1 %*% sqrt(V1)
    V2_structural = sqrt(V2) %*% R2 %*% sqrt(V2)
    V = rbind(cbind(V1_structural,
                    matrix(data = 0,nrow = nrow(V1),ncol = ncol(V1))),
              cbind(matrix(data = 0,nrow = nrow(V2),ncol = ncol(V2)),
                    V2_structural))
    diag(V)
    S = D %*% solve(V) %*% rbind((Index-mu1),(Y-mu2))
    I = D %*% solve(V) %*% t(D)
    beta01 = beta00 - solve(I) %*% S
    dif.beta = max(abs(beta01-beta00))
    #============================================
    if( dif.beta < stop.num )break
    if( (z-1) == m.max ){
      assign(x = "count.error.divergent",
             value = {count.error.divergent + 1},
             envir = .GlobalEnv)
      stop("\n","beta is divergent")
    }else{
      beta00 = beta01
      z = z + 1
    }
    #============================================
  }
  
  ## 判斷參數是否在設限區間外 ##
  if(D.break == "T" || D.break == "TRUE"){
    if( sum(abs(beta01)>abs(BB)*breaks) > 0 ){
      assign(x = "count.error.beta",
             value = {count.error.beta + 1},
             envir = .GlobalEnv)
      stop("beta is not under tolerance level")
    }
  }
  return(list(z = z,
              beta = beta01,
              alpha = alphahat,
              lambda = lambdahat,
              mu1 = mu1,mu2 = mu2,
              V1 = V1_structural,
              V2 = V2_structural))
}

spatial = function(par1,par2,result){
  
  result = result
  wq1 = c(log(par1))
  wq2 = c(log(par2))
  alphahat = result$alpha
  lambdahat = result$lambda
  
  # for alpha
  Index = ifelse(test = {Y == 0},yes = {1},no = {0})
  amu = result$mu1
  avar = amu*(1-amu)
  ares = (Index - amu)/sqrt(avar)
  mean(ares);var(ares)
  
  # for lambda
  lmu = result$mu2
  lvar = lmu*(1+lambdahat)-(lmu)^2
  lres = (Y-lmu)/sqrt(lvar)
  mean(lres);var(lres)
  
  # 算N(h)個數
  midpoint = c() # 組距離的組中點
  group = list() # 組距離內的點
  Hmax = max(distance) # 最大距離
  for( i in c(1:kk) ){ # 最大距離分kk組
    minh = (Hmax/kk)*(i-1)
    maxh = (Hmax/kk)*i
    midpoint  = c(midpoint,(minh+maxh)/2)
    group[[i]] = c( which( (distance > minh) & (distance <= maxh) ) )
  } # group[[r]][i] 是位置編號
  sum(sapply(group,length))
  
  site = list() # %/% 是商，%% 是餘數，表示對應位置
  for(r in 1:kk) site[[r]] = cbind(((group[[r]]-1)%/%length(Y) + 1),((group[[r]]-1)%%length(Y) + 1))
  
  # empirical semivariogram
  abc1 = abc2 = list()
  semivar1 = semivar2 = c() 
  for(r in 1:kk){
    abc1[[r]] = abc2[[r]] = numeric()
    if( length(group[[r]]) != 0 ){
      for(i in c(1:length(group[[r]])) ){ 
        abc1[[r]][i] = (ares[ site[[r]][i,1] ] - ares[ site[[r]][i,2] ])^2/(2*length(group[[r]]))
        abc2[[r]][i] = (lres[ site[[r]][i,1] ] - lres[ site[[r]][i,2] ])^2/(2*length(group[[r]]))
      } 
    }
    #加總 empirical semivariogram
    semivar1 = c(semivar1,sum(abc1[[r]])) 
    semivar2 = c(semivar2,sum(abc2[[r]]))
  }
  sapply(abc1,length)
  sum(semivar1);sum(semivar2)
  
  midpoint1 = rep(midpoint,times = period)
  semivar11 = rep(semivar1,times = period)
  semivar21 = rep(semivar2,times = period)
  k = rep(x = c(0:c(period-1)),each = kk)
  
  # 半變異函數 for twos spatial parameter
  RmRR1 = function(a){
    if( period == 1 ) a = c(1,a)
    ll = rep(sapply(abc1,length)/period,times = period)
    rm1 = ll*(semivar11 - (exp(a[2]) + exp(a[3])*(1-(1/(exp(a[1])*k**2+1))*exp(-(midpoint1/exp(a[4]))^(q1)))))^2
    return(sum(rm1))
  }
  RmRR2 = function(a){
    if( period == 1 ) a = c(1,a)
    ll = rep(sapply(abc2,length)/period,times = period)
    rm2 = ll*(semivar21 - (exp(a[2]) + exp(a[3])*(1-(1/(exp(a[1])*k**2+1))*exp(-(midpoint1/exp(a[4]))^(q2)))))^2
    return(sum(rm2))
  }
  
  # alpha(tau,kappa,phi)
  ws1 = optim(par = c(wq1),fn = RmRR1,gr = NULL,method = 'CG')
  ws2 = optim(par = c(wq2),fn = RmRR2,gr = NULL,method = 'CG')
  (s11 = exp(ws1$par))
  (s21 = exp(ws2$par))
  return(list(s11 = s11,s21 = s21))
}

ratio = function(x){
  # 大於1表示spatial參數比較好
  length(which(x > 1))/length(x)
}

newton.L = function(b,Y,X1,X2,R1,R2,stop.num,D.break = "T"){
  z = 1
  beta00 = matrix(data = c(b),nrow = length(c(b)),ncol = 1)
  repeat{
    
    #===============================================================
    beta1hat = matrix(data = c(beta00[ c(1:length(B1)), ]),ncol = 1)
    beta2hat = matrix(data = c(beta00[-c(1:length(B1)), ]),ncol = 1)
    alphahat = exp(X1 %*% beta1hat)/(1 + exp(X1 %*% beta1hat))
    lambdahat = exp(X2 %*% beta2hat)
    Index = ifelse(test = {Y == 0},yes = {1},no = {0})
    #==================================================================
    mu1 = alphahat + (1-alphahat)*exp(-lambdahat)
    mu2 = (1-alphahat)*lambdahat
    D11 = -t(X1) %*% diag(c((1-exp(-lambdahat))*alphahat*(1-alphahat)))
    D12 =  t(X1) %*% diag(c(alphahat*(1-alphahat)*lambdahat))
    D21 =  t(X2) %*% diag(c((1-alphahat)*exp(-lambdahat)*lambdahat))
    D22 = -t(X2) %*% diag(c((1-alphahat)*lambdahat))
    D = rbind(cbind(D11,D12),cbind(D21,D22))
    #================================================================     
    V1 = diag(x = c(mu1*(1-mu1)))
    V2 = diag(x = c(mu2*(1+lambdahat)-(mu2)^2))
    V1_structural = sqrt(V1) %*% R1 %*% sqrt(V1)
    V2_structural = sqrt(V2) %*% R2 %*% sqrt(V2)
    V = rbind(cbind(V1_structural,
                    matrix(data = 0,nrow = nrow(V1),ncol = ncol(V1))),
              cbind(matrix(data = 0,nrow = nrow(V2),ncol = ncol(V2)),
                    V2_structural))
    S = D %*% solve(V) %*% rbind(( Index-mu1 ),( L2%*%Y-mu2 ))
    I = D %*% solve(V) %*% t(D)
    beta01 = beta00 - solve(I) %*% S
    dif.beta = max(abs(beta01-beta00))
    #============================================
    if( dif.beta < stop.num )break
    if( (z-1) == m.max ){
      assign(x = "count.error.divergent",
             value = {count.error.divergent + 1},
             envir = .GlobalEnv)
      stop("\n","beta is divergent")
    }else{
      beta00 = beta01
      z = z + 1
    }
    #============================================
  }
  
  ## 判斷參數是否在設限區間外 ##
  if(D.break == "T" || D.break == "TRUE"){
    if( sum(abs(beta01)>abs(BB)*breaks) > 0 ){
      assign(x = "count.error.beta",
             value = {count.error.beta + 1},
             envir = .GlobalEnv)
      stop("beta is not under tolerance level")
    }
  }
  return(list(z = z,
              beta = beta01,
              alpha = alphahat,
              lambda = lambdahat,
              mu1 = mu1,mu2 = mu2,
              V2 = V2))
}





newton.b1 = function(b1,b2,Y,X1,X2,R){
  
  #=========================================================
  Index = ifelse(test = {Y == 0},yes = {1},no = {0})
  beta1hat = matrix(data = c(b1),ncol = 1)
  beta2hat = matrix(data = c(b2),ncol = 1)
  alphahat = exp(X1 %*% beta1hat)/(1 + exp(X1 %*% beta1hat))
  lambdahat = exp(X2 %*% beta2hat)
  #================================================================
  mu = alphahat + (1-alphahat)*exp(-lambdahat)
  D = -t(X1) %*% diag(c((1-exp(-lambdahat))*alphahat*(1-alphahat)))
  #================================================================     
  V1 = diag(x = c(mu*(1-mu)),nrow = n,ncol = n)
  V1_structural = sqrt(V1) %*% R %*% sqrt(V1)
  V = V1_structural
  S = D %*% solve(V) %*% (Index-mu)
  I = D %*% solve(V) %*% t(D)
  #==============================================================
  beta11 = beta10 - solve(I) %*% S
  #============================================
  # 判斷參數是否在設限區間外
  if( sum(abs(beta11) > abs(B1)*breaks) > 0 )
    stop("beta1 is not under tolerance level")
  return(list(beta = beta11,alpha = alphahat,lambda = lambdahat))
}
newton.b2 = function(b1,b2,Y,X1,X2,R){
  
  #============================================================
  beta1hat = matrix(data = c(b1),ncol = 1)
  beta2hat = matrix(data = c(b2),ncol = 1)
  alphahat = exp(X1 %*% beta1hat)/(1 + exp(X1 %*% beta1hat))
  lambdahat = exp(X2 %*% beta2hat)
  #======================================================================
  mu = (1-alphahat)*lambdahat
  D = -t(X2) %*% diag(c((1-alphahat)*lambdahat))
  #================================================================     
  V2 = diag(x = c(mu*(1+lambdahat)-(mu)^2),nrow = n,ncol = n)
  V2_structural = sqrt(V2) %*% R %*% sqrt(V2)
  V = V2_structural
  S = D %*% solve(V) %*% (Y-mu)
  I = D %*% solve(V) %*% t(D)
  #==============================================================
  beta21 = beta20 - solve(I) %*% S
  
  # 判斷參數是否在設限區間外
  if( sum( abs(beta21) > abs(B2)*breaks ) > 0 )
    stop("beta2 is not under tolerance level")
  
  return(list(beta = beta21,alpha = alphahat,lambda = lambdahat))
}
