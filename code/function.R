




estimation_alpha = function(Y,facs){
  N=dim(Y)[1]
  Ts=dim(Y)[2]
  p=dim(facs)[1]
  ######### Estimation ###########
  MF=diag(1,Ts)-t(facs)%*%solve(facs%*%t(facs))%*%facs
  temp=rep(1,Ts)
  ealpha=as.vector(solve(t(temp)%*%MF%*%temp))*Y%*%MF%*%temp
  eE=(Y-ealpha%*%t(temp))%*%MF
  eSigma=cov(t(eE))
  
  t_alpha=c()
  for (i in 1:N){
    t_alpha=c(t_alpha,ealpha[i]*sqrt(t(temp)%*%MF%*%temp)*sqrt(Ts-p-1)/sqrt(sum(eE[i,]^2)))
  }
  
  TMF = t(temp)%*%MF%*%temp
  return(list(ealpha=ealpha, eSigma=eSigma, TMF=TMF) )
}


estimate_stat = function(Y_in,facs_in,Y_out,facs_out,lambda){
  N=dim(Y_in)[1]
  
  ######## data split
  result_in = estimation_alpha(Y_in,facs_in)
  result_out = estimation_alpha(Y_out,facs_out)

  ealpha_in = result_in$ealpha
  eSigma_in = result_in$eSigma
  TMF_in =  result_in$TMF
  ealpha_out = result_out$ealpha
  eSigma_out = result_out$eSigma
  TMF_out =  result_out$TMF
  
  ######## SSA
  lambda_in = sum(diag(eSigma_in))/N
  lambda_out = sum(diag(eSigma_out))/N

  ridge_in = lambda*lambda_in
  w_GRS_in = solve(eSigma_in + ridge_in*diag(N))%*%ealpha_in
  
  stat_CGRS_out1 = try(as.vector(sqrt(TMF_out))*t(w_GRS_in)%*%ealpha_out%*%
                         (t(w_GRS_in)%*%eSigma_out%*%w_GRS_in)^{-1/2},silent=TRUE)
  
  ridge_out = lambda*lambda_out
  w_GRS_out = solve(eSigma_out + ridge_out*diag(N))%*%ealpha_out
  
  stat_CGRS_out2 = try(as.vector(sqrt(TMF_in))*t(w_GRS_out)%*%ealpha_in%*%
                         (t(w_GRS_out)%*%eSigma_in%*%w_GRS_out)^{-1/2},silent=TRUE)
  stat_CGRS_out_avg = (stat_CGRS_out1 + stat_CGRS_out2)/2


  return(stat_CGRS_out_avg)
}



boot_t = function(Y,facs,lambda,W, m){
  N=dim(Y)[1]
  Ts=dim(Y)[2]
  p=dim(facs)[1]
  ######### Estimation ###########
  MF=diag(1,Ts)-t(facs)%*%solve(facs%*%t(facs))%*%facs
  temp=rep(1,Ts)
  ealpha=as.vector(solve(t(temp)%*%MF%*%temp))*Y%*%MF%*%temp
  eE=(Y-ealpha%*%t(temp))%*%MF
  
  facs_T = rbind(t(temp), facs)
  beta_all = solve(facs_T%*%t(facs_T))%*%(facs_T%*%t(Y))

  set.seed(m)
  v1 <- matrix(2 * (runif(N * Ts) > 0.5) - 1, nrow = N, ncol = Ts)
  Y_boot = t(beta_all[-1,,drop = FALSE])%*%facs + v1*eE
    
    
  W_sample = sample(1:Ts, W)
  Y_in = Y_boot[,c(W_sample)]
  Y_out = Y_boot[,-c(W_sample)]
    
  facs_in = facs[,c(W_sample),drop = FALSE]
  facs_out = facs[,-c(W_sample),drop = FALSE]
    
  Y1 = cbind(Y_in, Y_out)
  facs1 = cbind(facs_in, facs_out)

  stat = try(estimate_stat(Y_in,facs_in,Y_out,facs_out,lambda),silent=TRUE)

  return(as.numeric(stat))
}

RGRS <- function(R, R_M, z) {
  
  P <- nrow(R)
  T <- ncol(R)
  K <- nrow(R_M)
  c <- P / T
  R_bar <- rowMeans(R)
  R_M_bar <- rowMeans(R_M)
  
  X <- t(cbind(1, t(R_M)))  
  beta_hat <- R %*% t(X) %*% solve(X %*% t(X))
  beta_hat <- beta_hat[, -1]  
  
  alpha_hat <- R - beta_hat %*% R_M
  alpha_bar <- rowMeans(alpha_hat)
  
  Sigma_hat <- (alpha_hat %*% t(alpha_hat)) / T - alpha_bar %*% t(alpha_bar)
  
  Sigma_hat_z <- z * diag(P) + Sigma_hat
  W <- solve(Sigma_hat_z)
  m_hat <- sum(diag(W)) / P
  
  xi_hat <- 1 / (1 - c + c * z * m_hat) - 1
  
  m_hat_prime <- sum(diag(W %*% W)) / P
  xi_hat_prime_num <- c * ( - m_hat + z * m_hat_prime ) 
  
  xi_hat_prime_denom <- (1 - c + c * z * m_hat)^2
  xi_hat_prime <- xi_hat_prime_num / xi_hat_prime_denom
  
  Delta_hat <- 2 * (xi_hat + z * xi_hat_prime) * (1 + xi_hat)^2
  
  Sigma_F <- cov(t(R_M))  
  theta_M_squared <- t(R_M_bar) %*% solve(Sigma_F) %*% R_M_bar
  
  theta_alpha_squared_z <- t(alpha_bar) %*% solve(Sigma_hat_z) %*% alpha_bar
  
  W_z <- theta_alpha_squared_z / (1 + theta_M_squared)
  
  test_stat <- sqrt(T) * (W_z - xi_hat) / sqrt(Delta_hat)
  
  return(list(test_stat = as.numeric(test_stat),
              W_z = as.numeric(W_z),
              xi_z = xi_hat,
              Delta_z = Delta_hat))
}


Alpha_test = function(facs,Y,lambda){
 # Y: N*T
 # facs: K*T
  
  N=dim(Y)[1]
  Ts=dim(Y)[2]
  p=dim(facs)[1]
  ######### Estimation ###########
  MF=diag(1,Ts)-t(facs)%*%solve(facs%*%t(facs))%*%facs
  temp=rep(1,Ts)
  ealpha=as.vector(solve(t(temp)%*%MF%*%temp))*Y%*%MF%*%temp
  eE=(Y-ealpha%*%t(temp))%*%MF
  eSigma=cov(t(eE))

  t_alpha=c()
  for (i in 1:N){
    t_alpha=c(t_alpha,ealpha[i]*sqrt(t(temp)%*%MF%*%temp)*sqrt(Ts-p-1)/sqrt(sum(eE[i,]^2)))
  }
  
  ######### Statistics ###########
  
    ### SSA test
    W = 1/2*Ts
    Y_in = Y[,1:W]
    Y_out = Y[,-c(1:W)]

    facs_in = facs[,1:W,drop = FALSE]
    facs_out = facs[,-c(1:W),drop = FALSE]
    lambda = 1
    SSA = estimate_stat(Y_in,facs_in,Y_out,facs_out,lambda)

    clnum <- parallel::detectCores()
    cl <- parallel::makeCluster(max(1, floor(clnum * 0.6)))
    doParallel::registerDoParallel(cl)

    parallel::clusterExport(cl, 
                            c("Y", "facs", "lambda", "W", 
                              "boot_t","estimation_alpha", "estimate_stat"),
                            envir = environment())

    # bootstrap
    boot_res <- foreach::foreach(m = 1:300, 
                                 .combine = "c",
                                 .packages = c("doParallel"),
                                 .errorhandling = "pass") %dopar% {
                                   tryCatch({
                                     boot_t(Y, facs, lambda, W, m)
                                   }, error = function(e) {
                                     message("Error in iteration ", m, ": ", e$message)
                                     return(NA)
                                   })
                                 }
    
    parallel::stopCluster(cl)

    p_SSA <- mean(na.omit(boot_res) > as.numeric(SSA))


  
  #### MAX
  stat_max=max(t_alpha^2) - 2*log(N)+log(log(N))
  cval_max=-log(pi)-2*log(-log(1-0.05))
  p_max = 1-exp(-pi^{-1/2}*exp(-stat_max/2))
  
  #### PY
  v=Ts-p-1
  rho2<-0
  Rut<-cor(t(eE))
  pn<-0.1/(N-1)
  thetan<-(qnorm(1-pn/2))^2
  rho2<-(sum((Rut[Rut^2*v>thetan])^2)-N)/2
  rhos<-rho2*2/N/(N-1)
  
  stats = t_alpha^2
  stat=N^(-1/2)*sum(stats-v/(v-2))/
    (v/(v-2)*sqrt(2*(v-1)/(v-4)*(1+N*rhos)))
  
  stat_PY = stat
  cval_PY =qnorm(0.95)
  p_PY<-1-pnorm(stat_PY)
  
  Tcom<-1-(p_PY>0.025)*(p_max>0.025)
  pcom=(p_PY+p_max)/2
  
  
  ############## CGRS
  stat_CGRS = RGRS(Y, facs, z = 10)
  stat_CGRS = stat_CGRS$test_stat
  p_CGRS = 2*(1 - pnorm(abs(stat_CGRS))) 
  cval_CGRS = qnorm(0.975)
  
  res1=c(p_SSA,p_max, p_PY,pcom,p_CGRS)

  
  return(res1)
}


