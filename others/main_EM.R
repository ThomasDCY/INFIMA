#### Full EM model implementation ####
# res <- ModelFitting(D.g, p.g, B.g.t, Y.g, PI)
ModelFitting <- function(.D.g, .p.g, .B.g.t, .Y.g, .PI){
  G <- length(.D.g)
  
  ####### The counts for the causal generative model
  # alpha <=> a1, beta <=> b1
  # preprocess the count data
  n0 <- vector('list', G)
  n1 <- vector('list', G)
  n2 <- vector('list', G)
  m0 <- vector('list', G)
  m1 <- vector('list', G)
  m2 <- vector('list', G)
  
  # the following are weighted sum counts with Z.g
  N0 <- rep(0,G)
  N1 <- rep(0,G)
  N2 <- rep(0,G)
  M0 <- rep(0,G)
  M1 <- rep(0,G)
  M2 <- rep(0,G)
  
  for(g in 1:G){
    if(.p.g[g] > 0){
      D1.g <- .D.g[[g]][,c(1,2,3,5,6)]
      D2.g <- .D.g[[g]][,c(4,7,8)]
      if(.p.g[g] == 1){
        D1.g <- t(as.matrix(D1.g))
        D2.g <- t(as.matrix(D2.g))
      }
      n0[[g]] <- apply(D1.g, 1, function(x){sum(x==0)})
      n1[[g]] <- apply(D1.g, 1, function(x){sum(x==1)})
      n2[[g]] <- apply(D1.g, 1, function(x){sum(x==2)})
      m0[[g]] <- apply(D2.g, 1, function(x){sum(x==0)})
      m1[[g]] <- apply(D2.g, 1, function(x){sum(x==1)})
      m2[[g]] <- apply(D2.g, 1, function(x){sum(x==2)})
    }
  }
  
  
  
  ####### The counts for the null generative model
  # alpha.0 <=> a0, beta.0 <=> b0
  N00 <- rep(0,G)
  N01 <- rep(0,G)
  N02 <- rep(0,G)
  M00 <- rep(0,G)
  M01 <- rep(0,G)
  M02 <- rep(0,G)
  
  for(g in 1:G){
    if(.p.g[g] > 0){
      tmp <- abs(.B.g.t[[g]] - .Y.g[[g]])
      tmp1 <- tmp[c(1,2,3,5,6)]
      tmp2 <- tmp[c(4,7,8)]
  
      N00[g] <- sum(tmp1 == 0)
      N01[g] <- sum(tmp1 == 1)
      N02[g] <- sum(tmp1 == 2)
      M00[g] <- sum(tmp2 == 0)
      M01[g] <- sum(tmp2 == 1)
      M02[g] <- sum(tmp2 == 2)
    }
  }
  
  
  ##### the reward terms for enforcing concordance
  lambda_0 <- 0.1
  lambda_1 <- 0.01
  lambda_2 <- 0
  
  ### initial values
  alpha.init <- c(0.6,0.3,0.1) # cannot contain any probability equal to zero
  beta.init <- c(0.5,0.4,0.1)
  alpha <- alpha.init
  beta <- beta.init
  alpha.0.init <- c(0.6,0.3,0.1) # cannot contain any probability equal to zero
  beta.0.init <- c(0.5,0.4,0.1)
  alpha.0 <- alpha.0.init
  beta.0 <- beta.0.init
  
  gamma <- 0.4
  
  Z.g <- vector('list', G)
  for(g in 1:G){
    if(.p.g[g] == 1){
      Z.g[[g]] <- 1
    }
    if(.p.g[g] > 1){
      Z.g[[g]] <- rep(1/.p.g[g], .p.g[g])
    }
  }
  
  theta.g <- Z.g
  V.g <- rep(0, G)
  
  alpha.prev <- alpha
  beta.prev <- beta
  for(em.iter in 1:100){
    # E-step
    for(g in 1:G){
      if(.p.g[g] > 0){
        for(k in 1:.p.g[g]){
          Z.g[[g]][k] <- gamma*theta.g[[g]][k]*prod(c(alpha,beta)^c(n0[[g]][k],n1[[g]][k],n2[[g]][k],m0[[g]][k],m1[[g]][k],m2[[g]][k]))
        }
        
        # null probability
        prob.null <- (1-gamma)*prod(c(alpha.0,beta.0)^c(N00[g],N01[g],N02[g],M00[g],M01[g],M02[g]))
        Z.g[[g]] <- Z.g[[g]]/(sum(Z.g[[g]]) + prob.null)
      }
    }
    
    # M-step
    for(g in 1:G){
      if(.p.g[g] > 0){
        V.g[g] <- sum(Z.g[[g]])
        theta.g[[g]] <- .PI[[g]] + Z.g[[g]]*V.g[g]
        theta.g[[g]] <- theta.g[[g]]/sum(theta.g[[g]])
        N0[g] <- sum(Z.g[[g]]*n0[[g]]) + lambda_0*.p.g[g]
        N1[g] <- sum(Z.g[[g]]*n1[[g]]) + lambda_1*.p.g[g]
        N2[g] <- sum(Z.g[[g]]*n2[[g]]) + lambda_2*.p.g[g]
        M0[g] <- sum(Z.g[[g]]*m0[[g]]) + lambda_0*.p.g[g]
        M1[g] <- sum(Z.g[[g]]*m1[[g]]) + lambda_1*.p.g[g]
        M2[g] <- sum(Z.g[[g]]*m2[[g]]) + lambda_2*.p.g[g]
      }
    }
    alpha <- c(sum(N0*V.g), sum(N1*V.g), sum(N2*V.g))
    beta <- c(sum(M0*V.g), sum(M1*V.g), sum(M2*V.g))
    alpha <- alpha/sum(alpha)
    beta <- beta/sum(beta)
    alpha.0 <- c(sum(N00*(1-V.g)), sum(N01*(1-V.g)), sum(N02*(1-V.g)))
    beta.0 <- c(sum(M00*(1-V.g)), sum(M01*(1-V.g)), sum(M02*(1-V.g)))
    alpha.0 <- alpha.0/sum(alpha.0)
    beta.0 <- beta.0/sum(beta.0)
    
    gamma <- mean(V.g)
    
    # print(em.iter)
    # print(sprintf('alpha = %.2f, %.2f, %.2f', alpha[1], alpha[2], alpha[3]))
    # print(sprintf('beta = %.2f, %.2f, %.2f', beta[1], beta[2], beta[3]))
    # print(sprintf('alpha.0 = %.2f, %.2f, %.2f', alpha.0[1], alpha.0[2], alpha.0[3]))
    # print(sprintf('beta.0 = %.2f, %.2f, %.2f', beta.0[1], beta.0[2], beta.0[3]))
    # print(sprintf('gamma = %.2f', gamma))
    
    if(sum((alpha.prev - alpha)^2) + sum((beta.prev - beta)^2) < 1e-6){
      break
    }
    
    alpha.prev <- alpha
    beta.prev <- beta
  }
  
  params <- list(alpha = alpha, beta = beta,
                 alpha.0 = alpha.0, beta.0 = beta.0,
                 gamma = gamma, em.iter = em.iter,
                 V.g = V.g, Z.g = Z.g)
  return(params)
}
