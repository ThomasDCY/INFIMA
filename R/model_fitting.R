#' @name model_fitting
#' @title INFIMA model fitting
#' 
#' @description The main function for INFIMA model fitting.
#' Taking inputs from founder ATAC-seq, founder RNA-seq, DO-eQTL
#' as well as the footprinting annotation data, INFIMA outputs SNP-level
#' quantifications to pinpoint the causal variants.
#' 
#' @param model_data The output from \code{\link{model_input_data}}.
#' @param prior The output from \code{\link{compute_prior}}.
#' @param alpha_init The initial parameters of alpha:
#' The three-dimensional multinomial parameters for the generative model 
#' of normal strains \code{129, AJ, B6, NOD, NZO}.
#' Default \code{alpha_init = c(0.6, 0.3, 0.1)}.
#' @param beta_init The initial parameters of beta:
#' The three-dimensional multinomial parameters for the generative model 
#' of wild strains \code{Cast, PWK, WSB}.
#' Default \code{alpha_init = c(0.5, 0.4, 0.1)}.
#' @param gamma_init The initial value of gamma: 
#' the grand probability of a DO gene to be causal.
#' Default \code{gamma_init = 0.5}.
#' @param max_iter The maximum iteration allowed. 
#' Default \code{max_iter = 100}.
#' @param verbose Print messages or not. 
#' Default \code{verbose = TRUE}.
#' @return An infima object containing the results after model fitting.
#' 
#' Parameters (see the paper for details):
#' \tabular{ll}{
#' \code{alpha} \tab The estimated \code{a1}. \cr
#' \code{beta} \tab The estimated \code{b1}. \cr
#' \code{alpha.0} \tab The estimated \code{a0}. \cr
#' \code{beta.0} \tab The estimated \code{b0}. \cr
#' \code{gamma} \tab The estimated \code{gamma}. \cr
#' }
#' The results of INFIMA model:
#' \tabular{ll}{
#' \code{em.iter} \tab The number of iterations used. \cr
#' \code{V.g} \tab The posterior probability for each DO gene to be causal. \cr
#' \code{Z.g} \tab The posterior probability of candidate SNPs for each DO gene. \cr
#' }
#' 
#' @examples
#' data('example-10-genes')
#' raw_data <- raw_input_data(dt1, dt2, dt3)
#' model_data <- model_input_data(raw_data)
#' prior <- compute_prior(raw_data, model_data)
#' infima <- model_fitting(model_data, prior)
#' 
#' @seealso \code{\link{raw_input_data}}, \code{\link{model_input_data}}, \code{\link{compute_prior}}.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @rawNamespace import(data.table, except = shift)
#' @export
model_fitting <- function(model_data = NULL,
                          prior = NULL,
                          alpha_init = c(0.6, 0.3, 0.1),
                          beta_init = c(0.5, 0.4, 0.1),
                          gamma_init = 0.5,
                          max_iter = 100,
                          verbose = TRUE) {
  if (is.null(model_data)) {
    stop('Error: please input the model_data.')
  }
  
  if (class(model_data) != 'model_data') {
    stop('Error: please input the S3 class model_data.')
  }
  
  if (is.null(prior)) {
    stop('Error: please input the prior.')
  }
  
  if (class(prior) != 'prior') {
    stop('Error: please input the S3 class prior.')
  }
  
  D.g <- model_data$D.g
  p.g <- model_data$p.g
  B.g <- model_data$B.g
  Y.g <- model_data$Y.g
  PI <- prior$PI
  
  
  stopifnot(exprs = {
    !is.null(D.g)
    ! is.null(p.g)
    ! is.null(B.g)
    ! is.null(Y.g)
    ! is.null(PI)
  })
  
  # cannot contain any probability equal to zero
  stopifnot(exprs = {
    sum(alpha_init) == 1
    all(alpha_init > 0)
    sum(beta_init) == 1
    all(beta_init > 0)
    (gamma_init > 0) && (gamma_init < 1)
  })
  
  G <- length(D.g)
  
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
  N0 <- rep(0, G)
  N1 <- rep(0, G)
  N2 <- rep(0, G)
  M0 <- rep(0, G)
  M1 <- rep(0, G)
  M2 <- rep(0, G)
  
  for (g in 1:G) {
    if (p.g[g] > 0) {
      D1.g <- D.g[[g]][, c(1, 2, 3, 5, 6)]
      D2.g <- D.g[[g]][, c(4, 7, 8)]
      if (p.g[g] == 1) {
        D1.g <- t(as.matrix(D1.g))
        D2.g <- t(as.matrix(D2.g))
      }
      n0[[g]] <- apply(D1.g, 1, function(x) {
        sum(x == 0)
      })
      n1[[g]] <- apply(D1.g, 1, function(x) {
        sum(x == 1)
      })
      n2[[g]] <- apply(D1.g, 1, function(x) {
        sum(x == 2)
      })
      m0[[g]] <- apply(D2.g, 1, function(x) {
        sum(x == 0)
      })
      m1[[g]] <- apply(D2.g, 1, function(x) {
        sum(x == 1)
      })
      m2[[g]] <- apply(D2.g, 1, function(x) {
        sum(x == 2)
      })
    }
  }
  
  
  
  ####### The counts for the null generative model
  # alpha.0 <=> a0, beta.0 <=> b0
  N00 <- rep(0, G)
  N01 <- rep(0, G)
  N02 <- rep(0, G)
  M00 <- rep(0, G)
  M01 <- rep(0, G)
  M02 <- rep(0, G)
  
  for (g in 1:G) {
    if (p.g[g] > 0) {
      tmp <- abs(B.g[[g]] - Y.g[[g]])
      tmp1 <- tmp[c(1, 2, 3, 5, 6)]
      tmp2 <- tmp[c(4, 7, 8)]
      
      N00[g] <- sum(tmp1 == 0)
      N01[g] <- sum(tmp1 == 1)
      N02[g] <- sum(tmp1 == 2)
      M00[g] <- sum(tmp2 == 0)
      M01[g] <- sum(tmp2 == 1)
      M02[g] <- sum(tmp2 == 2)
    }
  }
  
  
  ##### the reward terms for enforcing concordance
  ##### we do not set this as tuning parameter
  lambda_0 <- 0.1
  lambda_1 <- 0.01
  lambda_2 <- 0
  
  ### initial values
  alpha <- alpha_init
  beta <- beta_init
  alpha.0 <- alpha_init
  beta.0 <- beta_init
  
  gamma <- gamma_init
  
  Z.g <- vector('list', G)
  for (g in 1:G) {
    if (p.g[g] == 1) {
      Z.g[[g]] <- 1
    }
    if (p.g[g] > 1) {
      Z.g[[g]] <- rep(1 / p.g[g], p.g[g])
    }
  }
  
  theta.g <- Z.g
  V.g <- rep(0, G)
  
  alpha.prev <- alpha
  beta.prev <- beta
  for (em.iter in 1:max_iter) {
    # E-step
    for (g in 1:G) {
      if (p.g[g] > 0) {
        for (k in 1:p.g[g]) {
          Z.g[[g]][k] <-
            gamma * theta.g[[g]][k] * prod(c(alpha, beta) ^ c(n0[[g]][k], n1[[g]][k], n2[[g]][k], m0[[g]][k], m1[[g]][k], m2[[g]][k]))
        }
        
        # null probability
        prob.null <-
          (1 - gamma) * prod(c(alpha.0, beta.0) ^ c(N00[g], N01[g], N02[g], M00[g], M01[g], M02[g]))
        Z.g[[g]] <- Z.g[[g]] / (sum(Z.g[[g]]) + prob.null)
      }
    }
    
    # M-step
    for (g in 1:G) {
      if (p.g[g] > 0) {
        V.g[g] <- sum(Z.g[[g]])
        theta.g[[g]] <- PI[[g]] + Z.g[[g]] * V.g[g]
        theta.g[[g]] <- theta.g[[g]] / sum(theta.g[[g]])
        N0[g] <- sum(Z.g[[g]] * n0[[g]]) + lambda_0 * p.g[g]
        N1[g] <- sum(Z.g[[g]] * n1[[g]]) + lambda_1 * p.g[g]
        N2[g] <- sum(Z.g[[g]] * n2[[g]]) + lambda_2 * p.g[g]
        M0[g] <- sum(Z.g[[g]] * m0[[g]]) + lambda_0 * p.g[g]
        M1[g] <- sum(Z.g[[g]] * m1[[g]]) + lambda_1 * p.g[g]
        M2[g] <- sum(Z.g[[g]] * m2[[g]]) + lambda_2 * p.g[g]
      }
    }
    alpha <- c(sum(N0 * V.g), sum(N1 * V.g), sum(N2 * V.g))
    beta <- c(sum(M0 * V.g), sum(M1 * V.g), sum(M2 * V.g))
    alpha <- alpha / sum(alpha)
    beta <- beta / sum(beta)
    alpha.0 <-
      c(sum(N00 * (1 - V.g)), sum(N01 * (1 - V.g)), sum(N02 * (1 - V.g)))
    beta.0 <-
      c(sum(M00 * (1 - V.g)), sum(M01 * (1 - V.g)), sum(M02 * (1 - V.g)))
    alpha.0 <- alpha.0 / sum(alpha.0)
    beta.0 <- beta.0 / sum(beta.0)
    
    gamma <- mean(V.g)
    
    if (verbose) {
      print(em.iter)
      print(sprintf('alpha = %.2f, %.2f, %.2f', alpha[1], alpha[2], alpha[3]))
      print(sprintf('beta = %.2f, %.2f, %.2f', beta[1], beta[2], beta[3]))
      print(sprintf('alpha.0 = %.2f, %.2f, %.2f', alpha.0[1], alpha.0[2], alpha.0[3]))
      print(sprintf('beta.0 = %.2f, %.2f, %.2f', beta.0[1], beta.0[2], beta.0[3]))
      print(sprintf('gamma = %.2f', gamma))
    }
    
    if (sum((alpha.prev - alpha) ^ 2) + sum((beta.prev - beta) ^ 2) < 1e-6) {
      break
    }
    
    alpha.prev <- alpha
    beta.prev <- beta
  }
  
  infima <- list(
    alpha = alpha,
    beta = beta,
    alpha.0 = alpha.0,
    beta.0 = beta.0,
    gamma = gamma,
    em.iter = em.iter,
    V.g = V.g,
    Z.g = Z.g
  )
  
  class(infima) <- 'infima'
  return(infima)
}
