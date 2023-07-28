# Functions in SPSL simulation for exponential data

fisher = function(x)
  return(1 / 2 * log((1 + x) / (1 - x)))
inv.fisher = function(x)
  return((exp(2 * x) - 1) / (exp(2 * x) + 1))

get_data = function(n, nb, tau0, tau1, z){
  q = choose(nb, 2)
  
  #z = runif(n)
  rho = inv.fisher(t(t(outer(z, tau1)) + tau0))
  
  x = matrix(rep(0, nb * n), nrow = n)
  
  for (i in 1:n) {
    cor.mat = diag(rep(0.5, nb))
    cor.mat[lower.tri(cor.mat)] = rho[i, ]
    cor.mat = cor.mat + t(cor.mat)
    while (is.positive.definite(cor.mat) == F) {
      z[i] = z[i] + runif(1,-.0001,.0001)
      rho[i, ] = inv.fisher(t(t(outer(z[i], tau1)) + tau0))
      cor.mat = diag(rep(0.5, nb))
      cor.mat[lower.tri(cor.mat)] = rho[i, ]
      cor.mat = cor.mat + t(cor.mat)
    }
    x[i, ] = rnorm(nb) %*% chol(cor.mat) # Cholesky decomposition
  }
  x.scale = apply(x,2,scale)
  
  x = qexp(pnorm(x.scale), rate = 0.65)
  
  x = apply(x, 2, function(y,n){qnorm(rank(y)/(n+1))},n = n)
  x = apply(x, 2, scale)
  
  return(list(x, z))
}

################################################################################

# SPSL 

spsl =
  function(data,
           z,
           startvalue = rep(c(0, 0, 0, 0), choose(ncol(data), 2)),
           t0 = 200,
           v0 = 0.005,
           a1 = 5,
           a2 = 50,
           iterations = 10000,
           burnIn = 5000,
           chains = 3,
           e = 0.05,
           Sd = 10) {
    
    # data is the matrix of gene expression levels, nrows are sample size
    # startvalue is a vector of starting values of all parameters
    # t0 is the number of iterations in MCMC history to calculate covariance matrix for the implementation of adaptive MCMC
    # a1 and a2 are the shape and scale parameters in bimodal inverse gamma distributions for variable selection in SPSL and C-SPSL
    # e and Sd are step sizesin regular MCMC and adaptive MCMC
    
    nb = ncol(data) # number of gene
    n = nrow(data) # sample size
    q = choose(nb, 2) # size of combination
    getY <- function(x, nb) {
      # getY is the function to create Y=X1*X2
      y = c()
      col = nb
      for (j in 1:(col - 1)) {
        for (k in (j + 1):col) {
          y = cbind(y, x[, j] * x[, k])
        }
      }
      return(as.matrix(y))
    }
    y <- getY(scale(data), nb = nb)
    z = as.matrix(z)
    
    col = length(startvalue)
    tau0wh = seq(1, col, by = 4)
    tau1wh = seq(2, col, by = 4)
    v2wh = seq(3, col, by = 4)
    wwh = seq(4, col, by = 4)
    fisher <- function(x){
      return(1 / 2 * log((1 + x) / (1 - x)))
    }
    inv.fisher <- function(x){
      return((exp(2 * x) - 1) / (exp(2 * x) + 1))
    }
    logit = function(x) {
      return(log(x / (1 - x)))
    }
    inv.logit = function(x) {
      return(exp(x) / (1 + exp(x)))
    }
    dinv.logit = function(x) {
      return(exp(x) / (1 + exp(x)) ^ 2)
    }
    singlelikelihoods = function(rho, y) {
      return(-log(pi) - .5 * log(1 - rho ^ 2) + rho * y / (1 - rho ^ 2) + log(besselK(abs(y) /
                                                                                        (1 - rho ^ 2), nu = 0)))
    }
    f.invgamma = function(x, a, b, v0) {
      return(b ^ a / gamma(a) * (1 / x) ^ (a + 1) * exp(-b * v0 / x) * v0 ^ a)
    }
    likelihood <- function(tau0, tau1, z, y) {
      rho = inv.fisher(t(tau0 + t(t(outer(
        tau1, z
      )))))
      sumll = sum(unlist(mapply(singlelikelihoods, rho, y)))
      return(sumll)
    }
    prior.tau0 = function(tau0) {
      return(dnorm(
        tau0,
        mean = 0,
        sd = 1,
        log = T
      ))
    }
    
    update.tau0 <- function(param, z, y, ct) {
      tau1 = param[tau1wh]
      old.tau0 = param[tau0wh]
      proposal = c(rnorm(length(tau0wh)) %*% chol(ct)) + old.tau0
      new.tau0 = old.tau0
      for (j in 1:length(tau0wh)) {
        pro.rho = inv.fisher(proposal[j] + tau1[j] * z)
        old.rho = inv.fisher(old.tau0[j] + tau1[j] * z)
        post.ratio = exp(
          prior.tau0(proposal[j]) + sum(singlelikelihoods(pro.rho, y[, j])) - prior.tau0(old.tau0[j]) -
            sum(singlelikelihoods(old.rho, y[, j]))
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept),FALSE,accept)
        if (accept) {
          new.tau0[j] = proposal[j]
        } else {
          new.tau0[j] = old.tau0[j]
        }
      }
      return(new.tau0)
    }
    
    prior.tau1 = function(tau1, v2) {
      return((dnorm(
        tau1,
        mean = 0,
        sd = sqrt(exp(v2)),
        log = T
      )))
    }
    update.tau1 <- function(param, z, y, ct) {
      tau0 = param[tau0wh]
      v2 = param[v2wh]
      old.tau1 = param[tau1wh]
      proposal = c(rnorm(length(tau1wh)) %*% chol(ct)) + old.tau1
      new.tau1 = old.tau1
      for (j in 1:length(tau1wh)) {
        pro.rho = inv.fisher(tau0[j] + proposal[j] * z)
        old.rho = inv.fisher(tau0[j] + old.tau1[j] * z)
        post.ratio = exp(
          prior.tau1(proposal[j], v2[j]) + sum(singlelikelihoods(pro.rho, y[, j])) -
            prior.tau1(old.tau1[j], v2[j]) - sum(singlelikelihoods(old.rho, y[, j]))
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept),FALSE,accept)
        if (accept) {
          new.tau1[j] = proposal[j]
        } else {
          new.tau1[j] = old.tau1[j]
        }
      }
      return(new.tau1)
    }
    
    prior.v2 = function(v2, w, v0, a1, a2) {
      return((log(((1 - inv.logit(w)) * f.invgamma(exp(v2), a1, a2, v0) + inv.logit(w) *
                     dinvgamma(exp(v2), a1, a2)
      ) * exp(v2)
      )))
      
    }
    
    update.v2 <- function(param, z, y, ct) {
      old.v2 = param[v2wh]
      tau1 = param[tau1wh]
      w = param[wwh]
      proposal = c(rnorm(length(v2wh)) %*% chol(ct)) + old.v2
      new.v2 = old.v2
      for (j in 1:length(v2wh)) {
        post.ratio = exp(
          prior.v2(proposal[j], w[j], v0, a1, a2) + prior.tau1(tau1[j], proposal[j])
          - prior.v2(old.v2[j], w[j], v0, a1, a2) - prior.tau1(tau1[j], old.v2[j])
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept),FALSE,accept)
        if (accept) {
          new.v2[j] = proposal[j]
        } else {
          new.v2[j] = old.v2[j]
        }
      }
      return(new.v2)
    }
    
    prior.w = function(w) {
      return((log(
        dunif(inv.logit(w), 0, 1) * dinv.logit(w)
      )))
    }
    
    update.w <- function(param, z, y, ct) {
      old.w = param[wwh]
      v2 = param[v2wh]
      proposal = c(rnorm(length(wwh)) %*% chol(ct)) + old.w
      new.w = old.w
      for (j in 1:length(wwh)) {
        post.ratio = exp(
          prior.w(proposal[j]) + prior.v2(v2[j], proposal[j], v0, a1, a2)
          - prior.w(old.w[j]) - prior.v2(v2[j], old.w[j], v0, a1, a2)
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept),FALSE,accept)
        if (accept) {
          new.w[j] = proposal[j]
        } else {
          new.w[j] = old.w[j]
        }
      }
      return(new.w)
    }
    
    
    MH_MCMC <- function(startvalue, iterations, t0, z, y) {
      
      # e is the step size for regular MCMC
      # Sd is the step size for adaptive MCMC
      
      chain = array(dim = c(iterations + 1, col))
      chain[1, ] = startvalue
      c0 = diag(length(tau0wh)) * e  # covariance matrix for parameters
      for (i in 1:iterations) {
        if (i <= t0 | length(tau0wh)==1) {
          # When iteration is less than t0, use regular MCMC to update
          param = chain[i, ]
          param[wwh] = update.w(chain[i, ], z, y, c0)
          param[v2wh] = update.v2(chain[i, ], z, y, c0)
          param[tau1wh] = update.tau1(chain[i, ], z, y, c0)
          param[tau0wh] = update.tau0(chain[i, ], z, y, c0)
          chain[i + 1, ] = param
        } else if (i > t0 & length(tau0wh)>1) {
          # When iteration is greater than t0, use adaptive MCMC
          
          c_tau0 = (cov(chain[(i - t0 + 1):i, tau0wh]) + diag(tau0wh) * 0.00001) * Sd  # covariance matrix for tau0
          c_tau1 = (cov(chain[(i - t0 + 1):i, tau1wh]) + diag(tau0wh) * 0.00001) * Sd  # covariance matrix for tau1
          c_w = (cov(chain[(i - t0 + 1):i, wwh]) + diag(tau0wh) * 0.00001) * Sd  # covariance matrix for w
          c_v2 = (cov(chain[(i - t0 + 1):i, v2wh]) + diag(tau0wh) * 0.00001) * Sd  # covariance matrix for v square
          param = chain[i, ]
          param[tau0wh] = update.tau0(chain[i, ], z, y, c_tau0)
          param[wwh] = update.w(chain[i, ], z, y, c_w)
          param[v2wh] = update.v2(chain[i, ], z, y, c_v2)
          param[tau1wh] = update.tau1(chain[i, ], z, y, c_tau1)
          chain[i + 1, ] = param
        }
        #cat(i/iterations*100,"%","\n",sep="")
      }
      return(chain)
    }
    chain = list()
    for (i in 1:chains) {
      #cat("Chain-",i,"\n",sep="")
      result = MH_MCMC(startvalue, iterations, t0, z, y)[-(1:(burnIn + 1)), ]
      tau0 = result[, tau0wh]
      tau1 = result[, tau1wh]
      phi = exp(result[, v2wh])
      w = inv.logit(result[, wwh])
      resultlist = list(tau0, tau1, phi, w)
      names(resultlist) = c("tau0", "tau1", "phi", "w")
      chain[[i]] = resultlist
      #cat("SPSL acceptence rate:", mean(apply(result, 2, function(x)
      #  1 - mean(duplicated(x)))), "\n")
    }
    return(chain)
  }

################################################################################

# C-SPSL

cspsl =
  function(data,
           z,
           startvalue = c(rep(c(0, 0, 0, 0), choose(ncol(data), 2)),0),
           t0 = 200,
           v0 = 0.005,
           a1 = 5,
           a2 = 50,
           iterations = 10000,
           burnIn = 5000,
           chains = 3,
           e = 0.05,
           Sd = 10) {
    nb = ncol(data) # number of gene
    n = nrow(data)
    q = choose(nb, 2) # size of combination
    getY <- function(x, nb) {
      y = c()
      col = nb
      for (j in 1:(col - 1)) {
        for (k in (j + 1):col) {
          y = cbind(y, x[, j] * x[, k])
        }
      }
      return(as.matrix(y))
    }
    y <- getY(scale(x), nb = nb)
    z = as.matrix(z)
    
    col = length(startvalue)
    tau0wh = seq(1, col - 1, by = 4)
    tau1wh = seq(2, col - 1, by = 4)
    v2wh = seq(3, col - 1, by = 4)
    wwh = seq(4, col - 1, by = 4)
    fisher <- function(x){
      return(1 / 2 * log((1 + x) / (1 - x)))
    }
    inv.fisher <- function(x){
      return((exp(2 * x) - 1) / (exp(2 * x) + 1))
    }
    logit = function(x) {
      return(log(x / (1 - x)))
    }
    inv.logit = function(x) {
      return(exp(x) / (1 + exp(x)))
    }
    dinv.logit = function(x) {
      return(exp(x) / (1 + exp(x)) ^ 2)
    }
    singlelikelihoods = function(rho, y) {
      return(-log(pi) - .5 * log(1 - rho ^ 2) + rho * y / (1 - rho ^ 2) + log(besselK(abs(y) /
                                                                                        (1 - rho ^ 2), nu = 0)))
    }
    f.invgamma = function(x, a, b, v0) {
      return(b ^ a / gamma(a) * (1 / x) ^ (a + 1) * exp(-b * v0 / x) * v0 ^ a)
    }
    likelihood <- function(tau0, tau1, z, y) {
      rho = inv.fisher(t(tau0 + t(t(outer(
        tau1, z
      )))))
      sumll = sum(unlist(mapply(singlelikelihoods, rho, y)))
      return(sumll)
    }
    prior.tau0 = function(tau0) {
      return(dnorm(
        tau0,
        mean = 0,
        sd = 1,
        log = T
      ))
    }
    
    update.tau0 <- function(param, z, y, ct) {
      tau1 = param[tau1wh]
      old.tau0 = param[tau0wh]
      proposal = c(rnorm(length(tau0wh)) %*% chol(ct)) + old.tau0
      new.tau0 = old.tau0
      for (j in 1:length(tau0wh)) {
        pro.rho = inv.fisher(proposal[j] + tau1[j] * z)
        old.rho = inv.fisher(old.tau0[j] + tau1[j] * z)
        post.ratio = exp(
          prior.tau0(proposal[j]) + sum(singlelikelihoods(pro.rho, y[, j])) - prior.tau0(old.tau0[j]) -
            sum(singlelikelihoods(old.rho, y[, j]))
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept), FALSE, accept)
        if (accept) {
          new.tau0[j] = proposal[j]
        } else {
          new.tau0[j] = old.tau0[j]
        }
      }
      return(new.tau0)
    }
    
    prior.tau1 = function(tau1, v2) {
      return((dnorm(
        tau1,
        mean = 0,
        sd = sqrt(exp(v2)),
        log = T
      )))
    }
    update.tau1 <- function(param, z, y, ct) {
      tau0 = param[tau0wh]
      v2 = param[v2wh]
      old.tau1 = param[tau1wh]
      proposal = c(rnorm(length(tau1wh)) %*% chol(ct)) + old.tau1
      new.tau1 = old.tau1
      for (j in 1:length(tau1wh)) {
        pro.rho = inv.fisher(tau0[j] + proposal[j] * z)
        old.rho = inv.fisher(tau0[j] + old.tau1[j] * z)
        post.ratio = exp(
          prior.tau1(proposal[j], v2[j]) + sum(singlelikelihoods(pro.rho, y[, j])) -
            prior.tau1(old.tau1[j], v2[j]) - sum(singlelikelihoods(old.rho, y[, j]))
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept), FALSE, accept)
        if (accept) {
          new.tau1[j] = proposal[j]
        } else {
          new.tau1[j] = old.tau1[j]
        }
      }
      return(new.tau1)
    }
    
    prior.v2 = function(v2, w, v0, a1, a2) {
      return((log(((1 - inv.logit(w)) * f.invgamma(exp(v2), a1, a2, v0) + inv.logit(w) *
                     dinvgamma(exp(v2), a1, a2)
      ) * exp(v2)
      )))
      
    }
    
    update.v2 <- function(param, z, y, ct) {
      old.v2 = param[v2wh]
      tau1 = param[tau1wh]
      w = param[wwh]
      proposal = c(rnorm(length(v2wh)) %*% chol(ct)) + old.v2
      new.v2 = old.v2
      for (j in 1:length(v2wh)) {
        post.ratio = exp(
          prior.v2(proposal[j], w[j], v0, a1, a2) + prior.tau1(tau1[j], proposal[j])
          - prior.v2(old.v2[j], w[j], v0, a1, a2) - prior.tau1(tau1[j], old.v2[j])
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept), FALSE, accept)
        if (accept) {
          new.v2[j] = proposal[j]
        } else {
          new.v2[j] = old.v2[j]
        }
      }
      return(new.v2)
    }
    
    prior.w = function(w, c) {
      return((log(
        dbeta(
          inv.logit(w) + 0.0001 * (1 - inv.logit(w)) - 0.0001 * inv.logit(w),
          inv.logit(c),
          1 - inv.logit(c)
        ) * dinv.logit(w)
      )))
      
    }
    
    update.w <- function(param, z, y, ct) {
      old.w = param[wwh]
      c = param[col]
      v2 = param[v2wh]
      proposal = c(rnorm(length(wwh)) %*% chol(ct)) + old.w
      new.w = old.w
      for (j in 1:length(wwh)) {
        post.ratio = exp(
          prior.w(proposal[j], c) + prior.v2(v2[j], proposal[j], v0, a1 = 5, a2 =
                                               50)
          - prior.w(old.w[j], c) - prior.v2(v2[j], old.w[j], v0, a1 =
                                              5, a2 = 50)
        )
		accept = (runif(1) < post.ratio)
		accept = ifelse(is.na(accept), FALSE, accept)
        if (accept) {
          new.w[j] = proposal[j]
        } else {
          new.w[j] = old.w[j]
        }
      }
      return(new.w)
    }
    
    
    prior.c = function(c) {
      return((log(
        dunif(inv.logit(c), 0, 1) * dinv.logit(c)
      )))
    }
    
    update.c <- function(param, z, y, e) {
      old.c = param[col]
      w = param[wwh]
      proposal = rnorm(1, mean = 0, sd = sqrt(e)) + old.c
      post.ratio = exp(prior.c(proposal) + sum(prior.w(w, proposal)) - prior.c(old.c) -
                         sum(prior.w(w, old.c)))
	  accept = (runif(1) < post.ratio)
	  accept = ifelse(is.na(accept), FALSE, accept)
      if (accept) {
        return(proposal)
      } else {
        return(old.c)
      }
    }
    
    MH_MCMC <- function(startvalue, iterations, t0, z, y) {
      chain = array(dim = c(iterations + 1, col))
      chain[1, ] = startvalue
      c0 = diag(length(tau0wh)) * e
      for (i in 1:iterations) {
        if (i <= t0 | length(tau0wh)==1) {
          param = chain[i, ]
          param[wwh] = update.w(chain[i, ], z, y, c0)
          param[v2wh] = update.v2(chain[i, ], z, y, c0)
          param[tau1wh] = update.tau1(chain[i, ], z, y, c0)
          param[col] = update.c(chain[i, ], z, y, e)
          param[tau0wh] = update.tau0(chain[i, ], z, y, c0)
          chain[i + 1, ] = param
        } else if (i > t0 & length(tau0wh)>1) {
          c_tau0 = (cov(chain[(i - t0 + 1):i, tau0wh]) + diag(tau0wh) * 0.00001) *
            Sd
          c_tau1 = (cov(chain[(i - t0 + 1):i, tau1wh]) + diag(tau0wh) *
                      0.00001) * Sd
          c_w = (cov(chain[(i - t0 + 1):i, wwh]) + diag(tau0wh) * 0.00001) *
            Sd
          c_v2 = (cov(chain[(i - t0 + 1):i, v2wh]) + diag(tau0wh) * 0.00001) *
            Sd
          param = chain[i, ]
          param[tau0wh] = update.tau0(chain[i, ], z, y, c_tau0)
          param[wwh] = update.w(chain[i, ], z, y, c_w)
          param[v2wh] = update.v2(chain[i, ], z, y, c_v2)
          param[tau1wh] = update.tau1(chain[i, ], z, y, c_tau1)
          param[col] = update.c(chain[i, ], z, y, e)
          chain[i + 1, ] = param
        }
        #cat(i/iterations*100,"%","\n",sep="")
      }
      return(chain)
    }
    chain = list()
    for (i in 1:chains) {
      #cat("Chain-",i,"\n",sep="")
      result = MH_MCMC(startvalue, iterations, t0, z, y)[-(1:(burnIn + 1)),]
      tau0 = result[, tau0wh]
      tau1 = result[, tau1wh]
      phi = exp(result[, v2wh])
      w = inv.logit(result[, wwh])
      c = inv.logit(result[, col])
      resultlist = list(tau0, tau1, phi, w, c)
      names(resultlist) = c("tau0", "tau1", "phi", "w", "c")
      chain[[i]] = resultlist
      #cat("SPSL acceptence rate:", mean(apply(result, 2, function(x)
      #  1 - mean(duplicated(x)))), "\n")
    }
    return(chain)
  }

################################################################################

# FDR functions

summary.stats = function(tau1.post, e, alpha=0.05){
  bf = v = NULL
  n = nrow(tau1.post)
  d = ncol(tau1.post)
  
  fdr = function(c){
    sum(as.numeric(bf>c)*(1-v))/(sum(bf>c)+0.0001)
  }
  
  for(j in 1:d){
    v0 = 0.005
    v[j] = mean(abs(tau1.post[,j])>e)
    w = runif(n)
    v2 = NULL
    for (l in 1:n){
      u = runif(1)
      if (u < 1-w[l]){
        v2[l] = v0*rinvgamma(1,5,50)
      } else if (u < 1) {
        v2[l] = rinvgamma(1,5,50)
      }
    }
    tau1.prior = rnorm(n, 0, sd=sqrt(v2))
    bf[j] = v[j]/(1-v[j]) # /
      #(mean(abs(tau1.prior)>e)/(1-mean(abs(tau1.prior)>e)))
  }
  cv = seq(0,50,by=0.01)
  fdr.bar = NULL
  i = 0
  for (c in cv){
    i = i + 1
    fdr.bar[i] = fdr(c)
  }
  c = cv[which(fdr.bar<=alpha)[1]]
  return(list(bf = bf, r.opt = c))
}
