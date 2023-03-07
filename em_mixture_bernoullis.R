logsumexp <- function(x) return(log(sum(exp(x - max(x)))) + max(x))
## define prob
prob <- function(x,mu,return.log=FALSE) {
  l <- sum(log(mu[x==1]))+sum(log(1-mu[x==0]))
  if (return.log) {
    return(l)
  } else {
    return(exp(l))
  }
}

##compute loglikelihood
Loglikelihood <- function(xs,mus,lws,gammas) {
  ll <- 0
  n <- dim(xs)[1]
  K <- dim(mus)[1]
  for (i in 1:n) {
    for (k in 1:K) {
      if (gammas[i,k] > 0) {
        ll <- ll + gammas[i,k]*(lws[k]+prob(xs[i,],mus[k,],return.log=TRUE)-log(gammas[i,k]))
      }
    }
  }
  return(ll)
}


em_mix_bernoulli <- function(xs,K,start=NULL,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws), ws is mixture weights from lab4
  lws <- rep(log(1/K),K)
  if (is.null(start)) {
    mus <- .2 + .6*xs[sample(n,K),]
  } else {
    mus <- start
  }
  gammas <- matrix(0,n,K)
  
  converged <- FALSE
  numit <- 0
  ll <- -Inf
  while(!converged && numit < max.numit) {
    numit <- numit + 1
    mus.old <- mus
    ll.old <- ll
    
    ## E step - calculate gammas
    for (i in 1:n) {
      # the elements of lprs are log(w_k * p_k(x)) for each k in {1,...K}
      lprs <- rep(0,K)
      for (k in 1:K) {
        lprs[k] <- lws[k] + prob(xs[i,],mus[k,],return.log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- Log-likelihood(xs,mus,lws,gammas)
    
    # M step - update ws and mus
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
      }
    }
    # from lab4
    mus[which(mus > 1,arr.ind=TRUE)] <- 1 - 1e-15

    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,gammas=gammas,ll=ll, numit = numit))
}
