model{
  # Prior ----
    ninfo <- 1.0E-3
    U <- 30
    
    for(j in 1:Nsp){ b[j,1:5] ~ dmnorm(B[], TAU[,]) }
    for(m in 1:5){ B[m] ~ dnorm(0,ninfo) }
    
    TAU[1:5,1:5] ~ dwish(W[,], 6)
    SIGMA[1:5,1:5] <- inverse(TAU[,])
    
    for(t in 1:Nyear){ R[t] ~ dnorm(0,tau.R[2]) }
    for(i in 1:2){ tau.R[i] <- pow(sigma.R[i],-2); sigma.R[i] ~ dnorm(0,ninfo)T(0,U) }
    tau.eps <- pow(sigma.eps,-2); sigma.eps ~ dnorm(0,ninfo)T(0,U)
  
  # Likelihood----
    ## lambda.semi: Poisson mean for uncorrected timed-search data
    ## log.lambda.obs: Density with observation error (effort corrected)
    ## log.lambda: Density without observation error (effort corrected)
    
    for(n in 1:Nsample){
      Lik[n] <- log( dpois(Y[n], lambda.semi[ SITE[n], SP[n] ]) )
      Y[n] ~ dpois(lambda.semi[ SITE[n], SP[n] ])
    }
    
    for(i in 1:Nsite){
      scl.ELEV[i] <- (ELEV[i] - mean(ELEV[]))/sd(ELEV[])
      for(j in 1:Nsp){
        log(lambda.semi[i,j]) <- log.lambda.obs[i,j] - log(0.0712) + log(EFFORT[i]/60)
        log.lambda.obs[i,j] ~ dnorm(log.lambda[i,j], tau.eps)
        log.lambda[i,j] <- x[HUC12[i],j] + b[j,5]*scl.ELEV[i] + R[ YEAR[i] ]
      }
    }
    
    
    for(k in 1:K){
      scl.TSS[k] <- ( TSS[k] - mean(TSS[]) )/sd(TSS[])
      scl.TN[k] <- ( TN[k] - mean(TN[]) )/sd(TN[])
      scl.TP[k] <- ( TP[k] - mean(TP[]) )/sd(TP[])
      for(j in 1:Nsp){
        x[k,j] ~ dnorm(x_mu[k,j], tau.R[1])
        x_mu[k,j] <- b[j,1] + b[j,2]*scl.TSS[k] + b[j,3]*scl.TN[k] + b[j,4]*scl.TP[k]
      }
    }
    
    for(j in 1:Nsp){
      ## unscaled parameters
      Alpha[j] <- b[j,1] - b[j,2]*( mean(TSS[])/sd(TSS[]) )
      Beta[j] <- b[j,2]/sd(TSS[])
      
      ## scaled (standardized) parameters
      alpha[j] <- b[j,1]
      betaTSS[j] <- b[j,2]
      betaTN[j] <- b[j,3]
      betaTP[j] <- b[j,4]
      betaELEV[j] <- b[j,5]
    }
    
  # Bayesian p-value ----
    for(n in 1:Nsample){
      # ideal "replicated" data
      Y.new[n] ~ dpois(lambda.semi[ SITE[n], SP[n] ])
      # predicted values
      y.pred[n] <- lambda.semi[ SITE[n], SP[n] ]
    
      # residuals from data
      resid1[n] <- Y[n] - y.pred[n]
      # residuals from ideal data
      resid2[n] <- Y.new[n] - y.pred[n]
      
      # Discrepancy measure
      sq1[n] <- pow(resid1[n], 2)
      sq2[n] <- pow(resid2[n], 2)
    }
    fit <- sum(sq1[])
    fit.new <- sum(sq2[])
    test <- step(fit.new - fit)
    
}