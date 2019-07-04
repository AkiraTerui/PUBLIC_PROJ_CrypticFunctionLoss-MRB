model{
  ninfo <- 0.001

  # Prior ----
    ##for regression coefficients
    for(k in 1:2){ mean.raw.b[k] ~ dnorm(0, ninfo) }
    b[3] ~ dnorm(0,ninfo)
    ##for variances
    for(k in 1){ tau[k] <- pow(sigma[k], -2); sigma[k] ~ dnorm(0,ninfo)T(0, 100) }
    ##for random effect
    TAU[1:2,1:2] ~ dwish(W[,], 3)
    VAR[1:2,1:2] <- inverse(TAU[,])
    for(k in 1:2){ sigma.R[k] <- sqrt(VAR[k,k]) }
    
  # Likelihood ----
    for(i in 1:Nsample){
      ##variable transformation
      log.TSS[i] <- log(TSS[i])
  
      ##model formula
      log.Y[i] ~ dnorm(mu[i], tau[1])
      mu[i] <- raw.b[Species[i],1] + raw.b[Species[i],2]*Temp[i] + b[3]*log.TSS[i]
      
      ##log-likelihood
      loglik[i] <- logdensity.norm(log.Y[i], mu[i], tau[1])
    }
    ##random effect
    for(j in 1:Nsp){
      raw.b[j,1:2] ~ dmnorm(mean.raw.b[], TAU[,])
    }
    
    ##parameter transformation
    b[1] <- mean.raw.b[1]
    b[2] <- mean.raw.b[2]
    b1_20 <- b[1] + mean.raw.b[2]*20 #add an effect of water temperature at 20C
    
}