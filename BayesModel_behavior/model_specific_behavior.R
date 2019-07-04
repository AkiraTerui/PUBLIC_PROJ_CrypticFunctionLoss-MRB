model{
  ninfo <- 0.001
  
  # Prior ----
    ##for regression coefficients
    for(k in 1:3){ mean.raw.b[k] ~ dnorm(0, ninfo) }
    ##for variances
    for(k in 1){ tau[k] <- pow(sigma[k], -2); sigma[k] ~ dnorm(0,ninfo)T(0, 100) }
    ##for random effect
    TAU[1:3,1:3] ~ dwish(W[,], 4)
    VAR[1:3,1:3] <- inverse(TAU[,])
    for(k in 1:3){ sigma.R[k] <- sqrt(VAR[k,k]) }
    
  # Likelihood ----
    for(i in 1:Nsample){
      ##variable transformation
      log.TSS[i] <- log(TSS[i])
  
      ##model formula
      log.Y[i] ~ dnorm(mu[i], tau[1])
      mu[i] <- raw.b[Species[i],1] + raw.b[Species[i],2]*Temp[i] + raw.b[Species[i],3]*log.TSS[i]
      
      ##log-liklihood
      loglik[i] <- logdensity.norm(log.Y[i], mu[i], tau[1])
    }
    ##random effect
    for(j in 1:Nsp){
      raw.b[j,1:3] ~ dmnorm(mean.raw.b[], TAU[,])
    }
    
    ##parameter transformation
    for(k in 1:3){ b[k] <- mean.raw.b[k] }
    
}