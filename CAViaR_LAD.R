TAV_CAViaR <- function(y,Beta0,THETA,flag){
  ysort <- sort(y[1:300])
  empiricalQuantile <- ysort[round(300*THETA)]
  N <- length(y)
  if(flag==1){
    RQobjective <- function(Beta){
      VaR <- rep(NA,N)
      VaR[1] <- empiricalQuantile
      for(i in 2:N){
        VaR[i] <- (Beta[1]+Beta[2]*VaR[i-1]+Beta[3]*abs(y[i-1]))*ifelse(y[i-1]<0,1,0)+
          (Beta[4]+Beta[5]*VaR[i-1]+Beta[6]*abs(y[i-1]))*ifelse(y[i-1]>=0,1,0)
      }
      Hit <- (y<VaR)-THETA
      RQ = sum(-Hit*(y-VaR)/N)
      return(RQ)
    }
    ui <- rbind(c(0,1,0,0,0,0),c(0,-1,0,0,0,0),
                c(0,0,0,0,1,0),c(0,0,0,0,-1,0))
    ci <- c(0,-1,0,-1)
    m <- constrOptim(Beta0, RQobjective, NULL,ui=ui,ci=ci,method = "Nelder-Mead")
    Beta <- m$par
    VaR <- rep(NA,N)
    VaR[1] <- empiricalQuantile
    for(i in 2:N){
      VaR[i] <- (Beta[1]+Beta[2]*VaR[i-1]+Beta[3]*abs(y[i-1]))*ifelse(y[i-1]<0,1,0)+
        (Beta[4]+Beta[5]*VaR[i-1]+Beta[6]*abs(y[i-1]))*ifelse(y[i-1]>=0,1,0)
    }
    return(list(beta=Beta,VaR=VaR))
  }else{
    VaR <- rep(NA,N)
    VaR[1] <- empiricalQuantile
    Beta <- Beta0
    for(i in 2:N){
      VaR[i] <- (Beta[1]+Beta[2]*VaR[i-1]+Beta[3]*abs(y[i-1]))*ifelse(y[i-1]<0,1,0)+
        (Beta[4]+Beta[5]*VaR[i-1]+Beta[6]*abs(y[i-1]))*ifelse(y[i-1]>=0,1,0)
    }
    return(VaR)
  }
}

res <- TAV_CAViaR(y,rep(0.1,6),0.05,flag=1)
beta <- c(-1.6203174,0.6256671,-0.2932310,0.4477667,0.9999999,-0.2258167)
beta <- c(0.02086300,  0.99999997, -0.06658941, -0.06923799,  0.96447306,  0.02539998)

TIG_CAViaR <- function(y,Beta0,THETA,flag){
  ysort <- sort(y[1:300])
  empiricalQuantile <- ysort[round(300*THETA)]
  N <- length(y)
  if(flag==1){
    RQobjective <- function(Beta){
      VaR <- rep(NA,N)
      VaR[1] <- empiricalQuantile
      for(i in 2:N){
        VaR[i] <- -sqrt(Beta[1]+Beta[2]*VaR[i-1]^2+Beta[3]*y[i-1]^2)*ifelse(y[i-1]<0,1,0)-
          sqrt(Beta[4]+Beta[5]*VaR[i-1]^2+Beta[6]*y[i-1]^2)*ifelse(y[i-1]>=0,1,0)
      }
      Hit <- (y<VaR)-THETA
      RQ = sum(-Hit*(y-VaR)/N)
      return(RQ)
    }
    ui <- rbind(c(1,0,0,0,0,0),c(0,0,1,0,0,0),
                c(0,0,0,1,0,0),c(0,0,0,0,0,1),
                c(0,1,0,0,0,0),c(0,-1,0,0,0,0),
                c(0,0,0,0,1,0),c(0,0,0,0,-1,0))
    ci <- c(0,0,0,0,0,-1,0,-1)
    m <- constrOptim(Beta0, RQobjective, NULL,ui=ui,ci=ci,method = "Nelder-Mead")
    Beta <- m$par
    VaR <- rep(NA,N)
    VaR[1] <- empiricalQuantile
    for(i in 2:N){
      VaR[i] <- -sqrt(Beta[1]+Beta[2]*VaR[i-1]^2+Beta[3]*y[i-1]^2)*ifelse(y[i-1]<0,1,0)-
        sqrt(Beta[4]+Beta[5]*VaR[i-1]^2+Beta[6]*y[i-1]^2)*ifelse(y[i-1]>=0,1,0)
    }
    return(list(beta=Beta,VaR=VaR))
  }else{
    VaR <- rep(NA,N)
    VaR[1] <- empiricalQuantile
    Beta <- Beta0
    for(i in 2:N){
      VaR[i] <- -sqrt(Beta[1]+Beta[2]*VaR[i-1]^2+Beta[3]*y[i-1]^2)*ifelse(y[i-1]<0,1,0)-
        sqrt(Beta[4]+Beta[5]*VaR[i-1]^2+Beta[6]*y[i-1]^2)*ifelse(y[i-1]>=0,1,0)
    }
    return(VaR)
  }
}

res <- TIG_CAViaR(y,mu,0.01,flag=2)
beta <- c(6.3147033762, 0.7141486932, 0.7971350072, 0.0002700342,
          0.8661760917, 0.5042520775)
beta <-c(6.607129e-02 ,9.999998e-01, 8.130735e-02, 1.500087e-02,
         9.601042e-01, 8.394343e-08)


TIMP_CAViaR <- function(y,Beta0,THETA,flag){
  ysort <- sort(y[1:300])
  empiricalQuantile <- ysort[round(300*THETA)]
  N <- length(y)
  if(flag==1){
    RQobjective <- function(beta){
      VaR <- rep(NA,N)
      VaR[1] <- empiricalQuantile
      for(i in 2:N){
        if(y[i-1]<0){
          VaR[i] <- beta[1]+beta[2]*VaR[i-1]+(1-beta[2])*sqrt(beta[3]^2+
                   (1-beta[3])^2)/(1-beta[3])*abs(y[i-1]-mean(y))
        }else{
          VaR[i] <- beta[4]+beta[5]*VaR[i-1]+(1-beta[5])*sqrt(beta[6]^2+
                    (1-beta[6])^2)/(beta[6])*abs(y[i-1]-mean(y))
        }
      }
      Hit <- (y<VaR)-THETA
      RQ = sum(-Hit*(y-VaR)/N)
      return(RQ)
    }
    ui <- rbind(c(0,0,1,0,0,0),c(0,0,-1,0,0,0),
                c(0,0,0,0,0,1),c(0,0,0,0,0,-1),
                c(0,1,0,0,0,0),c(0,-1,0,0,0,0),
                c(0,0,0,0,1,0),c(0,0,0,0,-1,0))
    ci <- c(0,-1,0,-1,0,-1,0,-1)
    m <- constrOptim(Beta0, RQobjective, NULL,ui=ui,ci=ci,method = "Nelder-Mead")
    beta <- m$par
    VaR <- rep(NA,N)
    VaR[1] <- empiricalQuantile
    for(i in 2:N){
      if(y[i-1]<0){
        VaR[i] <- beta[1]+beta[2]*VaR[i-1]+(1-beta[2])*sqrt(beta[3]^2+
                                                              (1-beta[3])^2)/(1-beta[3])*abs(y[i-1]-mean(y))
      }else{
        VaR[i] <- beta[4]+beta[5]*VaR[i-1]+(1-beta[5])*sqrt(beta[6]^2+
                                                              (1-beta[6])^2)/(beta[6])*abs(y[i-1]-mean(y))
      }
    }
    return(list(beta=beta,VaR=VaR))
  }else{
    VaR <- rep(NA,N)
    VaR[1] <- empiricalQuantile
    beta <- Beta0
    for(i in 2:N){
      if(y[i-1]<0){
        VaR[i] <- beta[1]+beta[2]*VaR[i-1]+(1-beta[2])*sqrt(beta[3]^2+
                                                              (1-beta[3])^2)/(1-beta[3])*abs(y[i-1]-mean(y))
      }else{
        VaR[i] <- beta[4]+beta[5]*VaR[i-1]+(1-beta[5])*sqrt(beta[6]^2+
                                                              (1-beta[6])^2)/(beta[6])*abs(y[i-1]-mean(y))
      }
    }
    return(VaR)
  }
}
res <- TIMP_CAViaR(y,rep(0.1,6),0.05,flag=1)
beta <- c( -0.10989159,0.99663965,0.64075318,0.07742251 ,1.00000000 , 0.99037898)
beta <- c(-0.07163803,0.99999999,0.06793040,0.03892835,0.99353823,0.91642113)
ts.plot(ts(y),ts(res))
ts.plot(res)


