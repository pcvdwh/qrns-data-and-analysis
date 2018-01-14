# INITIALIZE 

  ddir <- '/put/data/directory/here/'
  setwd(ddir)

  library("survival")
  library("GGally")
  library("rstan")
  library("shinystan")
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

# FUNCTIONS

   # poisson confidence limits
   do.poissonlimits <- function(rate, min.ctr=5, max.ctr=500) {
      lambda <- seq(min.ctr, max.ctr, 1) * rate
     ci.99u <- qpois(0.995, lambda) 
     ci.95u <- qpois(0.975, lambda) 
     ci.50  <- qpois(0.500, lambda)
     ci.95l <- qpois(0.025, lambda) 
     ci.99l <- qpois(0.005, lambda) 
     uu <- cbind(lambda,ci.99u, c(0,ci.99u[1:(length(ci.99u)-1)]))
     uu <- uu[uu[,2] != uu[,3],]
     u  <- cbind(lambda,ci.95u, c(0,ci.95u[1:(length(ci.95u)-1)]))
     u  <- u[u[,2] != u[,3],]
     m  <- cbind(lambda,ci.50, c(0,ci.50[1:(length(ci.50)-1)]))
     m  <- m[m[,2] != m[,3],]
     l  <- cbind(lambda,ci.95l, c(0,ci.95l[1:(length(ci.95l)-1)]))
     l  <- l[l[,2] != l[,3],]
     ll <- cbind(lambda,ci.99l, c(0,ci.99l[1:(length(ci.99l)-1)]))
     ll <- ll[ll[,2] != ll[,3],]

     return(list(
      ll = cbind( E=ll[,1], SR=ll[,2]/ll[,1] ),
      l  = cbind( E= l[,1], SR= l[,2]/ l[,1] ),
      m  = cbind( E= m[,1], SR= m[,2]/ m[,1] ),
      u  = cbind( E= u[,1], SR= u[,2]/ u[,1] ),
      uu = cbind( E=uu[,1], SR=uu[,2]/uu[,1] )
      ))
  }

  # distance density equal
  dde.p <- function(r,p) {
    o <- order(r)
    r.o <- r[o]
    res <- unique(r.o)
    l <- length(r)
    ind <- c(1,p*(1:floor(l/p))) # index pts
    c <- 1
    for (i in 2:length(ind)) {
      res[res %in% r.o[ind[i-1]]:r.o[ind[i]]] <- c
      c=c+1
    }  
    res[res > r.o[ind[i]]] <- c    
    return(res)
  }


# READ DATA

  # VARIABLES
  # id       : hospital
  # yot      : year of treatment
  # sex      : gender (Vrouw=female, Man=male)
  # age      : age in years at surgery
  # age.n    : normalized age
  # kps      : Karnosfky perfroamce scale before surgery
  # kps.n    : normalized Karnofsky performance
  # tos      : type of surgery 
  # surv     : survival time in days
  # censor   : censor status (0=dead, 1=alive)
  # cause    : cause of death
  # causegrp : grouped death cause
  # deathloc : location of death
  # vol      : hospital volume in nrs of patients treated at a hospital in four years 
  # vol.n    : normalized volume
  # aca      : hospital academic setting
  # bxr      : hospital biopsy rate

  data <- read.csv("data.csv")
  data$id <- as.factor(data$id)
  data$yot <- as.factor(data$yot)
  str(data)

  # prepare complete cases dataset
  these.nona <- !(is.na(data$surv) | is.na(data$age.n) | is.na(data$kps.n))
  nrow(data) - sum(these.nona) # 101 excluded due to missingness
  censor <- !data$censor[these.nona] # 0:alive 1:dead 
  ctr   <- as.numeric(data$id[these.nona])
  age.n <- data$age.n[these.nona]
  kps.n <- data$kps.n[these.nona]
  yot   <- as.numeric(data$yot[these.nona])
  obs.t <- data$surv[these.nona]  # observation times with missings
  uq.t  <- unique(obs.t)          # unordered unique observation times
  t <- uq.t[order(uq.t)]          # ordered unique observation times
  ind.t <- match(obs.t,t)         # from observations to unique times
  t[length(t)] <- max(t) + 1      # add one to last time to distinguish dead from censored at last time point 
  period = dde.p(obs.t,112)       # equinr of pts per timeinterval, >100 pts
  nperiods = max(period)
  table(period)                   # nr of time intervals per period
  table(period[ind.t])            # nr of patients per period
  dt <- rep(NA,length(t))
  dt[1] = t[1];
  for (j in 2:length(t)) {
   dt[j] = t[j] - t[j - 1];  
  }
  t1m <- length(t)
  t2y <- length(t)
  for (j in 1:length(t)) {
    if (t[j] >= 30)  t1m = min(t1m, j); #  30d = 1m
    if (t[j] >= 730) t2y = min(t2y, j); # 730d = 2y
  }

  # list of data for input in stan model
  standata <- list(
    N  = length(obs.t),           # nr of pts
    M  = length(unique(data$id)), # nr of hospitals
    NT = length(t)-1,             # nr of unique failure times
    obs_t = obs.t,                # failure or censor time per patient
    ind = censor,                 # indicator: alive(censored) = 0, dead = 1
    ctr = ctr,                    # center per patient
    Zage  = age.n,                # normalized patient age
    Zkps  = kps.n,                # normalized patient kps
    yot   = yot,                  # year of treatment               
    ind_t = ind.t,                # corresponding unique time for observation time per patient
    t = t,                        # unique failure times + max followup time
    dt = dt,                      # length of interval before time t
    j1m = t1m,                    # time index corresponding to 1m
    j2y = t2y,                    # time index corresponding to 2y
    nperiods = nperiods,          # nr of periods
    period = period               # corresponding period to t
  )

  # inits for input in stan model
  nchains = 4
  initf <- function(chain_id = 1) {
    list( beta = rep(0.0, 2),
          eta  = rep(0.0, 4),
          zeta = rep(0.0, standata$M), 
          dL0  = rep(1.0, standata$nperiods), 
          tau  = 1 )
  }
  inits <- lapply(1:nchains, function(id) initf(chain_id=id))

  # run stan analysis, beware: takes quite some time
  stanresult <- stan(file = 'model.stan', 
              data = standata, 
              init = inits,
              pars = c(
                     "O_mg", "E_mg", "S_mg",
                     "O_1m", "E_1m", "S_1m",
                     "O_2y", "E_2y", "S_2y",
                     "O_1m_all", "E_1m_all", "S_1m_all", 
                     "O_2y_all", "E_2y_all", "S_2y_all", 
                     "tau_e","sigma_e","tau_z","sigma_z",
                     "dL0","beta","eta","zeta",
                     "V","V_1m","V_2y","R_1m","R_2y",
                     "Surv_all", "SurvM", "SurvMG", "Nrisk",
                     "Rhat_1m", "Rhat_2y"
                     ),
              iter = 5000, # iter = warmup + sample
              warmup = 1000, 
              chains = nchains
              )

  sum_stanresult <- summary(stanresult, pars=c(
                     "O_mg", "E_mg", "S_mg",
                     "O_1m", "E_1m", "S_1m",
                     "O_2y", "E_2y", "S_2y",
                     "O_1m_all", "E_1m_all", "S_1m_all", 
                     "O_2y_all", "E_2y_all", "S_2y_all", 
                     "tau_e","sigma_e","tau_z","sigma_z",
                     "dL0","beta","eta","zeta",
                     "V","V_1m","V_2y","R_1m","R_2y",
                     "Surv_all", "SurvM", "SurvMG", "Nrisk", 
                     "Rhat_1m","Rhat_2y"
                     ))

  stanresult_sso <- launch_shinystan(stanresult)

  # use distinctive colour scheme https://personal.sron.nl/~pault/
  mycols <- c( 
    rgb(136,46,114,maxColorValue=255),
    rgb(177,120,166,maxColorValue=255),
    rgb(214,193,222,maxColorValue=255),
    rgb(25,101,176,maxColorValue=255),
    rgb(82,137,199,maxColorValue=255),
    rgb(123,175,222,maxColorValue=255),
    rgb(78,178,101,maxColorValue=255),
    rgb(144,201,135,maxColorValue=255),
    rgb(202,224,171,maxColorValue=255),
    rgb(247,238,85,maxColorValue=255),
    rgb(246,193,65,maxColorValue=255),
    rgb(241,147,45,maxColorValue=255),
    rgb(232,96,28,maxColorValue=255),
    rgb(220,5,12,maxColorValue=255))

# EXTRACT RESULTS

  myt      <- standata$t[1:standata$NT]
  mySurv_all <- matrix(sum_stanresult$summary[grepl("\\bSurv_all\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE)
  mySurv_all_u <- matrix(sum_stanresult$summary[grepl("\\bSurv_all\\b",rownames(sum_stanresult$summary)),"97.5%"],ncol=(standata$NT),byrow=TRUE)
  mySurv_all_l <- matrix(sum_stanresult$summary[grepl("\\bSurv_all\\b",rownames(sum_stanresult$summary)),"2.5%"],ncol=(standata$NT),byrow=TRUE)
  mySurvM  <- matrix(sum_stanresult$summary[grepl("\\bSurvM\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE)
  mySurvMG <- matrix(sum_stanresult$summary[grepl("\\bSurvMG\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE)
  (V   <- sum_stanresult$summary[grepl("\\bV\\b",rownames(sum_stanresult$summary)),"mean"])
  myNr  <- matrix(sum_stanresult$summary[grepl("\\bNrisk\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE)
  dim(myNr)
  # nr at risk per 6m
  # corresponding observed days:  0 183,366,549,732,914,1095,1279,1464,1653,1827
  NRisk <- cbind(V, myNr[,which(myt %in% c(183,366,549,732,914,1095,1279,1464,1653,1827))])
  colnames(NRisk) <- seq(0,60,6)
  rownames(NRisk) <- 1:14

  (V1m <- sum_stanresult$summary[grepl("\\bV_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (V2y <- sum_stanresult$summary[grepl("\\bV_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (R1m <- sum_stanresult$summary[grepl("\\bR_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (R2y <- sum_stanresult$summary[grepl("\\bR_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (Rhat_1m <- sum_stanresult$summary[grepl("\\bRhat_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (O_1m_all <- sum_stanresult$summary[grepl("\\bO_1m_all\\b",rownames(sum_stanresult$summary)),"mean"])
  (Rhat_2y <- sum_stanresult$summary[grepl("\\bRhat_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (O_2y_all <- sum_stanresult$summary[grepl("\\bO_2y_all\\b",rownames(sum_stanresult$summary)),"mean"])

  (myMG_O <- matrix(sum_stanresult$summary[grepl("\\bO_mg\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE))
  (myMG_E <- matrix(sum_stanresult$summary[grepl("\\bE_mg\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE))
  (myMG_S <- matrix(sum_stanresult$summary[grepl("\\bS_mg\\b",rownames(sum_stanresult$summary)),"mean"],ncol=(standata$NT),byrow=TRUE))
  (myMG_OE <- myMG_O - myMG_E)

  (R_1m <- sum_stanresult$summary[grepl("\\bR_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (V_1m <- sum_stanresult$summary[grepl("\\bV_1m\\b",rownames(sum_stanresult$summary)),"mean"]) 
  (O_1m <- sum_stanresult$summary[grepl("\\bO_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (E_1m <- sum_stanresult$summary[grepl("\\bE_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (E_1m.lo <- sum_stanresult$summary[grepl("\\bE_1m\\b",rownames(sum_stanresult$summary)),"2.5%"])
  (E_1m.hi <- sum_stanresult$summary[grepl("\\bE_1m\\b",rownames(sum_stanresult$summary)),"97.5%"])
  (S_1m <- sum_stanresult$summary[grepl("\\bS_1m\\b",rownames(sum_stanresult$summary)),"mean"])
  (S_1m.lo <- sum_stanresult$summary[grepl("\\bS_1m\\b",rownames(sum_stanresult$summary)),"2.5%"])
  (S_1m.hi <- sum_stanresult$summary[grepl("\\bS_1m\\b",rownames(sum_stanresult$summary)),"97.5%"])

  (R_2y <- sum_stanresult$summary[grepl("\\bR_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (V_2y <- sum_stanresult$summary[grepl("\\bV_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (O_2y <- sum_stanresult$summary[grepl("\\bO_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (E_2y <- sum_stanresult$summary[grepl("\\bE_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (S_2y <- sum_stanresult$summary[grepl("\\bS_2y\\b",rownames(sum_stanresult$summary)),"mean"])
  (S_2y.lo <- sum_stanresult$summary[grepl("\\bS_2y\\b",rownames(sum_stanresult$summary)),"2.5%"])
  (S_2y.hi <- sum_stanresult$summary[grepl("\\bS_2y\\b",rownames(sum_stanresult$summary)),"97.5%"])

# TABLE 1

   library(survival)

   # calculate median OS per center
   mysurv <- survfit(Surv(d05$obs_t,d05$ind) ~ d05$ctr)
   printsurv <- read.table(textConnection(capture.output(mysurv)),skip=2,header=TRUE)
   mOS <- printsurv[,3]/30

   mymort_1m <- vector(length=M)
   mysurv_2y <- vector(length=M)
   for (k in 1:M) {
     obs_t_tmp <- d05$obs_t[d05$ctr==k]
     ind_tmp   <- d05$ind[d05$ctr==k]
     mysurv_tmp <- survfit(Surv(obs_t_tmp,ind_tmp) ~ 1)
     mymort_1m[k] <- 1 - mysurv_tmp$surv[min(which(mysurv_tmp$time>=30))]
     mysurv_2y[k] <- mysurv_tmp$surv[min(which(mysurv_tmp$time>=730))]
   }
  
  (tab1 <- cbind(
  rbind(
    table(data$id),
    table(data$id[these.nona]),
    table(data$sex,data$id),
    aggregate(data$age,list(ctr=data$id),function(x)mean(x,na.rm=TRUE))[,2],
    aggregate(data$age,list(ctr=data$id),function(x)sd(x,na.rm=TRUE))[,2],
    aggregate(data$age,list(ctr=data$id),function(x)sum(is.na(x)))[,2],
    table(data$kps,data$id,useNA='ifany')[c(10:1,11),],
    table(data$yot,data$id),
    c(rep("yes",7),rep("no",7)), #table(data$aca,data$id),
    table(data$tos,data$id)[c(2,1),],
    table(data$id,data$tos)[,"biopt"] / rowSums(table(data$id,data$tos)),
    mOS, 
    O_1m, mymort_1m,
    O_2y, mysurv_2y
    ),
  c(nrow(data),
    nrow(data[these.nona,]),
    table(data$sex),
    mean(data$age, na.rm=TRUE),
    sd(data$age, na.rm=TRUE),
    sum(is.na(data$age)),
    table(data$kps,useNA='ifany')[c(10:1,11)],
    table(data$yot),
    7,
    table(data$tos,useNA='ifany')[c(2,1)],
    table(data$tos)["biopt"] / sum(table(data$tos)),
    mOS_all, 
    O_1m_all, 1 - mySurv_all[min(which(myt>=30))], 
    O_2y_all, mySurv_all[min(which(myt>=730))]    
    )
  ))
  colnames(tab1) <- c(letters[1:14], "Overall")
  rownames(tab1) <- c(
    "N",
    "complete cases",
    "male","female",
    "mean age","s.d.","missing age",
    "KPS 100", "KPS 90", "KPS 80", "KPS 70", "KPS 60", "KPS 50", "KPS 40", "KPS 30", "KPS 20", "KPS 10", "KPS missing",
    "2011","2012","2013","2014",
    "academic setting",
    "resection","biopsy",
    "biopsy fraction",
    "median overall survival in months",
    "observed number of deaths at one month", 
    "mortality at one month",
    "observed number of survivors at two years", 
    "survival at two years",
  )
  write.table(tab1,"table1.txt",sep="\t")

# SUPPLEMENTAL FIGURE 2A & 2B

  myhospdat <- data.frame(
    aca = tab1["academic setting",1:14],
    vol = as.numeric(tab1["N",1:14]),
    bxr = as.numeric(tab1["biopsy fraction",1:14]),
    mos = as.numeric(tab1["median overall survival in months",1:14]))

  p <- ggpairs(myhospdat, 
    mapping = aes(col = aca,
                  alpha = 0.7),
    diag = list(continuous = wrap("blankDiag")),
    columnLabels = c("academic","volume in 4y","biopsy rate","median overall survival (months)"))
  col1 <- rgb(241,163, 64, maxColorValue = 255)
  col2 <- rgb(103,169,207, maxColorValue = 255)
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol){
      p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c(col1, col2)) +
        scale_color_manual(values=c(col1, col2))
    }
  }

  mypatdat <- data.frame(
    yot = data$yot,
    age = data$age,
    kps = data$kps * 10,
    surv = data$surv)
  q <- ggpairs(mypatdat, 
    mapping = aes(col = yot,
                  alpha = 0.1),
    diag = list(continuous = wrap("barDiag")),
    columnLabels = c("year","age","performance status","survival (days)"))
  fourcols <- 
    c(rgb(166,206,227, maxColorValue = 255),
      rgb(31,120,180, maxColorValue = 255),
      rgb(178,223,138, maxColorValue = 255),
      rgb(51,160,44, maxColorValue = 255))
  for(i in 1:q$nrow) {
    for(j in 1:q$ncol){
      q[i,j] <- q[i,j] + 
        scale_fill_manual(values=fourcols) +
        scale_color_manual(values=fourcols)
    }
  }

  cairo_ps('supp-fig-2a.ps')
    p
  dev.off()
  cairo_ps('supp-fig-2b.ps')
    q
  dev.off()

# FIGURE 1

  cairo_ps('fig-1.ps')
    par(mfrow=c(1, 1), mar = c(14,4,1,1))
    plot(survfit(Surv(standata$obs_t,standata$ind) ~ standata$ctr), col = mycols,
      xlab="time in months",xaxt='n',
      ylab="survival", las=1)
    axis(1,at=seq(0,11*182,182),labels=seq(0,11*6,6)) # time in months
    # survival function for all patients
    lines(c(0,myt), c(1,mySurv_all), col='black', lwd=2)     
    here <- seq(4,13,length.out=14)
    mtext("hospital", side=1, line=3, at = -100, cex=0.7, adj=1)
    mtext("number of patients at risk", side=1, line=3, at = -50, cex=0.7, adj=0)
    for (k in 1:14) {
      mtext(letters[k], side=1, line=here[k], at = -100, cex=0.7, adj=1, col=mycols[k])
      mtext(NRisk[k,], side=1, line=here[k], at=seq(0,10*182,182)+15, cex=0.7, adj=1)
    }
  dev.off()

# SUPPLEMENTAL FIGURE 1 

  cairo_ps('supp-fig-1.ps', width=11, height=11)

    par(mfrow=c(4, 4), mar = c(4,5,0.1,0.1))

    fitdata <- data.frame(obs.t=standata$obs_t, ind=standata$ind, Z=standata$ctr)

    # plot 1-14. centers
    for (k in 1:14) {
      fit <- survfit(Surv(obs.t,ind) ~ Z, data = fitdata, subset = (Z==k))
      #   kaplan meiers
      plot(fit, mark.time=TRUE, conf.int=FALSE,
         xlab="", xlim=c(0,3.1*365), xaxt='n',
         ylab="", las=1, lwd=2)
      abline(v=c(30,730),lty=3,col='grey',lwd=1)
      axis(1,at=seq(0,5*365,365),labels=seq(0,60,12),cex=1)
      legend("topright",legend=letters[k],lty=1,lwd=3,col=mycols[k], box.lwd=0, bg="white", cex=1.5)
      box()
      mtext("probability of survival", side=2, line=3, cex=0.8)
      mtext(NRisk[k,c(1,3,5,7)], side=1, line=2, at=seq(0,3*365,365), cex=0.6)
      mtext("nr at risk", side=1, line=2, at=-270, cex=0.6, adj=0)
      mtext("months", side=1, line=1, at=-270, cex=0.6, adj=0)
      lines(c(0,myt), c(1,mySurv_all), col='grey', lty=3, lwd=1)   # 'reference fit': survival function for all patients
      lines(c(0,myt), c(1,mySurvM[k,]), col='black', lty=1, lwd=1) # 'best fit': survival function for pts in center with specific age, kps, yot, ctr
      lines(c(0,myt), c(1,mySurvMG[k,]), col='blue', lty=1, lwd=2) # 'expected fit': survival function for pts in center with specific age, kps and average yot, ctr
    }

    # plot 15. 'martingale': O minus E
    plot(0,0,  
       xlab="",xlim=c(0,3.1*365),xaxt='n', las=1,
       ylab="", ylim=c(-22,22), col='white')
    mtext("time in months", side=1, line=3, cex=0.8)
    mtext("nr of deaths (observed - expected)", side=2, line=3, cex=0.8)
    axis(1,at=seq(0,5*365,365),labels=seq(0,60,12), las=1)
    abline(h=0,lty=1,col='black',lwd=1)
    for (k in 1:14) {
      lines(myt, myMG_OE[k,], col=mycols[k])
    }
    abline(v=c(30,730),lty=3,col='grey',lwd=1)

    # plot 16. 'martingale': O/E
    plot(0,0,  
       xlab="",xlim=c(0,3.1*365),xaxt='n', las=1,
       ylab="", ylim=c(-2,2),yaxt="n", col='white')
    mtext("time in months", side=1, line=3, cex=0.8)
    mtext("nr of deaths (observed / expected)", side=2, line=3, cex=0.8)
    axis(1,at=seq(0,5*365,365),labels=seq(0,60,12))
    myy <- c(1/8,1/4,1/2,1,2,4,8)
    axis(2,at=log(myy),labels=c("1/8","1/4","1/2","1","2","4","8"),las=1)
    abline(h=0,lty=1,col='black',lwd=1)
    for (k in 1:14) {
      lines(myt, log(myMG_S[k,]), col=mycols[k])
    }
    abline(v=c(30,730),lty=3,col='grey',lwd=1)


  dev.off()

# FIGURE 3

  dsmr <- data.frame(
    E  = E_1m,
    SR = S_1m,
    N  = V)
  lsmr <- do.poissonlimits(Rhat_1m)
  min(dsmr$E); max(dsmr$E)
  min(log(dsmr$SR)); max(log(dsmr$SR))

  dssr <- data.frame(
    E  = E_2y,
    SR = S_2y,
    N  = V)
  lssr <- do.poissonlimits(Rhat_2y)
  min(dssr$E); max(dssr$E)
  min(log(dssr$SR)); max(log(dssr$SR))

  cairo_ps('fig-3.ps')

  # 4 panels

     par(mfrow=c(2, 2))
     myx <- c(1/4,1/3,2/5,1/2,2/3,1,1.5,2,2.5,3,4)
     myy <- c(1/4,1/3,2/5,1/2,2/3,1,1.5,2,2.5,3,4)

  # panel 1 standardized mortality ratio
     par(mar=c(0,4,1,0))
     plot(log(dsmr$SR),dsmr$E,
          xlab='standardized mortality ratio',xlim=c(-1,+1),xaxt='n',yaxt='n',
          ylab='expected nr of early deaths at one month',ylim=c(0,20),
          col="white")
     axis(2, cex.axis=0.8, las=1)
     lines(log(lsmr$ll[,"SR"]),lsmr$ll[,"E"], col="black", lty=2)
     lines(log(lsmr$l[,"SR"]),lsmr$l[,"E"], col="black", lty=1)
     lines(log(lsmr$u[,"SR"]),lsmr$u[,"E"], col="black", lty=1)
     lines(log(lsmr$uu[,"SR"]),lsmr$uu[,"E"], col="black", lty=2)
     abline(v=log(1),col="grey", las=1)
     for (i in 1:nrow(dsmr)) {
       lines(c(log(dsmr$SR[i]),log(dsmr$SR[i])),
           c(dsmr$E[i],-5),
           col='grey',lty=3)
     }
     dsmr.cols <- rep('grey',14)
     dsmr.cols[2] <- 'black'
     points(log(dsmr$SR),dsmr$E,col=dsmr.cols,pch=19)
     # center labels
     text(log(dsmr$SR) + 0.05, dsmr$E + 0.5,labels=letters[1:14],col='black',cex=0.8)
     # panel label
     mtext("A",side=3, line=-2, at=log(2.5), cex=1.5) 
 
  # panel 2
     par(mar=c(0,0,0,0))
     plot(c(-0.7,0.7),c(-0.7,0.7),col='white',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
   
  # panel 3 smr vs ssr
     par(mar=c(4,4,0,0))
     plot(log(dsmr$SR),log(dssr$SR),
          xlab='risk-standardized early mortality ratio',xlim=c(-1,+1),xaxt='n',
          ylab='riak-standardized late survival ratio',ylim=c(-1,+1),yaxt='n',
          col="grey",pch='.')
     transp <- 70
     polygon(c(0,+1.2,+1.2,0),c(0,0,-1.6,-1.6),col=rgb( 222,45,38,transp,maxColorValue=255),border=NA) # red
     polygon(c(-1.2,-1.2,0,0),c(0,+1.6,+1.6,0),col=rgb( 49,163, 84,transp,maxColorValue=255),border=NA) # green
     abline(h=log(1),col="grey")
     abline(v=log(1),col="grey")
     axis(1,at=log(myx),labels=c("1/4","1/3","2/5","1/2","2/3","1","1.5","2","2.5","3","4"), las=1, cex.axis=0.8)
     axis(2,at=log(myy),labels=c("1/4","1/3","2/5","1/2","2/3","1","1.5","2","2.5","3","4"), las=1, cex.axis=0.8)
     # vertical
     for (i in 1:nrow(dsmr)) {
       lines(c(log(dsmr$SR[i]),log(dsmr$SR[i])),
           c(log(dssr$SR[i]),3),
           col='grey',lty=3)
     }
     # horizontal
     for (i in 1:nrow(dssr)) {
       lines(
         c(log(dsmr$SR[i]),3),
         c(log(dssr$SR[i]),log(dssr$SR[i])),
           col='grey',lty=3)
     }
     symbols(log(dsmr$SR),log(dssr$SR),circles=sqrt(dsmr$N/4)/100,inches=FALSE,add=TRUE)
     text(   log(dsmr$SR),log(dssr$SR),labels=letters[1:14],cex=0.8)

     legend.loc    <- c(-0.80,-0.90)
     legend.breaks <- c(20,50,100)*4
     r <- sqrt(legend.breaks/4)/100
     brk <- legend.breaks
     for (i in 1:length(r)) {
       symbols(legend.loc[1],
               legend.loc[2]+1*r[i],
               circle=r[i],inches=F,add=T,
               fg="grey")
       lines(c(legend.loc[1],legend.loc[1]+1.4*max(r)),
             rep(legend.loc[2]+1.95*r[i],2),
             col="grey")
       text(legend.loc[1]+1.6*max(r),
            legend.loc[2]+1.95*r[i], 
            paste(brk[i]/4,"per year"),adj=c(0,.5),
            cex=0.4,col="black")
     }   
     # panel label
     mtext("C",side=3, line=-2, at=log(2.5), cex=1.5) 

  # panel 4 standardized survival ratio
     par(mar=c(4,0,0,1))
     plot(dssr$E,log(dssr$SR),
          xlab='expected nr of late survivors at two years',xlim=c(0,60),
          ylab="",ylim=c(-1,+1),xaxt='n', yaxt='n',
          col="white")
     axis(1, cex.axis=0.8, las=1)
     lines(lssr$ll[,"E"],log(lssr$ll[,"SR"]), col="black", lty=2)
     lines(lssr$l[,"E"],log(lssr$l[,"SR"]), col="black", lty=1)
     lines(lssr$u[,"E"],log(lssr$u[,"SR"]), col="black", lty=1)
     lines(lssr$uu[,"E"],log(lssr$uu[,"SR"]), col="black", lty=2)
     abline(h=log(1),col="grey")
     # horizontal
     for (i in 1:nrow(dssr)) {
       lines(
         c(-5, dssr$E[i]),
         c(log(dssr$SR[i]),log(dssr$SR[i])),
           col='grey',lty=3)
     }
     dssr.cols <- rep('grey',14)
     dssr.cols[c(1,9,11,14)] <- 'black'
     points(dssr$E,log(dssr$SR),col=dssr.cols,pch=19)
     # center labels
     text(dssr$E - 1.5, log(dssr$SR) + 0.05,labels=letters[1:14],col='black',cex=0.8)
     # panel label
     mtext("B",side=3, line=-2, at=55, cex=1.5) 

dev.off()

# HOSPITAL CHARACTERISTICS RESULTS

# logistic regression with volume as in Spiegelhalter pmid 15568194 p1192

  # glm 
  # as example: http://ww2.coastal.edu/kingw/statistics/R-tutorials/logistic.html
  # Logistic Regression: One Numeric Predictor

  d1m.glm <- matrix( c(O1m,V1m-O1m), nc=2)
  fit.1m <- glm(d1m.glm ~ log(V1m), family=binomial(logit))  
  summary(fit.1m)
  # odds ratios with 95%ci
  exp(cbind(coef(fit.1m), confint(fit.1m)))
  newV <- 75:350
  newdata <- data.frame(V1m=newV)

  # changepoint
  exp((-log((1 / Rhat_1m) - 1) - coef(fit.1m)[1])/coef(fit.1m)[2])

  d2y.glm <- matrix( c(O2y,V2y-O2y), nc=2)
  fit.2y <- glm(d2y.glm ~ log(V2y), family=binomial(logit))  
  summary(fit.2y)
  # odds ratios with 95%ci
  exp(cbind(coef(fit.2y), confint(fit.2y)))
  newV <- 50:350
  newdata <- data.frame(V2y=newV)

  # changepoint
  exp((-log((1 / Rhat_2y) - 1) - coef(fit.2y)[1])/coef(fit.2y)[2])

# cox regression with volume

   library(survival)

   str(data)
   data$vol.b <- data$vol < 160

   fit <- survfit(Surv(surv,!censor) ~ age.n + kps.n + yot + id, data = data)

   summary(s00 <- coxph(Surv(surv,!censor) ~ frailty(id), data))

   summary(s01 <- coxph(Surv(surv,!censor) ~ age.n + kps.n + yot + frailty(id), data))
   pred.s01 <- predict(s01, newdata=data, type="expected")

   # regression coefficients from coxph very similar to Bayesian model
   sum_stanresult$summary[grepl("\\bbeta\\b",rownames(sum_stanresult$summary)),c("mean","2.5%","97.5%")]
   sum_stanresult$summary[grepl("\\beta\\b",rownames(sum_stanresult$summary)),c("mean","2.5%","97.5%")]
   sum_stanresult$summary[grepl("\\bzeta\\b",rownames(sum_stanresult$summary)),c("mean","2.5%","97.5%")]

   summary(s02 <- coxph(Surv(surv,!censor) ~ age.n + kps.n + yot, data))
   summary(s02 <- coxph(Surv(surv,!censor) ~ log(vol), data))
   summary(s02 <- coxph(Surv(surv,!censor) ~ log(vol) + age.n + kps.n + yot, data))

# logistic regression with academic setting

  aca <- c(1,1,1,1,1,1,1,0,0,0,0,0,0,0)
  fit.1m.aca <- glm(d1m.glm ~ aca, family=binomial(logit))  
  summary(fit.1m.aca)
  # odds ratios with 95%ci
  exp(cbind(coef(fit.1m.aca), confint(fit.1m.aca)))

  fit.2y.aca <- glm(d2y.glm ~ aca, family=binomial(logit))  
  summary(fit.2y.aca)
  # odds ratios with 95%ci
  exp(cbind(coef(fit.2y.aca), confint(fit.2y.aca)))

# cox regression with academic setting

  summary(s04 <- coxph(Surv(surv,!censor) ~ aca, data))
  summary(s04 <- coxph(Surv(surv,!censor) ~ aca + age.n + kps.n + yot, data))

# logistic regression with biopsy rate

  table(data$id,data$tos,useNA = 'always')
  BxR <- table(data$id,data$tos)[,"biopt"] / rowSums(table(data$id,data$tos))

  fit.1m.bxr <- glm(d1m.glm ~ BxR, family=binomial(logit))  
  summary(fit.1m.bxr)
  # odds ratios with 95%ci
  exp(cbind(coef(fit.1m.bxr), confint(fit.1m.bxr)))

  fit.2y.bxr <- glm(d2y.glm ~ BxR, family=binomial(logit))  
  summary(fit.2y.bxr)
  # odds ratios with 95%ci
  exp(cbind(coef(fit.2y.bxr), confint(fit.2y.bxr)))

  summary(s05 <- coxph(Surv(surv,!censor) ~ bxr, data))
  summary(s05 <- coxph(Surv(surv,!censor) ~ bxr + age.n + kps.n + yot, data))


# FIGURE 2

  mycexsm  = 0.7
  mycex    = 1
  mycexbig = 2
  myline   = 3

  cairo_ps('fig-2.ps', width=9, height=9)

    par(mfrow=c(2, 2), mar = c(5,5,0.1,0.1))

    # plot 1. volume vs mortality
    plot(V1m, R1m, col="white",
       xlim = c(0,350), xlab = "hospital volume in four years",  
       ylim = c(0,0.12), ylab = "observed early mortality rate at one month",
       las = 1, cex=mycex
    )
    # center points & labels
    mtext("A", side=2, line=myline, cex=mycexbig, las=1, at=0.12)
    points(V1m, R1m, pch=19, col=c(rep("black",7),rep("grey",7)))
    text(V1m+10, R1m+0.0025, labels=letters[1:14],cex=mycexsm)
    abline(Rhat_1m,0,col="black",lty=3)
    # lr line
    newV <- 60:350
    newdata <- data.frame(V1m=newV)
    preds <- predict(fit.1m, newdata=newdata, type="response")
    lines(newV, preds, type="l", col="black")

    # plot2. volume vs survival
    plot(V2y, R2y, col="white",
       xlim = c(0,350), xlab = "hospital volume in four years",  
       ylim = c(0,0.25), ylab = "observed late survival rate at two years",
       las = 1, cex=mycex
    )
    # center labels
    mtext("B", side=2, line=myline, cex=mycexbig, las=1, at=0.25)
    points(V2y, R2y, pch=19, col=c(rep("black",7),rep("grey",7)))
    text(V2y+10, R2y+0.005, labels=letters[1:14],cex=mycexsm)
    abline(Rhat_2y,0,col="black",lty=3)
    # lr line
    newV <- 30:330
    newdata <- data.frame(V2y=newV)
    preds <- predict(fit.2y, newdata=newdata, type="response")
    lines(newV, preds, type="l", col="grey")

    # plot3. biopsy vs mortality
    plot(BxR, R1m, col="white",
       xlim = c(0,1), xlab = "biopsy rate",  
       ylim = c(0,0.12), ylab = "observed early mortality rate at one month",
       las = 1, cex=mycex
    )
    # center labels
    mtext("C", side=2, line=myline, cex=mycexbig, las=1, at=0.12)
    points(BxR, R1m, pch=19, col=c(rep("black",7),rep("grey",7)))
    text(BxR+0.025, R1m+0.0025, labels=letters[1:14],cex=mycexsm)
    abline(Rhat_1m,0,col="black",lty=3)
    # lr line
    newBxR <- seq(0.1,0.8,0.01)
    newdata <- data.frame(BxR=newBxR)
    preds <- predict(fit.1m.bxr, newdata=newdata, type="response")
    lines(newBxR, preds, type="l", col="grey")

    # plot 4. biopsy vs survival
    plot(BxR, R2y, col="white",
       xlim = c(0,1), xlab = "biopsy rate",  
       ylim = c(0,0.25), ylab = "observed late survival rate at two years",
       las = 1, cex=mycex
    )
    # center labels
    mtext("D", side=2, line=myline, cex=mycexbig, las=1, at=0.25)
    points(BxR, R2y, pch=19, col=c(rep("black",7),rep("grey",7)))
    text(BxR+0.025,R2y+0.005, labels=letters[1:14],cex=mycexsm)
    abline(Rhat_2y,0,col="black",lty=3)
    # lr line
    newBxR <- seq(0.1,0.8,0.01)
    newdata <- data.frame(BxR=newBxR)
    preds <- predict(fit.2y.bxr, newdata=newdata, type="response")
    lines(newBxR, preds, type="l", col="grey")

  dev.off()

# cause of death & location

  dead_1m <- (data$surv <= 30 & data$censor == 0 & these.nona == TRUE)
  sum(dead_1m, na.rm=TRUE) # 119

  table(data$cause,useNA='always')
  (tabcause <- table(data$cause[dead_1m]))
  sum(tabcause) # 119

  table(data$causegrp,useNA='always')
  (tabcausegrp <- table(data$causegrp[dead_1m]))
  sum(tabcausegrp) # 119
  table(data$cause[data$causegrp == 'death directly related to intervention'])
  table(data$cause[data$causegrp == 'death indirectly related to intervention'])
  table(data$cause[data$causegrp == 'death unrelated'])

  table(data$deathloc,useNA='always')
  (tabloc <- table(data$deathloc[dead_1m]))
  

