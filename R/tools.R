# tools.R

# The Kaplan-Meier distribution.
#
# Estimate various functions of the univariate distributions.
#
# ttilde The censored event time.
# deltatilde The event indicator.
# at The time point where the win-loss probabilities are evaluated.
unidistribution <- function(ttilde, deltatilde, at=NULL) {
  n <- length(ttilde)
  fit <- prodlim::prodlim(prodlim::Hist(ttilde,deltatilde)~1)
  time <- fit$time
  index <- time<=at
  time <- time[index]
  m <- length(time)
  S <- fit$surv[index]
  dN=fit$n.event[index]
  Y <- fit$n.risk[index]
  # Save only time, p, dH, H, dLambda on the interval [0,at].
  # Derived.
  p <- c(1-S[1],S[1:(m-1)]-S[2:m])    # Here we use that it is a distribution(?)
  H <- Y/n                            # Needed to compute dH. Note H[1]=1.
  dH <- c(H[2:m]-H[1:(m-1)], 0)       # Move the jump one to the left.

  dLambda <- (1/Y)*dN
  data.frame(time, p, H, dH, dLambda)
}

# Extending distribution.
#
# Extending distribution.
#
# @param dist The distribution.
# @param addtime The added times.
distextendtime <- function(dist, addtime) {
  time <- dist$time
  p <- dist$p
  H <- dist$H
  dH <- dist$dH
  dLambda <- dist$dLambda

  H_at0 <- H[1]

  # The addtime that is not in time.
  timedif <- setdiff(addtime,time)
  mc <- length(timedif)
  timec <- c(time,timedif)
  time_extend <- sort(timec)
  m <- length(time_extend)
  index <- match(time_extend, timec)
  p <- c(p, rep(0, mc))[index]
  dH <- c(dH, rep(0, mc))[index]
  # Moving it back.
  dH <- c(H_at0,dH[1:(m-1)])
  dLambda <- c(dLambda, rep(0, mc))[index]
  # addtimes <- c(rep(FALSE, length(time)), rep(TRUE, mc))[index]

  # Derived
  S <- 1-cumsum(p)
  Sm <- c(1, S[1:(m-1)])
  dS <- -p
  H <- cumsum(dH)
  invHLambda <- 1/(H*(1-dLambda))

  data.frame(time=time_extend, p, dH, dLambda, S, Sm, dS, H, invHLambda)
}


# The Kaplan-Meier distribution for the censoring distribution.
#
# Estimate various functions of the univariate distributions.
#
# @param ctilde The censored event time.
# @param cdeltatilde The event indicator.
# @param at The time point where the win-loss probabilites are evaluated.
unidistribution0 <- function(ctilde, cdeltatilde, at=NULL) {

  n <- length(ctilde)
  anycensoring <- (sum(cdeltatilde)>0)
  # Indicator for event.
  if(anycensoring) {
    fit <- prodlim::prodlim(prodlim::Hist(ctilde,cdeltatilde==0)~1, reverse=TRUE)
    G <- fit$surv  # reverse is only affecting the surv object.
  } else {
    fit <- prodlim::prodlim(prodlim::Hist(ctilde,cdeltatilde==0)~1)
    G <- rep(1, length(time))  # Defining the censoring survival function.
  }
  # Defining censoring specific.
  time <- fit$time
  index <- time<=at
  time <- time[index]
  m <- length(time)
  G <- G[index]
  dN <- fit$n.event[index]
  dN0 <- fit$n.lost[index]
  Y <- fit$n.risk[index]
  Y0 <- Y-dN

  # The idea is to save only time, p, dH, H, dLambda on the interval [0,at].
  # Derived.
  p <- c(1-G[1],G[1:(m-1)]-G[2:m])
  H <- Y0/n
  # Here we move the jump one to the left.
  dH <- c(H[2:m]-H[1:(m-1)],0)
  dLambda <- (1/Y0)*dN0
  data.frame(time, p, H, dH, dLambda)
}




# The recurrent event distribution.
#
# The recurrent event distribution.
#
# @param ttilde The selection from higher prioritized events.
# @param Ntilde2 The number of recurrent events.
# @param at The time point where the win-loss probabilities are evaluated.
# -> We carry on the Gt.
recurrent_distribution <- function(selection, Ntilde2, Gt, at=NULL) {

  # Recurrent distribution.
  Nk <- seq(0,max(Ntilde2[selection]))
  pN <- rep(NA, times=length(Nk))
  for(k in 1:length(Nk)) {
    pN[k] <- (1/Gt)*mean( (Ntilde2==Nk[k])*(selection) )
  }

  list(Nk=Nk, pN=pN, Gt=Gt)
}

recurrent_distribution_extend <- function(dist, addNk) {
  Nk <- dist$Nk
  pN <- dist$pN
  Gt <- dist$Gt
  Nkmax <- max(Nk)
  Nkmax_add <- max(addNk)
  if(Nkmax_add>Nkmax) {
    Nk <- c(Nk,seq(Nkmax+1,Nkmax_add))
    pN <- c(pN,rep(0,Nkmax_add-Nkmax))
  }
  list(Nk=Nk, pN=pN, Gt=Gt)
}


second_distribution <- function(selection, ttilde2, Gt, at=NULL) {

  n <- length(ttilde2)
  time <- sort(ttilde2[selection])
  pH_total <- sum(1*(selection)) / n
  pH <- tapply( time, time, length)
  pH <- unname(pH)
  pH <- pH / n
  index <- time<=at
  time <- time[index]
  pH <- pH[index]

  list(time=time, pH=pH, pH_total=pH_total, Gt=Gt)
}

second_distribution_extend <- function(dist, addtime) {

  time <- dist$time
  pH <- dist$pH
  pH_total <- dist$pH_total
  Gt <- dist$Gt

  # The addtime that is not in time.
  timedif <- setdiff(addtime,time)
  mc <- length(timedif)
  timec <- c(time,timedif)
  time_extend <- sort(timec)
  index <- match(time_extend, timec)
  pH <- c(pH, rep(0, mc))[index]
  m <- length(pH)
  # Derived
  H <- pH_total-cumsum(pH)
  Hm <- c(pH_total, H[1:(m-1)])
  p <- (1/Gt)*pH
  S <- (1/Gt)*H
  dS <- -(1/Gt)*pH

  list(time=time_extend, pH=pH, pH_total=pH_total, H=H, Hm=Hm, p=p, S=S, dS=dS, Gt=Gt)
}


#' Estimate win-loss parameters and their variance
#'
#' Estimate win-loss parameters and their variance.
#'
#' @param id The id variable.
#' @param time The time variable.
#' @param status The status variable coded: 0=censoring, 1=death, 2=secondary events, ...
#' @param group The group variable. Assuming groups are coded 0,1.
#' @param at The time point where the win-loss probabilities are evaluated.
#' @param type The type of analysis for each event type:
#'        1=comparing event times, 2=comparing number of recurrent events at time point at.
#' @param conf.level The level of the confidence interval.
#' @return A list containing: the win-loss parameter *wl* ordered as: first win, first loss,
#'         second win, second loss, ect. The asymptotic variance *sigma*.
#'         For all win-loss parameters, stub say, the fit contains the standard error,
#'         *se_stub*, lower and upper 95% confidence interval, *l_stub* and *u_stub*.
#'         The coverage of the confidence interval can be changed in the win-loss function.
#'         The parameters are:
#'         win-loss of each event type (*wl*);
#'         total win ratio and log win ratio (*wr* and *logwr*);
#'         total win difference (*wd*);
#'         event specific win ratio and log win ratio (*wrk* and *logwrk*);
#'         event specific win difference (*wrk* and *logwrk*);
#'         total win and loss (*w* and *l*);
#'         event specific ranking probability (*rankedk*);
#'         total ranking probability (*ranked*).
#' @export
#' @examples
#' winloss(id=hf_action$id,
#'         time=hf_action$time,
#'         status=hf_action$status,
#'         group=hf_action$group,
#'         type=c(1,2),
#'         at=47)
winloss <- function(id, time, status, group, at, type=NULL, conf.level = 0.95) {

  data <- data.frame(id=id, time=time, status=status, group=group)

  time_max <- max(time)

  # type=1,2. 1=event, 2=count.
  event_list <- sort(unique(subset(status, status!= 0 )))
  no_event_types <- length(event_list)
  no_parameters <- 2*no_event_types
  if(is.null(type)) {
    type <- rep(1, no_event_types)
  }

  # Check that all event except the last is single event.
  test_type <- length(unique(type[1:(no_event_types-1)]))==1
  if(!test_type) {
    stop("Only the last event is allowed to be recurrent.")
  }

  # survival event.
  data1 <- subset(data, status %in% c(0,1) )
  ttilde1 <- data1$time
  deltatilde1 <- data1$status
  # Group 1.
  ttilde11 <- ttilde1[data1$group==0]
  deltatilde11 <- deltatilde1[data1$group==0]
  n1 <- length(ttilde11)
  # Group 2.
  ttilde21 <- ttilde1[data1$group==1]
  deltatilde21 <- deltatilde1[data1$group==1]
  n2 <- length(ttilde21)

  # Defining phi and dphi.
  # Perhaps one should treat win-loss and dphi for each event type at a time.
  # The index is in the order i=individ, j=group, k=prioritized order.
  phi <- rep(NA, no_parameters)
  dphi1 <- matrix(NA, nrow=no_parameters, ncol=n1)
  dphi2 <- matrix(NA, nrow=no_parameters, ncol=n2)

  # Estimating distributions.
  # Survival distribution.
  dist1 <- unidistribution(ttilde11, deltatilde11, at=at)
  dist2 <- unidistribution(ttilde21, deltatilde21, at=at)
  dist1 <- distextendtime(dist1, dist2$time)
  dist2 <- distextendtime(dist2, dist1$time)
  # Group 1.
  time1 <- dist1$time
  m1 <- length(time1)
  p1 <- dist1$p
  S1 <- dist1$S
  Sm1 <- dist1$Sm
  dS1 <- dist1$dS
  dLambda1 <- dist1$dLambda
  H1 <- dist1$H
  invHLambda1 <- 1/(H1*(1-dLambda1))
  # Group 2.
  p2 <- dist2$p
  S2 <- dist2$S
  dS2 <- dist2$dS
  Sm2 <- dist2$Sm
  dLambda2 <- dist2$dLambda
  H2 <- dist2$H
  invHLambda2 <- 1/(H2*(1-dLambda2))

  # Censoring distribution.
  ctilde1 <- ttilde11
  cdeltatilde1 <- 1-deltatilde11
  ctilde2 <- ttilde21
  cdeltatilde2 <- 1-deltatilde21
  cens1 <- unidistribution0(ctilde1, cdeltatilde1, at=at)
  cens2 <- unidistribution0(ctilde2, cdeltatilde2, at=at)
  cens1 <- distextendtime(cens1, c(cens2$time, at) )
  cens2 <- distextendtime(cens2, c(cens1$time, at) )
  time0 <- cens1$time
  m0 <- length(time0)
  G1 <- cens1$S
  Gm1 <- cens1$Sm
  dLambda10 <- cens1$dLambda
  H10 <- cens1$H
  Gtm1 <- utils::tail(Gm1, 1)
  # Group 2.
  G2 <- cens2$S
  Gm2 <- cens2$Sm
  dLambda20 <- cens2$dLambda
  H20 <- cens2$H
  Gtm2 <- utils::tail(Gm2, 1)

  # Win and loss of the survival outcome.
  phi[1] <- sum(S2*p1) # win1
  phi[2] <- sum(S1*p2) # loss1

  # dphi(klj):
  # k=1,2 primary(survival)/secondary outcome.
  # l=1,2 win, loss.
  # j=1,2 group j component.

  # Win and loss of survival: group 1 component.
  # dphi111 and dphi121.
  dphi111 <- rep(NA,length=n1)
  dphi121 <- rep(NA,length=n1)
  for(i in 1:n1) {
    dNi11 <- 1*(ttilde11[i]==time1 & deltatilde11[i]==1)
    Yi11 <- 1*(time1<=ttilde11[i])
    Sprime1 <- -S1 *( cumsum(invHLambda1*dNi11) - cumsum(Yi11*invHLambda1*dLambda1) )
    dSprime1 <- c(Sprime1[1], Sprime1[2:m1]-Sprime1[1:(m1-1)])
    dphi111[i] <- -sum(S2*dSprime1)  # win1
    dphi121[i] <- -sum(Sprime1*dS2)  # win2
  }

  # Win and loss of survival: group 2 component.
  # dphi112 and dphi122.
  dphi112 <- rep(NA,length=n2)
  dphi122 <- rep(NA,length=n2)
  for(i in 1:n2) {
    dNi21 <- 1*(ttilde21[i]==time1 & deltatilde21[i]==1)
    Yi21 <- 1*(time1<=ttilde21[i])
    Sprime2 <- -S2 *( cumsum(invHLambda2*dNi21) - cumsum(Yi21*invHLambda2*dLambda2) )
    dSprime2 <- c(Sprime2[1], Sprime2[2:m1]-Sprime2[1:(m1-1)])
    dphi112[i] <- -sum(Sprime2*dS1)
    dphi122[i] <- -sum(S1*dSprime2)
  }

  dphi1[1,] <- dphi111
  dphi1[2,] <- dphi121
  dphi2[1,] <- dphi112
  dphi2[2,] <- dphi122

  # Define the selection after the survival event.
  selection1 <- ttilde11>at | (ttilde11==at & deltatilde11==0)
  selection2 <- ttilde21>at | (ttilde21==at & deltatilde21==0)

  for(k in 2:no_event_types) {
    # Single event.
    if(type[k]==1) {
      # Defining the time-to-event.
      data2_events <- subset(data, status %in% c(k) )
      data2 <- data.frame(id=data1$id, ttilde1, deltatilde1, group=data1$group)
      data2 <- merge(data2, data2_events, all=TRUE)
      data2$ttilde2 <- data2$time
      data2$deltatilde2 <- 1
      data2$deltatilde2[is.na(data2$time)] <- 0
      data2$ttilde2[is.na(data2$time) & data2$deltatilde1==1] <- time_max + 1
      data2$ttilde2[is.na(data2$time) & data2$deltatilde1==0] <- data2$ttilde1[is.na(data2$time) & data2$deltatilde1==0]
      ttilde2 <- data2$ttilde2
      deltatilde2 <- data2$deltatilde2
      ttilde12 <- ttilde2[data1$group==0]
      deltatilde12 <- deltatilde2[data1$group==0]
      ttilde22 <- ttilde2[data1$group==1]
      deltatilde22 <- deltatilde2[data1$group==1]

      # Event time distribution.
      second_dist1 <- second_distribution(selection1, ttilde12, deltatilde12, Gtm=Gtm1, at=at)
      second_dist2 <- second_distribution(selection2, ttilde22, deltatilde22, Gtm=Gtm2, at=at)
      second_dist1 <- second_distribution_extend(second_dist1, second_dist2$time)
      second_dist2 <- second_distribution_extend(second_dist2, second_dist1$time)
      # Group 1.
      time2 <- second_dist1$time
      m2 <- length(time2)
      H12 <- second_dist1$H
      H12m <- second_dist1$Hm
      pbiv1 <- second_dist1$p
      Sbiv1 <- second_dist1$S
      dSbiv1 <- second_dist1$dS
      # Group 2.
      H22 <- second_dist2$H
      H22m <- second_dist2$Hm
      pbiv2 <- second_dist2$p
      Sbiv2 <- second_dist2$S
      dSbiv2 <- second_dist2$dS

      phi[ 2*(k-1)+1 ] <- sum(Sbiv2*pbiv1)
      phi[ 2*(k-1)+2 ] <- sum(Sbiv1*pbiv2)

      # Win and loss of second outcome: group 1 component.
      # dphi211 and dphi121.
      dphi211 <- rep(NA,length=n1)
      dphi221 <- rep(NA,length=n1)

      for(i in 1:n1) {
        N12i <- 1*(selection1[i]==1 & ttilde12[i]>time2)
        N12im <- 1*(selection1[i]==1 & ttilde12[i]>=time2)
        N10i <- 1*(ctilde1[i]<=time0 & cdeltatilde1[i]==1)
        Y10i <- 1*(time0<ctilde1[i]) + 1*(ctilde1[i]==time0 & cdeltatilde1[i]==1)
        M10i <- N10i - cumsum(Y10i*dLambda10)
        dM10i <- c(M10i[1],M10i[2:m0]-M10i[1:(m0-1)])
        invLambdaHdM10i <- (1/((1-dLambda10)*H10))*dM10i
        invLambdaHdM10i[m0] <- 0
        intinvLambdaHdM10i <- sum( invLambdaHdM10i )
        Sbivprime1 <- (1/Gtm1)*(N12i-H12) + (1/Gtm1)*intinvLambdaHdM10i*H12
        Sbivprime1m <- (1/Gtm1)*(N12im-H12m) + (1/Gtm1)*intinvLambdaHdM10i*H12m
        dSbivprime1 <- Sbivprime1-Sbivprime1m
        dphi211[i] <- -sum(Sbiv2*dSbivprime1) # win2
        dphi221[i] <- -sum(Sbivprime1*dSbiv2) # loss2
      }

      # Win and loss of second outcome: group 2 component.
      # dphi211 and dphi121.
      dphi212 <- rep(NA,length=n2)
      dphi222 <- rep(NA,length=n2)
      for(i in 1:n2) {
        N22i <- 1*(selection2[i]==1 & ttilde22[i]>time2)
        N22im <- 1*(selection2[i]==1 & ttilde22[i]>=time2)
        N20i <- 1*(ctilde2[i]<=time0 & cdeltatilde2[i]==1)
        Y20i <- 1*(time0<ctilde2[i]) + 1*(ctilde2[i]==time0 & cdeltatilde2[i]==1)
        M20i <- N20i - cumsum(Y20i*dLambda20)
        dM20i <- c(M20i[1],M20i[2:m0]-M20i[1:(m0-1)])
        invLambdaHdM20i <- (1/((1-dLambda20)*H20))*dM20i
        invLambdaHdM20i[m0] <- 0
        intinvLambdaHdM20i <- sum( invLambdaHdM20i )
        Sbivprime2 <- (1/Gtm2)*(N22i-H22) +  (1/Gtm2)*intinvLambdaHdM20i*H22
        Sbivprime2m <- (1/Gtm2)*(N22im-H22m) +  (1/Gtm2)*intinvLambdaHdM20i*H22m
        dSbivprime2 <- Sbivprime2-Sbivprime2m
        dphi212[i] <- -sum(Sbivprime2*dSbiv1) # win2
        dphi222[i] <- -sum(Sbiv1*dSbivprime2) # loss2
      }
      dphi1[2*(k-1)+1,] <- dphi211
      dphi1[2*(k-1)+2,] <- dphi221
      dphi2[2*(k-1)+1,] <- dphi212
      dphi2[2*(k-1)+2,] <- dphi222

      # Update the selection
      selection1 <- selection1 & (ttilde12>at | (ttilde12==at & deltatilde12==0))
      selection2 <- selection2 & (ttilde22>at | (ttilde22==at & deltatilde22==0))
    }

    # Secondary count.
    if(type[k]==2) {
      Ntilde <- with(data, tapply( 1*(status==2 & time<=at), id, sum))
      Ntilde12 <- Ntilde[data1$group==0]
      Ntilde22 <- Ntilde[data1$group==1]

      # Recurrent event distribution.
      recurrent_dist1 <- recurrent_distribution(selection1, Ntilde12, Gtm=Gtm1, at=at)
      recurrent_dist2 <- recurrent_distribution(selection2, Ntilde22, Gtm=Gtm2, at=at)
      recurrent_dist1 <- recurrent_distribution_extend(recurrent_dist1, recurrent_dist2$Nk)
      recurrent_dist2 <- recurrent_distribution_extend(recurrent_dist2, recurrent_dist1$Nk)
      Nk <- recurrent_dist1$Nk
      mN <- length(Nk)
      # Group 1.
      pN1 <- recurrent_dist1$pN
      Gtm1 <- recurrent_dist1$Gtm
      FN1 <- cumsum(pN1)
      FmN1 <- c(0,FN1[1:(mN-1)])
      pH1 <- Gtm1*pN1
      FH1 <- cumsum(pH1)
      FmH1 <- c(0,FH1[1:(mN-1)])
      # Group 2.
      pN2 <- recurrent_dist2$pN
      Gtm2 <- recurrent_dist2$Gtm
      FN2 <- cumsum(pN2)
      FmN2 <- c(0,FN2[1:(mN-1)])
      pH2 <- Gtm2*pN2
      FH2 <- cumsum(pH2)
      FmH2 <- c(0,FH2[1:(mN-1)])

      phi[ 2*(k-1)+1 ] <- sum(FmN2*pN1)
      phi[ 2*(k-1)+2 ] <- sum(FmN1*pN2)

      # Win and loss of second recurrent outcome: group 1 component.
      # dphi211 and dphi121.
      dphi211 <- rep(NA,length=n1)
      dphi221 <- rep(NA,length=n1)
      for(i in 1:n1) {
        pN1i <- 1*(selection1[i]==1 & Ntilde12[i]==Nk)
        FN1i <- 1*(selection1[i]==1 & Ntilde12[i]<=Nk)
        FmN1i <- 1*(selection1[i]==1 & Ntilde12[i]<Nk)
        N10i <- 1*(ctilde1[i]<=time0 & cdeltatilde1[i]==1)
        Y10i <- 1*(time0<ctilde1[i]) + 1*(ctilde1[i]==time0 & cdeltatilde1[i]==1)
        M10i <- N10i - cumsum(Y10i*dLambda10)
        dM10i <- c(M10i[1],M10i[2:m0]-M10i[1:(m0-1)])
        invLambdaHdM10i <- (1/((1-dLambda10)*H10))*dM10i
        invLambdaHdM10i[m0] <- 0
        intinvLambdaHdM10i <- sum( invLambdaHdM10i )
        pN1prime <- (1/Gtm1)*(pN1i-pH1) +  (1/Gtm1)*intinvLambdaHdM10i*pH1
        FN1prime <- (1/Gtm1)*(FN1i-FH1) +  (1/Gtm1)*intinvLambdaHdM10i*FH1
        FmN1prime <- (1/Gtm1)*(FmN1i-FmH1) +  (1/Gtm1)*intinvLambdaHdM10i*FmH1
        dphi211[i] <- sum(FmN2*pN1prime)  # win2
        dphi221[i] <- sum(FmN1prime*pN2)  # loss2
      }

      # Win and loss of second recurrent outcome: group 2 component.
      # dphi211 and dphi121.
      dphi212 <- rep(NA,length=n2)
      dphi222 <- rep(NA,length=n2)
      i <- 1
      for(i in 1:n2) {
        pN2i <- 1*(selection2[i]==1 & Ntilde22[i]==Nk)
        FN2i <- 1*(selection2[i]==1 & Ntilde22[i]<=Nk)
        FmN2i <- 1*(selection2[i]==1 & Ntilde22[i]<Nk)
        N20i <- 1*(ctilde2[i]<=time0 & cdeltatilde2[i]==1)
        # xxx
        Y20i <- 1*(time0<ctilde2[i]) + 1*(ctilde2[i]==time0 & cdeltatilde2[i]==1)
        M20i <- N20i - cumsum(Y20i*dLambda20)
        dM20i <- c(M20i[1],M20i[2:m0]-M20i[1:(m0-1)])
        invLambdaHdM20i <-  (1/((1-dLambda20)*H20))*dM20i
        invLambdaHdM20i[m0] <- 0
        intinvLambdaHdM20i <- sum( invLambdaHdM20i )
        pN2prime <- (1/Gtm2)*(pN2i-pH2) +  (1/Gtm2)*intinvLambdaHdM20i*pH2
        FN2prime <- (1/Gtm2)*(FN2i-FH2) +  (1/Gtm2)*intinvLambdaHdM20i*FH2
        FmN2prime <- (1/Gtm2)*(FmN2i-FmH2) +  (1/Gtm2)*intinvLambdaHdM20i*FmH2
        dphi212[i] <- sum(FmN2prime*pN1)  # win2
        dphi222[i] <- sum(FmN1*pN2prime)  # loss2
      }

      dphi1[2*(k-1)+1,] <- dphi211
      dphi1[2*(k-1)+2,] <- dphi221
      dphi2[2*(k-1)+1,] <- dphi212
      dphi2[2*(k-1)+2,] <- dphi222
    }

  }

  # Sigma1.
  #           win1  loss1   win2   loss2
  #           n1*(2*number of prioritized events)
  sigma1 <- dphi1 %*% t(dphi1) / n1

  # Sigma2.
  #           win1  loss1   win2   loss2
  #           n2*(2*number of prioritized events)
  sigma2 <- dphi2 %*% t(dphi2) / n1

  # Sigma.
  sigma <- sigma1/n1 + sigma2/n2

  # Derived quantaties.
  factor <- stats::qnorm(1-(1-conf.level)/2 )
  no_parameters <- length(phi)
  no_event_types <- no_parameters/2
  index_win_list <- 1+(0:(no_event_types-1))*2
  index_loss_list <- 2+(0:(no_event_types-1))*2
  index_win <- rep(FALSE,no_parameters)
  index_win[index_win_list] <- TRUE
  index_loss <- rep(FALSE,no_parameters)
  index_loss[index_loss_list] <- TRUE

  # Win-loss.
  wl <- phi
  var_wl <- sigma
  se_wl <- sqrt(sigma[cbind(1:no_parameters,1:no_parameters)])
  l_wl <- wl - factor*se_wl
  u_wl <- wl + factor*se_wl

  # log(wr).
  wr <- sum(phi[index_win])/sum(phi[index_loss])
  logwr <-  log(wr)
  gprime <- rep(0, no_parameters)
  gprime[index_win] <- 1/sum(wl[index_win])
  gprime[index_loss] <- -1/sum(wl[index_loss])
  gprime <- matrix(gprime, ncol=1, nrow=no_parameters)
  var_logwr <- t(gprime) %*% sigma %*% gprime
  se_logwr <- c(sqrt(var_logwr))
  l_wr <- exp(logwr - factor*se_logwr)
  u_wr <- exp(logwr + factor*se_logwr)

  # wd.
  wd <- sum(phi[index_win]) - sum(phi[index_loss])
  gprime <- rep(0, no_parameters)
  gprime[index_win] <- 1
  gprime[index_loss] <- -1
  gprime <- matrix(gprime, ncol=1, nrow=no_parameters)
  var_wd <- t(gprime) %*% sigma %*% gprime
  se_wd <- sqrt(var_wd)
  l_wd <- wd - factor*se_wd
  u_wd <- wd + factor*se_wd

  # log(wr1), log(wr2) ect.
  wrk <- phi[index_win]/phi[index_loss]
  logwrk <- log(wrk)
  gprime <- matrix(0, ncol=no_event_types, nrow=no_parameters)
  for(i in 1:no_event_types) {
    gprime[index_win_list[i],i] <- 1/phi[index_win_list[i]]
    gprime[index_loss_list[i],i] <- -1/phi[index_loss_list[i]]
  }
  var_logwrk <- t(gprime) %*% sigma %*% gprime
  se_logwrk <- sqrt(diag(var_logwrk))
  l_wrk <- exp(logwrk - factor*se_logwrk)
  u_wrk <- exp(logwrk + factor*se_logwrk)

  # wd1, wd2, ect.
  wdk <- phi[index_win]-phi[index_loss]
  gprime <- matrix(0, ncol=no_event_types, nrow=no_parameters)
  for(i in 1:no_event_types) {
    gprime[index_win_list[i],i] <- 1
    gprime[index_loss_list[i],i] <- -1
  }
  var_wdk <- t(gprime) %*% sigma %*% gprime
  se_wdk <- sqrt(diag(var_wdk))
  l_wdk <- wdk - factor*se_wdk
  u_wdk <- wdk + factor*se_wdk

  # w and l.
  w <- sum(phi[index_win])
  gprime <- rep(0, no_parameters)
  gprime[index_win] <- 1
  gprime <- matrix(gprime, ncol=1, nrow=no_parameters)
  var_w <- t(gprime) %*% sigma %*% gprime
  se_w <- sqrt(var_w)
  l_w <- w - factor*se_w
  u_w <- w + factor*se_w
  l <- sum(phi[index_loss])
  gprime <- rep(0, no_parameters)
  gprime[index_loss] <- 1
  gprime <- matrix(gprime, ncol=1, nrow=no_parameters)
  var_l <- t(gprime) %*% sigma %*% gprime
  se_l <- sqrt(var_l)
  l_l <- l - factor*se_l
  u_l <- l + factor*se_l

  # ranked.
  rankedk <- rep(NA, no_event_types)
  se_rankedk <- rep(NA, no_event_types)
  l_rankedk <- rep(NA, no_event_types)
  u_rankedk <- rep(NA, no_event_types)
  for(i in 1:no_event_types) {
    gprime <- rep(0, no_parameters)
    index <- 1:2+2*(i-1)
    rankedk[i] <- sum(wl[index])
    gprime[index] <- 1
    gprime <- matrix(gprime, ncol=1, nrow=no_parameters)
    var <- t(gprime) %*% sigma %*% gprime
    se <- sqrt(var)
    se_rankedk[i] <- se
    l_rankedk[i] <- rankedk[i] - factor*se
    u_rankedk[i] <- rankedk[i] + factor*se
  }
  ranked <- sum(wl)
  gprime <- matrix(rep(1,no_parameters),
                   ncol=1, nrow=no_parameters)
  var_ranked <- c(t(gprime) %*% sigma %*% gprime)
  se_ranked=sqrt(var_ranked)
  l_ranked <- ranked - factor*se_ranked
  u_ranked <- ranked + factor*se_ranked

  result <- list(wl=wl,se_wl=se_wl,l_wl=l_wl,u_wl=u_wl, sigma=sigma,
                 wr=wr,logwr=logwr,se_logwr=se_logwr,l_wr=l_wr,u_wr=u_wr,
                 wd=wd,se_wd=se_wd,l_wd=l_wd,u_wd=u_wd,
                 wrk=wrk,logwrk=logwrk,se_logwrk=se_logwrk,l_wrk=l_wrk,u_wrk=u_wrk,
                 wdk=wdk,se_wdk=se_wdk,l_wdk=l_wdk,u_wdk1=u_wdk,
                 w=w,se_w=se_w,l_w=l_w,u_w=u_w,
                 l=l,se_l=se_l,l_l=l_l,u_l=u_l,
                 rankedk=rankedk,se_rankedk=se_rankedk,l_rankedk=l_rankedk,u_rankedk=u_rankedk,
                 ranked=ranked,se_ranked=se_ranked,l_ranked=l_ranked,u_ranked=u_ranked)
  class(result) <- "winloss"
  result
}

