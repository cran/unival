unival<-function(y, FP, fg, PHI, FA_model = "Lineal", type, SEP, SEG, relip, relig, percent = 90, display = TRUE){

  ######################################################################
  #  y : Related external variable
  ######################################################################

  if (missing(y)){
    stop("The argument y is not optional, please provide a valid related external variable")
  }

  y=as.matrix(y)

  ######################################################################
  #  FP : Primary factor score estimates
  ######################################################################

  if (missing(FP)){
    stop("The argument FP is not optional, please provide a valid primary factor score estimates")
  }

  FP=as.matrix(FP)


  ######################################################################
  #  PHI : Inter-Factor correlation matrix
  ######################################################################

    #compute gamma using PHI (correlation between factors)
    if (missing(PHI)==FALSE){
      PHI <- as.matrix(PHI)
      #Compute them using PHI
      k<-size(PHI)[2]
      if (k>3){
        #4 or more, factorizing minres
        gam <- psych::fa(r = PHI, nfactors = 1, fm = "minres", min.err = 0.0001)$loadings
        gam <- t(as.numeric(gam))
        #check if there are more positive than negative values
        pos <- sum(gam>0)
        neg <- sum(gam<0)
        if (neg>pos){
          gam <- -gam
        }
      }
      if (k==3){
        gam <- matrix(0,1,k)
        #triadas
        gam[1] <- sqrt(PHI[1,2]*PHI[1,3]/PHI[2,3])
        gam[2] <- sqrt(PHI[1,2]*PHI[2,3]/PHI[1,3])
        gam[3] <- sqrt(PHI[2,3]*PHI[1,3]/PHI[1,2])
      }
      if (k==2||k==1){ #if the user provides PHI like 2x2 matrix or as a single value
        # relip<-as.matrix(relip)
        # gam <- matrix(0,1,k)
        # tmp1 <- cor(FP[,1],fg)
        # tmp1 <- tmp1 / sqrt(relip[1]*relig)
        # tmp2 <- cor(FP[,2],fg)
        # tmp2 <- tmp2 / sqrt(relip[2]*relig)
        # if (k==2){
        #   tmp3 <- PHI[1,2]
        # }
        # else {
        #   tmp3 <- PHI
        # }
        #
        # beta1 <- (tmp1-(tmp2*tmp3)) / (1-(tmp3*tmp3))
        # beta2 <- (tmp2-(tmp1*tmp3)) / (1-(tmp3*tmp3))
        # gam[1] <- beta1 / (sqrt(1+(beta1*beta1)))
        # gam[2] <- beta2 / (sqrt(1+beta2*beta2))
      }
    }
  else {
    stop("The argument PHI is not optional, please provide a valid inter-factor correlation matrix")
  }

  ######################################################################
  #  FA_model : Which FA model was used for calibration and scoring. Available options are:
  #             - "Lineal" (by default)
  #             - "Graded"
  ######################################################################

  if (FA_model=="Lineal" || FA_model=="Graded") {
    if (missing(type)){
      warning("The type of factor scores was not provided, ML scores were assumed.")
      type="ML"
    }

    if (type=="EAP" || type=="ML"){}
    else {
      stop("The argument type has to be ML or EAP.")
    }
  }
  else{
      stop("The argument FA_model has to be Lineal of Graded.")
  }


  ######################################################################
  #  type : Which type of factor score estimates were used in FP and fg. If not specified, ML will be assumed
  #             - "ML"
  #             - "EAP"
  ######################################################################




  ######################################################################
  #  percent : Width of the confidence interval (by default 90 for 90% confidence interval)
  ######################################################################

  if (percent>99.99 || percent<=0){
    stop("percent argument has to be between 1 and 99")
  }
  else {
    tmp=(100-percent)/2
    cent1=tmp/100
    cent2=(100-tmp)/100
  }

  ######################################################################
  #  display : Determines if the output will be displayed in the console (TRUE by default)
  ######################################################################

  if (display!=F && display!=T){
    stop("display argument has to be logical (TRUE or FALSE, 0 or 1)")
  }


  ######################################################################
  #  relip : A vector containing the marginal reliabilities of the primary factor scores estimates. It is optional except when the number of factors is 2
  ######################################################################

  if (missing(relip)){

    if (k==2){
      stop("The argument relip is not optional when 2 primary factors are provided, please provide a vector containing the reliability estimates of the primary factor scores")
    }

    if (type =='EAP'){
      relip = var(FP) # EAP scores, the reliability estimation is the variance of the factor scores
    }
    else {
      relip = 1 / (var(FP)) # Bartlett scores
    }


    # if (max(var(FP)) < 1){
    #   relip = var(FP) # EAP scores, the reliability estimation is the variance of the factor scores
    # }
    # else {
    #   relip = 1 / (var(FP)) # Bartlett scores
    # }

    relip=diag(relip)
  }

  ######################################################################
  #  relig : The marginal reliability of the general factor (optional)
  ######################################################################

  if (missing(relig)){
    if (missing(fg)==FALSE){
      if (type == 'EAP'){
        relig = var(fg) # EAP scores, the reliability estimation is the variance of the factor scores
      }
      else {
        relig = 1 / (var(fg)) # Bartlett scores
      }
      # if (var(fg) < 1){
      #   relig = var(fg) # EAP scores, the reliability estimation is the variance of the factor scores
      # }
      # else {
      #   relig = 1 / (var(fg)) # Bartlett scores
      # }
    }
  }

  ######################################################################
  #  fg : General or second-order factor score estimates
  ######################################################################

  #Compute general factor scores
  if (missing(fg)){

    #Check if there are at least 3 primary factors
    k<-dim(FP)[2]
    if (k<3){
      stop("The argument fg is not optional when two factors were retained, please provide a valid first-order factor score estimates")
    }
    else {
      #compute fg & relig
      f1<-dim(gam)[1]
      f2<-dim(gam)[2]

      PSI <- matrix(1,k,1) - transpose(gam^2)
      PSI <- diag(c(PSI))

      if (type=='ML'){
        N<-dim(FP)[1]
        XR <- FP - matrix(1,N,1) %*% colMeans(FP)
        W <- solve(PSI)
        term <- solve(gam %*% W %*% transpose(gam))
        relig <- 1 / (1+term)
        fg <- XR %*% W %*% transpose(gam) %*% term

      }
      else {
        N <- dim(FP)[1]
        ML <- matrix(0,N,k)
        for (i in 1:k){
          ML[,i] <- FP[,i] / relip[i]
        }

        XR <- ML - matrix(1,N,1) %*% colMeans(ML)
        W <- solve(PSI)
        term <- solve((gam %*% W %*% transpose(gam)) +1)
        SEG <- sqrt(term) #psd
        relig <- 1 / (1+term)
        fg <- XR %*% W %*% transpose(gam) %*% term
      }

    }
  }

  fg=as.matrix(fg)

  ######################################################################
  #  SEP : Standard Errors (ML scores) or PSDs (EAP scores) for primary factor scores (only required when using graded model)
  ######################################################################

  ######################################################################
  #  SEG : Standard Errors (ML scores) or PSDs (EAP scores) for the general factor (only required when when using graded model)
  ######################################################################

  check_SE=F
  if (FA_model=="Graded"){
    if (missing(SEP)||(missing(SEG))){
      check_SE=F #if the scores are graded but the user did not provide errors
    }
    else {
      SEP <- as.matrix(SEP)
      SEG <- as.matrix(SEG)
      check_SE=T
    }
  }

  #################################    Checks     #################################

  #check if the relationship between the factors and the criterion is direct
  for (i in 1:k){
    if (cor(FP[,i],y)<0){
      #negative correlation
      FP[,i]=-FP[,i]
    }
  }
  if (cor(fg,y)<0){
    fg=-fg
  }

  #check the size of the inputs
  f1<-size(fg)[1]
  f2<-size(fg)[2]

  if (f2>f1){
    fg=t(fg)
  }
  if (f2>1){
    stop("The fg argument has to be a vector containing the factor scores of the general factor")
  }
  g1<-size(FP)[1]
  g2<-size(FP)[2]

  if (g2>g1){
    FP=t(FP)
  }

  if ((f1==g1)==FALSE){
    stop("The arguments fg and FP should have the same number of observations")
  }

  h1<-size(relip)[1]
  h2<-size(relip)[2]
  if (h1>h2){
    relip<-t(relip)
    h1<-size(relip)[1]
    h2<-size(relip)[2]
  }
  if ((h2==k)==FALSE){
    stop("The number of reliabilities provided in relip should be consistent with the number of retained factors")
  }
  j1<-size(relig)[1]
  j2<-size(relig)[2]
  if ((j1==1&&j2==1)==FALSE){
    stop("General reliability (relig) should be a single value")
  }

  if (check_SE==T){
    k1<-size(SEP)[1]
    k2<-size(SEP)[2]
    l1<-size(SEG)[1]
    l2<-size(SEG)[2]

    if (k2>k1){
      SEP<-t(SEP)
      k1<-size(SEP)[1]
      k2<-size(SEP)[2]
    }

    if ((k1==g1||l1==f1)==FALSE){
      stop('The arguments SEP, SEG, FP and fg should have the same number of observations')
    }

    if ((k2==k)==FALSE){
      stop("The number of columns of SEP should be consistent with the number of retained factors")
    }

    if (l2>1){
      stop("SEG should be a vector")
    }
    remove(k1,k2,l1,l2)
  }

  # Check if y is standarized
  tmp1<-round(mean(y), digits = 1)
  tmp2<-round(sd(y),digits = 1)
  if ((tmp1==0.0) & (tmp2==1.0)){
    #standarized, everything ok
  }
  else {
    y <- (y-mean(y))/sd(y)
  }

  remove(f1,f2,g1,g2,h1,h2,j1,j2,tmp1,tmp2)


  ################################# Everything  OK #################################
  ################################# Begin Analysis #################################


  N<-size(FP)[1]
  k<-size(FP)[2]

  if (k==2){

    dife<-matrix(0,1,k)
    difeg<-matrix(0,1,k)
    phi=PHI
    if (dim(PHI)[2]==2){
      phi=PHI[1,2]
    }

    #Disattenuated correlations: error corrections
    #in the primary factor score estimates

    for (i in 1:k){
      tmp1<-cor(FP[,i],y)
      tmp1<-tmp1/sqrt(relip[i])
      dife[i]<-tmp1
      remove(tmp1)
      tmp2<-cor(FP[,i],fg)
      tmp2<-tmp2/sqrt(relip[i])
      difeg[i]<-tmp2
      remove(tmp2)
    }

    #Error-in-variables multiple regressions

    betag1<-(difeg[1]-(difeg[2]*phi))/(1-(phi*phi))
    betag2<-(difeg[2]-(difeg[1]*phi))/(1-(phi*phi))
    betay1<-(dife[1]-(dife[2]*phi))/(1-(phi*phi))
    betay2<-(dife[2]-(dife[1]*phi))/(1-(phi*phi))

    # test for differential validity
    difev<-dife/difeg

    median_d<-abs(median(difev))

    extreme<-max(abs(difev))

    max_dif<-extreme-median_d

    #test for incremental validity: covariances
    R2g <- betag1 * dife[1] + betag2 * dife[2]
    R2y <- betay1 * dife[1] + betay2 * dife[2]

    betag <- c(betag1, betag2)
    betay <- c(betay1, betay2)

    #PHI like 2x2?
    denog <- sqrt(betag%*%PHI%*%transpose(betag))
    denoy <- sqrt(betay%*%PHI%*%transpose(betay))
    R2g <- R2g/denog
    R2y <- R2y/denoy
    contrast2 <- c(R2y,R2g)
    incre <- R2y-R2g

    #Dife Bootstrap

    difeboot <- matrix(0,500,k)
    difegboot <- matrix(0,500,k)
    for (i in 1:500){
      I <- round(runif(n = N,min = 1,max = N))
      outboot <- univalboot2(y[I],FP[I,],fg[I],relip)
      difeboot[i,]<-outboot$dife
      difegboot[i,]<-outboot$difeg
    }

    #difevboot<-sweep(difeboot,2,difegboot,FUN = '/',check.margin=FALSE)
    difevboot<-difeboot/difegboot

    CI<-matrix(0,2,k)
    for (i in 1:k){
      tmp <- quantile(difevboot[,i],c(cent1, cent2))
      CI[,i] <- t(tmp)
    }

    In<-which(colMeans(difevboot)==max(colMeans(difevboot)))
    tmp<-difevboot[,In] - median(colMeans(difevboot))
    CI_dif<-matrix(0,1,2)
    CI_dif[1]<-quantile(tmp,cent1)
    if (CI_dif[1]>max_dif){
      CI_dif[1]=max_dif
    }
    h1_dif=FALSE
    if (CI_dif[1]>0){
      h1_dif=TRUE
    }
    CI_dif[2]<-quantile(tmp,cent2)


    #Incremental bootstrap
    R2gboot<-matrix(0,500,1)
    R2yboot<-matrix(0,500,1)

    for (i in 1:500){

      I <- round(runif(n = N,min = 1,max = N))

      #test for incremental validity: covariances
      R2gboot[i] <- betag1 * difeboot[i,1] + betag2 * difeboot[i,2]
      R2yboot[i] <- betay1 * difeboot[i,1] + betay2 * difeboot[i,2]

      betag <- c(betag1, betag2)
      betay <- c(betay1, betay2)

      #PHI like 2x2?
      denog <- sqrt(betag%*%PHI%*%transpose(betag))
      denoy <- sqrt(betay%*%PHI%*%transpose(betay))
      R2gboot[i] <- R2gboot[i]/denog
      R2yboot[i] <- R2yboot[i]/denoy

    }

    contrast2boot <- cbind(R2yboot, R2gboot)
    increboot <- R2yboot - R2gboot

    # Confidence intervals
    CIc <- matrix(0,2,2)
    CIc[1,1] <- quantile(contrast2boot[,1],cent1)
    CIc[2,1] <- quantile(contrast2boot[,2],cent1)
    CIc[1,2] <- quantile(contrast2boot[,1],cent2)
    CIc[2,2] <- quantile(contrast2boot[,2],cent2)

    CI_incre <- matrix(0,2,1)
    CI_incre[1] <- quantile(increboot,cent1)
    if (CI_incre[1]>incre){
      CI_incre[1]=incre
    }
    h1_incre=FALSE
    if (CI_incre[1]>0){
      h1_incre=TRUE
    }
    CI_incre[2] <- quantile(increboot,cent2)


    dife=difev
    #output
    OUT<-list('differential_validity'=dife,'dife_CI'=CI,'max_dife'=max_dif,'max_dice_CI'=CI_dif,'contrast2'=contrast2,'contrast2_CI'=CIc,'incremental_validity'=incre,'incre_CI'=CI_incre)




  } #fi de if (k==2)
  else {

    dife<-matrix(0,1,k)

    if (check_SE==F){ #Lineal model, or graded when the user did not provide errors

      # Differential validity

      for (i in 1:k){
        tmp<-cor(FP[,i],y)
        tmp<-tmp/sqrt(relip[i])
        dife[i]<-tmp
        remove(tmp)
      }

      vrdis <- dife

      dife<-dife/gam

      median_d<-abs(median(dife))

      extreme<-max(abs(dife))

      max_dif<-extreme-median_d


      # Incremental validity
      gam=t(gam)

      tmpu <- matrix(1,k)
      uvec <- sqrt(tmpu - (gam*gam))

      CM <- gam%*%t(gam)
      CMT <- (CM - diag(diag(CM))) + diag(1,k)

      # Atenuation factor for the general-factor scores
      v <- gam / (uvec*uvec)
      den <- sum((gam*gam) / (uvec*uvec))
      v <- v/den
      nume <- t(v)%*%CM%*%v
      deno <- t(v)%*%CMT%*%v
      correcg <- sqrt(nume/deno)

      # Error-corrected correlation between the general factor and the criterion

      rdis <- cor(fg,y)
      rdis <- (rdis / sqrt(relig)) * correcg
      rdis <- abs(rdis)

      # Error-in-variables multiple correlation

      C <- cov(FP)
      m <- dim(C)[2]
      CT <- (C-diag(diag(C))) + diag(1,m)
      C2 <- cbind(y,FP)
      c2 <- cov(C2)
      c2 <- c2[,1]
      c2 <- c2[2:(k+1)]
      bt <- solve(CT)%*%c2
      rmult <- sqrt(vrdis%*%bt)

      # Test for incremental validity

      contrast2 <- c(rdis, rmult)
      incre <- rmult-rdis

      # Dife Bootstrap

      difeboot <- matrix(0,500,k)
      for (i in 1:500){
        I <- round(runif(n = N,min = 1,max = N))
        difeboot[i,] <- univalboot(y[I],FP[I,],relip)
      }

      vrdisboot<-difeboot

      difeboot<-sweep(difeboot,2,gam,FUN = '/')

      CI<-matrix(0,2,k)
      for (i in 1:k){
        tmp <- quantile(difeboot[,i],c(cent1, cent2))
        CI[,i] <- t(tmp)
      }

      In<-which(colMeans(difeboot)==max(colMeans(difeboot)))
      tmp<-difeboot[,In] - median(colMeans(difeboot))
      CI_dif<-matrix(0,1,2)
      CI_dif[1]<-quantile(tmp,cent1)
      if (CI_dif[1]>max_dif){
        CI_dif[1]=max_dif
      }
      h1_dif=FALSE
      if (CI_dif[1]>0){
        h1_dif=TRUE
      }
      CI_dif[2]<-quantile(tmp,cent2)




      #Incremental bootstrap
      rdisboot<-matrix(0,500,1)
      rmultboot<-matrix(0,500,1)

      for (i in 1:500){

        I <- round(runif(n = N,min = 1,max = N))

        #Error-corrected correlation between the general factor and the criterion

        tmp <- cor(fg[I],y[I])
        rdisboot[i] <- (tmp/sqrt(relig))*correcg
        rdisboot[i] <- abs(rdisboot[i])

        #Error-in-variables multiple correlation

        C <- cov(FP[I,])
        m <- dim(C)[2]
        CT <- (C-diag(diag(C))) + diag(1,m)
        C2 <- cbind(y[I],FP[I,])
        c2 <- cov(C2)
        c2 <- c2[,1]
        c2 <- c2[2:(k+1)]
        bt <- solve(CT)%*%c2
        if (vrdisboot[i,]%*%bt<0){ #Will produce a NaN when sqrt a negative value
          tmp <- 0
        }
        else {
          tmp <- sqrt(vrdisboot[i,]%*%bt)
        }


        rmultboot[i] <-tmp

      }

      contrast2boot <- cbind(rdisboot, rmultboot)
      increboot <- rmultboot - rdisboot

      # Confidence intervals
      CIc <- matrix(0,2,2)
      CIc[1,1] <- quantile(contrast2boot[,1],cent1)
      CIc[2,1] <- quantile(contrast2boot[,2],cent1)
      CIc[1,2] <- quantile(contrast2boot[,1],cent2)
      CIc[2,2] <- quantile(contrast2boot[,2],cent2)

      CI_incre <- matrix(0,2,1)
      CI_incre[1] <- quantile(increboot,cent1)
      if (CI_incre[1]>incre){
        CI_incre[1]=incre
      }
      h1_incre=FALSE
      if (CI_incre[1]>0){
        h1_incre=TRUE
      }
      CI_incre[2] <- quantile(increboot,cent2)



    }

    if (check_SE==T){ #graded model with standard errors or PSD

      # Differential validity

      PSDp <- SEP
      psdg <- SEG

      SPSDp <- PSDp*PSDp
      spsdg <- psdg*psdg
      tmpp <- mean(SPSDp)
      tmpg <- mean (spsdg)
      tmpp <- sqrt(matrix(1,k,1)+tmpp)
      tmpg <- sqrt(1+tmpg)

      for (i in 1:k){
        tmp <- cor(FP[,i],y)
        tmp <- tmp*tmpp[i]
        dife[i] <- tmp
        remove(tmp)
      }

      vrdis <- dife

      dife <- dife / gam

      median_d<-abs(median(dife))

      extreme<-max(abs(dife))

      max_dif<-extreme-median_d


      # Incremental validity
      gam=t(gam)

      tmpu <- matrix(1,k,1)
      uvec <- sqrt(tmpu-(gam*gam))

      CM <- gam%*%t(gam)
      CMT <- (CM-diag(diag(CM))) + diag(1,k)

      # Atenuation factor for the general-factor scores
      v <- gam / (uvec*uvec)
      den <- sum((gam*gam) / (uvec*uvec))
      v <- v/den
      nume <- t(v)%*%CM%*%v
      deno <- t(v)%*%CMT%*%v
      correcg <- sqrt(nume/deno)

      # Error-corrected correlation between the general factor and the criterion

      seap <- sd(fg) # apply(fg,2,sd)
      rdis <- cor(fg,y)
      rdis <- abs(rdis*tmpg*tmpg*seap*correcg)

      # Error-in-variables multiple correlation

      C2 <- cbind(y,FP)
      c2 <- cov(C2)
      c2 <- c2[,1]
      c2 <- c2[2:(k+1)]*tmpp*tmpp
      bt <- solve(CMT)%*%c2
      rmult <- sqrt(vrdis%*%bt)

      # Test for incremental validity

      contrast2 <- c(rdis, rmult)
      incre <- rmult-rdis


      # Dife Bootstrap

      difeboot <- matrix(0,500,k)
      for (i in 1:500){
        I <- round(runif(n = N,min = 1,max = N))
        difeboot[i,] <- univalboot(y[I],FP[I,],relip)
      }

      vrdisboot<-difeboot

      difeboot<-sweep(difeboot,2,gam,FUN = '/')

      CI<-matrix(0,2,k)
      for (i in 1:k){
        tmp <- quantile(difeboot[,i],c(cent1, cent2))
        CI[,i] <- t(tmp)
        if (CI[1,i]>dife[i]){
          CI[1,i]=dife[i]
        }
      }

      In<-which(colMeans(difeboot)==max(colMeans(difeboot)))
      tmp<-difeboot[,In] - median(colMeans(difeboot))
      CI_dif<-matrix(0,1,2)

      CI_dif[1]<-quantile(tmp,cent1)
      if (CI_dif[1]>max_dif){
        CI_dif[1]=max_dif
      }
      h1_dif=FALSE
      if (CI_dif[1]>0){
        h1_dif=TRUE
      }
      CI_dif[2]<-quantile(tmp,cent2)




      #Incremental bootstrap
      rdisboot<-matrix(0,500,1)
      rmultboot<-matrix(0,500,1)

      for (i in 1:500){


        I <- round(runif(n = N,min = 1,max = N))

        #Error-corrected correlation between the general factor and the criterion

        seap <- sd(fg[I]) # apply(fg,2,sd)
        rdisboot[i] <- cor(fg[I],y[I])
        rdisboot[i] <- abs(rdisboot[i]*tmpg*tmpg*seap*correcg)


        #Error-in-variables multiple correlation

        C2 <- cbind(y[I],FP[I,])
        c2 <- cov(C2)
        c2 <- c2[,1]
        c2 <- c2[2:(k+1)]*tmpp*tmpp
        bt <- solve(CMT)%*%c2
        if ((vrdisboot[i,]%*%bt)<0){ #Check
          rmultboot[i]=0
        }
        else {
          rmultboot[i] <- sqrt(vrdisboot[i,]%*%bt)
        }

        # C <- cov(FP[I,])
        # m <- dim(C)[2]
        # CT <- (C-diag(diag(C))) + diag(1,m)
        # C2 <- cbind(y[I],FP[I,])
        # c2 <- cov(C2)
        # c2 <- c2[,1]
        # c2 <- c2[2:(k+1)]
        # bt <- solve(CT)%*%c2
        # rmultboot[i] <- sqrt(vrdisboot[i,]%*%bt)

      }

      contrast2boot <- cbind(rdisboot, rmultboot)
      increboot <- rmultboot - rdisboot

      # Confidence intervals
      CIc <- matrix(0,2,2)
      CIc[1,1] <- quantile(contrast2boot[,1],cent1)
      if (CIc[1,1]>contrast2[1]){
        CIc[1,1]=contrast2[1]
      }
      CIc[2,1] <- quantile(contrast2boot[,2],cent1)
      if (CIc[2,1]>contrast2[2]){
        CIc[2,1]=contrast2[2]
      }
      CIc[1,2] <- quantile(contrast2boot[,1],cent2)
      CIc[2,2] <- quantile(contrast2boot[,2],cent2)

      CI_incre <- matrix(0,2,1)
      CI_incre[1] <- quantile(increboot,cent1)
      if (CI_incre[1]>incre){
        CI_incre[1]=incre
      }
      h1_incre=FALSE
      if (CI_incre[1]>0){
        h1_incre=TRUE
      }
      CI_incre[2] <- quantile(increboot,cent2)


    }

    #output
    OUT<-list('differential_validity'=dife,'differential_CI'=CI,'max_diffe'=max_dif,'max_diffe_CI'=CI_dif,'contrast2'=contrast2,'contrast2_CI'=CIc,'incremental_validity'=incre,'incremental_CI'=CI_incre)

  } # fi del else (k==2)

  if (display==TRUE){
    cat('\n')
    cat('Unival: Assessing essential unidimensionality using external validity information\n\n')

    cat('Differential validity assessment:\n\n')
    for (i in 1:k){
      cat(sprintf('%.4f (%.4f - %.4f) \n',dife[i],CI[1,i],CI[2,i]))
    }
    cat('\n')
    cat('Maximum difference\n\n')
    cat(sprintf('%.4f (%.4f - %.4f) ', max_dif,CI_dif[1],CI_dif[2]))
    if (h1_dif==T){
      cat('*')
    }
    cat('\n\n')

    cat('Incremental validity assessment:\n\n')
    cat(sprintf('%.4f (%.4f - %.4f) \n',contrast2[1],CIc[1,1],CIc[1,2]))
    cat(sprintf('%.4f (%.4f - %.4f)\n\n',contrast2[2],CIc[2,1],CIc[2,2]))
    cat('Incremental value estimate \n\n')
    cat(sprintf('%.4f (%.4f - %.4f) ', incre, CI_incre[1],CI_incre[2]))
    if (h1_incre==T){
      cat('**')
    }
    cat('\n\n')

    if (h1_dif==T){ #which.max(dife)
      cat('* Some factors are more strongly or weakly related to the criterion that can be predicted from their relations to the general factor\n')
    }
    if (h1_incre==T){
      cat('** There is a significant increase in accuracy between the prediction based on the primary factor score estimates and that based on the general factor score estimates.\n')
    }
    invisible(OUT)

  }
  else {
    return(OUT)
  }


}
