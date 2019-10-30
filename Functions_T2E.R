##############################################################################################
# Title: Functions for composite time to event endpoints in randomized controlled trials 
# Authors: Moises Gomez Mateu, Jordi Cortes Martinez, Marta Bofill Roig, Guadalupe Gomez Melis
# Last updated: 2019-07-09
# 
# Software version: R version 3.6.1 (2019-10-28)
# Depends:
# reshape
# data.table 
# ggplot2
# copula
# numDeriv
# rootSolve
# DiceKriging
##############################################################################################

#######################################################################################
# Function: ARE 
#
#######################################################################################
# Description: It computes ARE (Assymptotic Relative Efficiency)
# Parameters:
# rho0	                Spearman's coefficient between T1 and T2 in control group
# rho1	                Spearman's coefficient between T1 and T2 in treatment group
# beta1	                Shape parameter for a Weibull law for the relevant event
# beta2                 Shape parameter for a Weibull law for the additional event 
# HR1                   Hazard Ratio for a Weibull law for the relevant event
# HR2                   Hazard Ratio for a Weibull law for the additional event
# p1                    Proportion of the relevant event expected in group zero
# p2                    Proportion of the additional event expected in group zero
# case                  Censoring case (1,2,3 or 4)
# copula                Copula used:
#                          Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#                          Elliptical: "Normal"
#                          Other: "FGM"
# rho_thype             Type of correlation (Spearman or Kendall)
#######################################################################################

ARE <- function(rho0,rho1=rho0, beta1, beta2, HR1, HR2, p1, p2, case = case, copula = copula, rhoType='Spearman'){ 
  
  ############################################################
  ##-- 1. SELECTION OF THE COPULA
  copula0 <- CopulaSelection(copula,rho0,rhoType)
  theta <- copula0[[2]]   
  which.copula0 <- copula0[[1]]
  which.copula1 <- CopulaSelection(copula,rho1,rhoType)[[1]]  
  ############################################################

  ############################################################
  ##-- 2. SELECTION OF THE MARGINAL DISTRIBUTIONS
  MarginSelec <- MarginalsSelection(beta1,beta2,HR1,HR2,p1,p2,case,theta,copula=copula)
  T1dist <- MarginSelec[[1]]
  T2dist <- MarginSelec[[2]]
  T1pdist <- MarginSelec[[3]]
  T2pdist <- MarginSelec[[4]]
  T10param <- MarginSelec[[5]]
  T20param <- MarginSelec[[6]]
  T11param <- MarginSelec[[7]]
  T21param <- MarginSelec[[8]]
  ############################################################
  
  ############################################################
  ###### 3. ARE EXPRESSION FOLLOWING: 
  # Gomez G, Lagakos SW. Statistical considerations when using a composite endpoint for comparing treatment groups. 
  # Stat Med. 2013; 32:719–738 (pages 2 and 3).
  
  # Bivariate distribution in control and treatment groups
  distribution0 <- mvdc(copula = which.copula0, margins = c(T1dist, T2dist), paramMargins = list(T10param, T20param))
  distribution1 <- mvdc(copula = which.copula1, margins = c(T1dist, T2dist), paramMargins = list(T11param, T21param))
  
  if(case==1|case==3) {
    
    inside_integral <- function(t){
      
      Sstar0 <- Sstar(x = t,dist1 = T1pdist,dist2 = T2pdist, param1 = T10param,param2 = T20param, dist_biv = distribution0)
      Sstar1 <- Sstar(x = t,dist1 = T1pdist,dist2 = T2pdist, param1 = T11param,param2 = T21param, dist_biv = distribution1)
      
      fstar0 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T10param, param2 = T20param,dist_biv = distribution0))
      fstar1 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T11param, param2 = T21param,dist_biv = distribution1))
      
      Lstar0 <- (fstar0/Sstar0)
      Lstar1 <- (fstar1/Sstar1)
      
      HRstar <- (Lstar1/Lstar0)
      
      logHRstar <- log(HRstar)
      
      ##-- Numerical issues
      fstar0[fstar0<0] <- 0
      logHRstar[is.na(logHRstar) | logHRstar==Inf | logHRstar== -Inf] <- 0
      
      return(logHRstar*fstar0)
    }
    
    # Numerator
    integral <- integrate(inside_integral,lower=0,upper=1,subdivisions=10000,stop.on.error = FALSE) 
    numerator <-(integral$value)^2
    
    # Denominator
    Sstar0_1 <- Sstar(x=1,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0) 
    ST10_1 <- 1-do.call(T1pdist,c(q=1,T10param)) 
    denominator <- ((log(HR1))^2)*(1-Sstar0_1)*(1-ST10_1)
    
    # ARE value
    AREstarT <- (numerator/denominator)
    
    # If the integral is not computed, we assign a missing value
    if(integral$message!="OK") {AREstarT <- NA}
    
  } else
    if(case==2|case==4) {
      
      # Computation of the scale parameter values: b10, b20 
      if(case==2) {
        
        # Compute b20
        b20 <- 1/(-log(1-p2))^(1/beta2)
        
        # Compute b10
        Fb10 <- function(b10,p1){
          integral<-integrate(function(u) {
            sapply(u, function(u) {
              integrate(function(v) (	(theta*(1-exp(-theta))*exp(-theta*(u+v)))/ (exp(-theta)+  exp(-theta*(u+v))  - exp(-theta*u)-exp(-theta*v))^2 )    , lower=0, upper=   exp((b10*(-log(u))^(1/beta1))^beta2*log(1-p2)) 	)$value
            })
          }, lower= exp(-1/b10^beta1), upper=1)$value
          return(integral-p1)
        }
        
        limits <- c(0.00001,10000)                        # The first and the last values must be in opposite signs for the function
        b10 <- uniroot(Fb10, interval=limits,p1=p1)$root  # Find the root (value which equals the function to zero)
        
      }
      
      if(case==4) {
        # We need to create x[1] and x[2] to run 'multiroot' function (library: rootSolve) (NA's initially assigned)
        x<-NA
        y<-NA
        x[1]<-x
        x[2]<-y
        
        # We need to change the name of variables as (b10=x[1],b20=[2]) to execute 'multiroot'
        # Compute b10
        Fb10 <- function(b10,b20,p1){
          b10 -> x[1]
          b20 -> x[2]
          integral<-integrate(function(u) {
            sapply(u, function(u) {
              integrate(function(v) (	(theta*(1-exp(-theta))*exp(-theta*(u+v)))
                                      /(exp(-theta)+  exp(-theta*(u+v))  - exp(-theta*u)-exp(-theta*v))^2 )    , lower=0,
                        upper=   exp((x[1]*(-log(u))^(1/beta1))^beta2*(  -1/(x[2]^beta2)       )) 	)$value
            })
          }, lower= exp(-1/x[1]^beta1), upper=1)$value
          return(integral-p1)
        }
        
        # Compute b20
        Fb20 <- function(b10,b20,p2) {
          b10 -> x[1]
          b20 -> x[2]
          integral<-integrate(function(v) {
            sapply(v,function(v) {
              integrate(function(u)((theta*(1-exp(-theta))*exp(-theta*(u+v)))
                                    /(exp(-theta)+exp(-theta*(u+v))-exp(-theta*u)-exp(-theta*v))^2),lower=0,
                        upper=exp(-((((-log(v))^(1/beta2))*x[2])/x[1])^beta1))$value
            })
          },
          lower= exp(-(1/x[2])^beta2), upper=1)$value
          return(integral-p2)
        }
        
        model <- function(x){
          c(Fb10(x[1],x[2],p1), Fb20(x[1],x[2],p2))
        }
        
        (sol <- multiroot(f = model, start = c(1, 1)))
        
        sol<-as.data.frame(sol[1])
        b10<-sol[1,]
        b20<-sol[2,]
        
      }
      
      ## Computation of the numerator
      # Note: Only marginal Weibull distributions for fT10, fT20, ST10, ST20.
      
      fT10 <- function(t) (beta1/b10) * ( (t/b10)^(beta1-1) ) * (exp(-(t/b10)^beta1))
      ST10 <- function(t) exp(-(t/b10)^beta1)
      
      fT20 <- function(t) (beta2/b20) * ( (t/b20)^(beta2-1) ) * (exp(-(t/b20)^beta2))
      ST20 <- function(t) exp(-(t/b20)^beta2)
      
      
      # Sstar0 and fstar0 for any copula
      Sstar0 <- function(t) Sstar0 <- Sstar(x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
      fstar0 <- function(t) fstar0 <-(-grad(Sstar,x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
      
      #########################################
      
      aux21 <- function(t,y) theta*exp(-theta*(ST10(t)+y))*(1-exp(-theta))/(exp(-theta)-exp(-theta*ST10(t))-exp(-theta*y)+exp(-theta*(ST10(t)+y)))^2
      
      aux22 <- function(u) {integrate(aux21,0, ST20(u),t=u,subdivisions=10000)$value} # t=u indicates that we are derivating respect to the other variable in aux21(t,y). That is, respect to y.
      
      lambdaC10 <- function(t) aux22(t)*fT10(t)/Sstar0(t)
      
      lambdaC11 <- function(t) HR1*lambdaC10(t)
      
      aux23 <- function(x,t) theta*exp(-theta*(x+ST20(t)))*(1-exp(-theta))/(exp(-theta)-exp(-theta*x)-exp(-theta*ST20(t))+exp(-theta*(x+ST20(t))))^2
      aux24 <- Vectorize(function(u){integrate(aux23,0,ST10(u),t=u,subdivisions=10000)$value}) #t=u indicates that we are derivating respect to the other variable in aux23(x,t). That is, respect to x.
      
      lambdaC20 <- function(t) aux24(t)*fT20(t)/Sstar0(t)
      
      lambdaC21 <- function(t) HR2*lambdaC20(t)
      
      # EVALUATION OF LambdaC20 BEFORE COMPUTATION (IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 BECAUSE IT IS NOT ALWAYS EVALUABLE AT T=0)
      LambdaC20_check <- tryCatch(LambdaC20 <- function(t) integrate(lambdaC20,lower=0,upper=t,subdivisions=10000)$value,error = function(e) e)
      
      # WHENEVER LambdaC20 FAILS, WE INCREASE THE LOWER LIMIT OF INTEGRATION
      lower_LambdaC20 <- 0
      while(inherits(LambdaC20_check, "error")=="TRUE" ){
        lower_LambdaC20 <- lower_LambdaC20 + 0.001
        LambdaC20_check <- tryCatch(LambdaC20 <- function(t) integrate(lambdaC20,lower=0+lower_LambdaC20,upper=t,subdivisions=10000)$value,error = function(e) e)
      }
      
      # Double integral directly to calculate LambdaC20.
      LambdaC20 <- Vectorize(function(t) integrate(lambdaC20,lower=0+lower_LambdaC20,upper=t,subdivisions=10000)$value)
      
      # Computation of the hazards for both groups
      Lstar0 <- function(t) lambdaC10(t) + lambdaC20(t)
      Lstar1 <- function(t) lambdaC11(t) + lambdaC21(t)
      
      # Computation of HRstar
      HRstar <- function(t) Lstar1(t)/Lstar0(t)
      logHRstar <- function(t) log(Lstar1(t)/Lstar0(t))
      
      # temp3 <- function(t) logHRstar(t)*fstar0(t)
      # https://stackoverflow.com/questions/43189512/double-integral-in-r
      temp3 <- Vectorize(function(t) logHRstar(t)*fstar0(t))
      
      # EVALUATION OF temp4 BEFORE COMPUTATION (IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 BECAUSE IT IS NOT ALWAYS EVALUABLE AT T=0)
      temp4_check <- tryCatch(temp4 <- integrate(temp3,0,1,subdivisions=10000)$value,error = function(e) e)
      
      # WHENEVER temp4 FAILS, WE INCREASE THE LOWER LIMIT OF INTEGRATION
      lower_temp4 <- 0
      while(inherits(temp4_check, "error")=="TRUE" ){
        lower_temp4 <- lower_temp4 + 0.001
        temp4_check<-tryCatch(temp4 <- integrate(temp3, lower_temp4, 1, subdivisions=10000)$value,error = function(e) e)
      }
      
      # Double integral directly
      temp4 <- integrate(temp3,0+lower_temp4, 1, subdivisions=10000)$value
      numerator <- (temp4)^2

      ## Computation of PROBT1UNC
      PROBT1UNC_temp_num <- function(t) exp(-HR2*LambdaC20(t))*Sstar0(t)*lambdaC10(t)
      PROBT1UNC_temp_den <- function(t) exp(-LambdaC20(t))*1/2 + exp(-HR2*LambdaC20(t))*1/2
      PROBT1UNC_temp <- Vectorize(function(t){PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t)})
      PROBT1UNC_int_check <- tryCatch(integrate(PROBT1UNC_temp,lower=0, upper=1,subdivisions=10000)$value, error = function(e) e)
      
      ############################################
      # WE EVALUATE THE FUNCTION PROBT1UNC_int BECAUSE IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 PROBABLY DUE TO THE LOWER LIMITS OF THE INTERGATES lambdaC20 AND temp4.
      # WHEN IT FAILS, WE SEARCH FOR THE MINIMUM EVALUABLE LIMIT OF INTEGARTION FOR lambdaC20 AND temp4; AND WE SET A LOWER LIMITS OF 0.001 FOR THE REST OF INTEGRATES
      # TO ENSURE CONVERGENCE. 
      
      lower_PROBT1UNC_int <- 0
      inc_lower <- 0
      
      while(inherits(PROBT1UNC_int_check, "error")=="TRUE"){
        lower_PROBT1UNC_int<-0.001
        inc_lower<-inc_lower+0.001
        
        aux22 <- function(u) integrate(aux21,0.001, ST20(u),t=u,subdivisions=10000)$value # t=u indicates that we are derivating respect to the other variable in aux21(t,y). That is, respect to y.
        
        lambdaC10 <- function(t) aux22(t)*fT10(t)/Sstar0(t)
        lambdaC11 <- function(t) HR1*lambdaC10(t)
        
        aux23 <- function(x,t) theta*exp(-theta*(x+ST20(t)))*(1-exp(-theta))/(exp(-theta)-exp(-theta*x)-exp(-theta*ST20(t))+exp(-theta*(x+ST20(t))))^2
        aux24 <- function(u) integrate(aux23,0.001,ST10(u),t=u,subdivisions=10000)$value #t=u indicates that we are derivating respect to the other variable in aux23(x,t). That is, respect to x.
        
        lambdaC20 <- function(t) aux24(t)*fT20(t)/Sstar0(t)
        lambdaC21 <- function(t) HR2*lambdaC20(t)
        
        LambdaC20 <- function(t) integrate(lambdaC20,lower=lower_LambdaC20 + inc_lower,upper=t,subdivisions=10000)$value
        
        Lstar0 <- function(t) lambdaC10(t) + lambdaC20(t)
        Lstar1 <- function(t) lambdaC11(t) + lambdaC21(t)
        
        HRstar <- function(t) Lstar1(t)/Lstar0(t)
        logHRstar <- function(t) log(Lstar1(t)/Lstar0(t))
        
        temp3 <- Vectorize(function(t) logHRstar(t)*fstar0(t))
        
        # Double integral
        temp4 <- integrate(temp3,lower_temp4 + inc_lower,1,subdivisions=10000)$value
        numerator <- (temp4)^2
        
        PROBT1UNC_temp_num <- function(t) exp(-HR2*LambdaC20(t))*Sstar0(t)*lambdaC10(t)
        PROBT1UNC_temp_den <- function(t) exp(-LambdaC20(t))*1/2 + exp(-HR2*LambdaC20(t))*1/2
        PROBT1UNC_temp <- Vectorize(function(t) PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t))
        PROBT1UNC_int_check<- tryCatch(integrate(PROBT1UNC_temp,lower=0.001, upper=1,subdivisions=10000)$value, error = function(e) e)
      } 

      PROBT1UNC_int <- integrate(PROBT1UNC_temp,lower=lower_PROBT1UNC_int, upper=1,subdivisions=10000)$value
      
      ############################################
      # ARE VALUE:
      AREstarT <- numerator/((log(HR1)^2) * PROBT1UNC_int * (1-Sstar0(1)))
      AREstarT 
    }
  
  return(AREstarT)
}


#######################################################################################
# Function: COMBINED_HR 
#
#######################################################################################
# Description: It computes the Combined HR(t) (HR_Star) funtion 
# Parameters:
# rho0	  Spearman's coefficient between T1 and T2 in control group
# rho1	  Spearman's coefficient between T1 and T2 in treatment group
# beta1	  Shape parameter for a Weibull law for the relevant event
# beta2   Shape parameter for a Weibull law for the additional event 
# HR1     Hazard Ratio for a Weibull law for the relevant event
# HR2     Hazard Ratio for a Weibull law for the additional event
# p1      Proportion of the relevant event expected in group zero
# p2      Proportion of the additional event expected in group zero
# case    Censoring case -- > 1 (default) or 3
# copula  Copula used:
#            Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#            Elliptical: "Normal" or "T"
#            Extreme Value: "Galambos", "HuslerReiss", "Gumbel", "Tawn" or "Tev"
#            Other: "FGM" or "Plackett"
#######################################################################################

COMBINED_HR <- function(t,rho0,rho1=rho0, beta1, beta2, HR1, HR2, p1, p2, case = 1, copula="Frank",rhoType='Spearman'){ 
  ############################################################
  ###### 0. WARNINGS AND ERRORS
  if(rho0 > 1)   stop("rho too big",call.=FALSE)
  if(rho0 < -1)  stop("rho too small",call.=FALSE)
  ############################################################
  
  ############################################################ 
  ###### 1. ELECTION OF THE COPULA
  copula0 <- CopulaSelection(copula,rho0,rhoType)
  theta <- copula0[[2]]   
  which.copula0 <- copula0[[1]]
  which.copula1 <- CopulaSelection(copula,rho1,rhoType)[[1]]  
  ############################################################
  
  ############################################################
  ###### 2. ELECTION OF THE MARGINAL DISTRIBUTIONS
  MarginSelec <- MarginalsSelection(beta1,beta2,HR1,HR2,p1,p2,case,theta,copula=copula)
  T1dist <- MarginSelec[[1]]
  T2dist <- MarginSelec[[2]]
  T1pdist <- MarginSelec[[3]]
  T2pdist <- MarginSelec[[4]]
  T10param <- MarginSelec[[5]]
  T20param <- MarginSelec[[6]]
  T11param <- MarginSelec[[7]]
  T21param <- MarginSelec[[8]]
  ############################################################
  
  ############################################################
  ###### 3. ARE EXPRESSION FOLLOWING: 
  # Gomez G, Lagakos SW. Statistical considerations when using a composite endpoint for comparing treatment groups. 
  # Stat Med. 2013; 32:719–738 (pages 2 and 3).
  
  # Bivariate distribution in control and treatment groups
  distribution0 <- mvdc(copula = which.copula0, margins = c(T1dist, T2dist),paramMargins = list(T10param, T20param))
  distribution1 <- mvdc(copula = which.copula1, margins = c(T1dist, T2dist),paramMargins = list(T11param, T21param))
  
  # Inside the integral in the numerator
  inside_integral <- function(t){
    
    Sstar0 <- Sstar(x = t,dist1 = T1pdist,dist2 = T2pdist, param1 = T10param,param2 = T20param, dist_biv = distribution0)
    Sstar1 <- Sstar(x = t,dist1 = T1pdist,dist2 = T2pdist, param1 = T11param,param2 = T21param, dist_biv = distribution1)
    
    fstar0 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T10param, param2 = T20param,dist_biv = distribution0))
    fstar1 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T11param, param2 = T21param,dist_biv = distribution1))
    
    Lstar0 <- (fstar0/Sstar0)
    Lstar1 <- (fstar1/Sstar1)
    
    HRstar <- (Lstar1/Lstar0)
    
    logHRstar <- log(HRstar)
    
    ##-- Numerical issues
    fstar0[fstar0<0] <- 0
    logHRstar[is.na(logHRstar) | logHRstar==Inf | logHRstar== -Inf] <- 0
    
    return(logHRstar*fstar0)
  }
  
  # Integral in the numerator
  integral <-integrate(inside_integral,lower=0,upper=1,subdivisions=1000,stop.on.error = FALSE) 
  numerator<-(integral$value)^2
  
  # Denominator
  Sstar0_1 <- Sstar(x=1,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0) 
  ST10_1 <- 1 - do.call(T1pdist,c(q=1,T10param)) 
  denominator <- ((log(HR1))^2)*(1-Sstar0_1)*(1-ST10_1)
  
  # ARE value
  AREstarT <- (numerator/denominator)
  
  # If the integral is not computed, we assign a missing value
  if(integral$message!="OK") {AREstarT <- NA}
  
  ############################################################ 
  ###### COMBINED HAZARD RATIO FUNCTION
  Sstar0<-Sstar(x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0) 
  Sstar1<-Sstar(x=t,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)   
  fstar0<-(-grad(Sstar,x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  fstar1<-(-grad(Sstar,x=t,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))   
  
  if (rho0==0.001 & HR1==HR2){     
    HRstar_function<-function(t) {
      HRstar <-HR1*((t+1)/(t+1)) 
      return(c(Sstar0,Sstar1,HRstar) *((t+1)/(t+1))   ) 
    } 
  } else
    if (rho0==0.001 & beta1==beta2){    
      HRstar_function<-function(t) {
        b10<-as.numeric(T10param[2])
        b20<-as.numeric(T20param[2])
        HRstar <-(HR1+HR2*(b10/b20)^beta1)/(1+(b10/b20)^beta1)
        return(c(Sstar0,Sstar1,HRstar) *((t+1)/(t+1))   ) 
      } 
    } else {
      HRstar_function <- function(t) {
        Lstar0 <- (fstar0/Sstar0)     
        Lstar1 <- (fstar1/Sstar1)                                                             
        HRstar <- (Lstar1/Lstar0)                                                                                                                   
        return(c(Sstar0,Sstar1,HRstar,Lstar0,Lstar1))
      }
    }
  
  return(HRstar_function(t))
  
}


#######################################################################################
# Function: CopulaSelection 
#
#######################################################################################
# Description: Constructs a copula class object from the family given and the
#              the corresponding dependence parameter grom the given correlation
# Parameters:
# copula  Copula given:
#            Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#            Elliptical: "Normal" or "T"
#            Extreme Value: "Galambos", "HuslerReiss", "Gumbel", "Tawn" or "Tev"
#            Other: "FGM" or "Plackett"
# rho	  Spearman's coefficient between the 2 marginal distributions
#######################################################################################

CopulaSelection <- function(copula,rho,rhoType='Spearman'){
  
  if(rhoType=='Spearman'){
    theta <- switch(copula,
                    Frank =       iRho(frankCopula(1),rho),
                    Gumbel =      iRho(gumbelCopula(2),rho),
                    Clayton =     iRho(claytonCopula(1),rho),
                    FGM =         iRho(fgmCopula(1),rho),
                    Normal =      iRho(normalCopula(0.5),rho),
                    'T' =         iRho(normalCopula(0.5),rho), # iRho(tCopula(0.5),rho), --> see Details in ?iRho
                    Galambos =    iRho(galambosCopula(0.5),rho),
                    HuslerReiss = iRho(huslerReissCopula(0.5),rho),
                    Tawn =        iRho(tawnCopula(0.5),rho),
                    Tev =         iRho(tevCopula(0.5),rho),
                    Plackett =    iRho(plackettCopula(0.5),rho))
  }else{
    theta <- switch(copula,
                    Frank =       iTau(frankCopula(1),rho),
                    Gumbel =      iTau(gumbelCopula(2),rho),
                    Clayton =     iTau(claytonCopula(1),rho),
                    FGM =         iTau(fgmCopula(1),rho),
                    Normal =      iTau(normalCopula(0.5),rho),
                    'T' =         iTau(normalCopula(0.5),rho), # iRho(tCopula(0.5),rho), --> see Details in ?iRho
                    Galambos =    iTau(galambosCopula(0.5),rho),
                    HuslerReiss = iTau(huslerReissCopula(0.5),rho),
                    Tawn =        iTau(tawnCopula(0.5),rho),
                    Tev =         iTau(tevCopula(0.5),rho),
                    Plackett =    iTau(plackettCopula(0.5),rho))
  }
  

  which.copula <- switch(copula,
                         Frank =       archmCopula(family = "frank", dim = 2, param = theta),
                         Gumbel =      archmCopula(family = "gumbel", dim = 2, param = theta),
                         Clayton =     archmCopula(family = "clayton", dim = 2, param = theta),
                         FGM =         fgmCopula(dim = 2, param = theta),
                         Normal =      normalCopula(dim = 2, param = theta),
                         'T' =         tCopula(dim = 2, param = theta),
                         Galambos =    galambosCopula(param = theta),
                         HuslerReiss = huslerReissCopula(param = theta),
                         Tawn =        tawnCopula(param = theta),
                         Tev =         tevCopula(param = theta),
                         Plackett =    plackettCopula(param = theta ))
  
  

  return(c(which.copula,theta))
}


#######################################################################################
# Functions: several copula densities
# dFrank: Frank copula
# dClayton: Clayton copula
# dGumbel: Gumbel copula
# dFGM: FGM copula
# dNormal: Normal copula
#
# Reference: Table 3.3 of Trivedi, Pravin K., and David M. Zimmer. Copula Modeling: An Introduction for
#            Practitioners. Foundations and Trends in Econometrics, 1 (2007), 1–111.
#
#######################################################################################
# Description: It computes the Combined HR*(t) (HR_Star) function 
# Parameters:
# u	      value between 0 and 1 for endpoint 1 u=F(t_1)
# v	      value between 0 and 1 for endpoint 2 v=F(t_2)
# theta	  Copula parameter
#######################################################################################

dFrank <- function(u,v,theta){(theta*(1-exp(-theta))*exp(-theta*(u+v)))/ (exp(-theta) +  exp(-theta*(u+v))- exp(-theta*u)-exp(-theta*v))^2}

dClayton <- function(u,v,theta){ (u*v)^(-theta-1) * (theta+1) * (u^(-theta) + v^(-theta) - 1)^(-2 - 1/theta)}

dGumbel <- function(u,v,theta){
  u1 <- -log(u)
  u2 <- -log(v)
  num1 <- exp(-(u1^theta + u2^theta)^(1/theta))
  num2 <- (u*v)^(-1)
  num3 <- (u1*u2)^(theta-1)
  num4 <- (u1^theta + u2^theta)^(1/theta) + theta - 1
  num <- num1 * num2 * num3 * num4
  den <- (u1^theta + u2^theta)^(2-1/theta)
  num/den
}

dFGM <- function(u,v,theta) 1 + theta * (1- 2*u) * (1 - 2*v)

dNormal <- function(u,v,theta){
  x <- qnorm(u)
  y <- qnorm(v)
  
  (1 - theta^2)^(-1/2) * exp(-(x^2 + y^2 - 2*theta*x*y)/(2*(1-theta^2))) * exp((x^2 + y^2)/2)
}


#######################################################################################
# Function: MarginalsSelection 
#
#######################################################################################
# Description: Returns the family distribution and parameters of the marginals
#              (ONLY WEIBULL DISTRIBUTIONS SO FAR) 
# Parameters:
# beta1	  Shape parameter for a Weibull law for the relevant event
# beta2   Shape parameter for a Weibull law for the additional event 
# HR1     Hazard Ratio for a Weibull law for the relevant event
# HR2     Hazard Ratio for a Weibull law for the additional event
# p1      Proportion of the relevant event expected in group zero
# p2      Proportion of the additional event expected in group zero
# case    Censoring case -- > 1 (default) or 3
# theta   Dependence parameter for the bivariate distribution in control group
#######################################################################################

MarginalsSelection <- function(beta1,beta2,HR1,HR2,p1,p2,case,theta,copula='Frank') 
{
  # Scale parameters for group 0 b10,b20
  
  dcopula <- get(paste0('d',copula))

  ## -- Case 1 -----------------------
  if(case==1) {
    b10 <- 1/((-log(1-p1))^(1/beta1))
    b20 <- 1/((-log(1-p2))^(1/beta2))
  
  ## -- Case 2 -----------------------  
  } else if (case==2) {  
    Fb10 <- function(b10,p1){
      integral<-integrate(function(u) {
        sapply(u, function(u) {
          integrate(function(v){dcopula(u,v,theta)}, lower=1-exp((b10*(-log(1-u))^(1/beta1))^beta2*log(1-p2)), upper=1)$value
        })
      }, lower=0 , upper=1-exp(-1/b10^beta1))$value
      return(integral-p1) 
    }
    limits <- c(0.00001,10000)                                         # The first and the last values must be in opposite signs for the function
    b10 <- try(uniroot(Fb10, interval=limits,p1=p1)$root,silent=TRUE)  # Find the root (value which equals the function zero)
    b20 <- 1/(-log(1-p2))^(1/beta2)
    if(class(b10)=='try-error'){
      dcopula <- dFrank
      b10 <- uniroot(Fb10, interval=limits,p1=p1)$root
      dcopula <- get(paste0('d',copula))
      limits <- c(0.8,1.2)*b10 
      b10 <- uniroot(Fb10, interval=limits,p1=p1)$root
    }
    
  ## -- Case 3 -----------------------  
  } else if (case==3) {
    b10 <- 1/((-log(1-p1))^(1/beta1))
    Fb20 <- function(b20,p2) {
      integral<-integrate(function(v) {
        sapply(v,function(v) { 
          integrate(function(u){dcopula(u,v,theta)},lower=1-exp((b20*(-log(1-v))^(1/beta2))^beta1*log(1-p1)),upper=1)$value
        })
      }, 
      lower=0 , upper=1-exp(-1/b20^beta2))$value
      return(integral-p2)
    }
    limits <- c(0.00001,10000) 
    b20 <- try(uniroot(Fb20, interval=limits,p2=p2)$root,silent=TRUE)
    if(class(b20)=='try-error'){
      dcopula <- dFrank
      b20 <- uniroot(Fb20, interval=limits,p2=p2)$root
      dcopula <- get(paste0('d',copula))
      limits <- c(0.8,1.2)*b20 
      b20 <- uniroot(Fb20, interval=limits,p2=p2)$root
    }
    
  ## -- Case 4 -----------------------
  } else if (case==4) {
    
    # We need to create x[1] and x[1] (we assign NA's)
    x <- c(NA,NA)
    
    # We need to change the name of variables as (b10=x[1],b20=[2]) to execute 'multiroot'
    # Compute b10
    Fb10 <- function(b10,b20,p1){
      b10-> x[1]
      b20-> x[2]
      integral<-integrate(function(u) {
        sapply(u, function(u) {
          integrate(function(v){dcopula(u,v,theta)}, lower=1-exp(-(x[1]*(-log(1-u))^(1/beta1)/x[2])^beta2), upper=1)$value
        })
      }, lower= 0 , upper=1-exp(-1/x[1]^beta1))$value
      return(integral-p1)
    }
    
    # Compute b20
    Fb20<-function(b10,b20,p2) {
      b10-> x[1]
      b20-> x[2]
      integral<-integrate(function(v) {
        sapply(v,function(v) {
          integrate(function(u){dcopula(u,v,theta)},lower=1-exp(-(x[2]*(-log(1-v))^(1/beta2)/x[1])^beta1), upper=1)$value
        })
      },
      lower= 0, upper=1-exp(-1/x[2]^beta2))$value
      return(integral-p2)
    }
    
    model <- function(x){
      c(Fb10(x[1],x[2],p1), Fb20(x[1],x[2],p2))
    }
    
    (sol <- multiroot(f = model, start = c(1, 1)))
    
    sol<-as.data.frame(sol[1])
    b10<-sol[1,]
    b20<-sol[2,]
    
  }
  
  # Scale parameters for group 1 b11,b21 (Although we do not need to compute the scale parameters for group 1 (b11, b21) to calculate the ARE)
  b11 <- b10/HR1^(1/beta1)
  b21 <- b20/HR2^(1/beta2)
  
  # Probabilities p11,p21 (Although we do not need to calculate the ARE)
  p11 <- 1-exp(-(1/b11)^beta1)  # Need to be reviewed taking into account the competing risk
  p21 <- 1-exp(-(1/b21)^beta2)  # Need to be reviewed taking into account the competing risk
  
  T1dist <- "weibull"
  T2dist <- "weibull"
  
  T1pdist <- pweibull
  T2pdist <- pweibull
  
  T10param <- list(shape = beta1, scale = b10)
  T20param <- list(shape = beta2, scale = b20)
  
  T11param <- list(shape = beta1, scale = b11)
  T21param <- list(shape = beta2, scale = b21)
  
  return(list(T1dist,T2dist,T1pdist,T2pdist,T10param,T20param,T11param,T21param,p11,p21))	
}


#######################################################################################
# Function: Sstar 
#
#######################################################################################
# Description: Returns the value of the survival function of S* at point x given the
#              marginal distributions and the bivariate distributions via copula 
# Parameters:
# x	        Point in which to be evaluated
# dist1     Distribution function of the marginal T1 (pweibull) 
# dist2     Distribution function of the marginal T2 (pweibull) 
# param1    Parameters of the marginal distribution function T1 (pweibull) 
# param2    Parameters of the marginal distribution function T2 (pweibull) 
# dist_biv  Distribution function of the bivariate distribution via copula
#######################################################################################

Sstar<-function(x,dist1,dist2,param1,param2,dist_biv) { 
  y <- if(length(x) == 1) c(x,x) else cbind(x,x)
  return(
    1
    - do.call(dist1,c(list(q=x),param1))
    - do.call(dist2,c(list(q=x),param2))
    + (pMvdc(y, dist_biv))
  )
}


#######################################################################################
# Function: schoendfeld.formula 
#
#######################################################################################
# Description: Returns the NUMBER OF EVENTS according to Schoenfeld formula
# Parameters:
# alpha     Prob(Error type I)
# power     1 - Prob(Error type II)
# HR        Hazard Ratio
#######################################################################################

schoendfeld.formula <- function(alpha,power,HR) E <- 4*(qnorm(1-alpha) +  qnorm(power))^2 / (log(HR))^2


#######################################################################################
# Function: freedman.formula 
#
#######################################################################################
# Description: Returns the NUMBER OF EVENTS according to Freedman formula
# Parameters:
# alpha     Prob(Error type I)
# power     1 - Prob(Error type II)
# HR        Hazard Ratio
#######################################################################################

freedman.formula <- function(alpha,power,HR) E <- (HR+1)^2 * (qnorm(1-alpha) +  qnorm(power))^2/(HR-1)^2


#######################################################################################
# Function: trapezoidal.integration
#
#######################################################################################
# Description: It calculates an integral using trapezoides
# Parameters:
# x	      Variable of integration
# f       Image of the function for each point x
#######################################################################################

trapezoidal.integration = function(x, f)
{
  ### 3 checks to ensure that the arguments are numeric and of equal lengths
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable of integration "x" is not numeric.')
  }
  
  # check if the integrand is numeric
  if (!is.numeric(f))
  {
    stop('The integrand "f" is not numeric.')
  }
  
  # check if the variable of integration and the integrand have equal lengths 
  # WARNING AVOIDED BECAUSE IT FAILS WHEN USING THIS FUNCTION INSIDE ANOTHER FUNTION  
  ##if (length(x) != length(f))
  ##{
  ##       stop('The lengths of the variable of integration and the integrand do not match.')
  ##}
  
  ### finish checks
  
  # obtain length of variable of integration and integrand
  n = length(x)
  
  # integrate using the trapezoidal rule
  integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))
  
  # print the definite integral
  return(integral)
}


#######################################################################################
# Function: sample_size_single_events 
#
#######################################################################################
# Description: It computes the sample size for each endpoint 
#              (taking into account Cause Specific Hazard Ratios)
# Parameters:
# rho0	                Spearman's coefficient between T1 and T2 in control group
# rho1	                Spearman's coefficient between T1 and T2 in treatment group
# beta1	                Shape parameter for a Weibull law for the relevant event
# beta2                 Shape parameter for a Weibull law for the additional event 
# HR1                   Hazard Ratio for a Weibull law for the relevant event
# HR2                   Hazard Ratio for a Weibull law for the additional event
# p1                    Proportion of the relevant event expected in group zero
# p2                    Proportion of the additional event expected in group zero
# case                  Censoring case -- > 1 (default) or 3
# copula                Copula used:
#                          Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#                          Elliptical: "Normal"
#                          Other: "FGM"
# alpha_error           Probability of type I error
# power                 Power for test
# formula_sample_size   Formulae to perform sample size calculations
#######################################################################################

sample_size_single_events <- function(rho0,rho1=rho0, beta1, beta2, HR1, HR2, p1, p2, case, copula, alpha_error, power, formula_sample_size,rhoType='Spearman'){
  
  
  ##-- Copula selection
  # copula0 <- CopulaSelection(copula,rho0)
  # theta <- copula0[[2]]   
  # which.copula0 <- copula0[[1]]
  # which.copula1 <- CopulaSelection(copula,rho1)[[1]]
  theta <- CopulaSelection(copula,rho0,rhoType)[[2]]
  
  ##-- Computation of probabilities of observed events in treated arm for both endpoints
  p_treated <- unlist(MarginalsSelection(beta1=beta1,beta2=beta2,HR1=HR1,HR2=HR2,p1=p1,p2=p2,case=case,theta=theta)[9:10])
  
  ##-- Number of events
  if (formula_sample_size == "Schoenfeld"){
    E1 <- schoendfeld.formula(alpha_error,power,HR1)   # Endpoint 1
    E2 <- schoendfeld.formula(alpha_error,power,HR2)   # Endpoint 2
  }
  if (formula_sample_size == "Freedman"){
    E1 <- freedman.formula(alpha_error,power,HR1)      # Endpoint 1
    E2 <- freedman.formula(alpha_error,power,HR2)      # Endpoint 2
  }
  
  ##-- Survival proportions
  # Endpoint 1
  pi10 <- 1 - p1                                       # probability of survial for endpoint 1 in control arm       
  pi11 <- 1 - p_treated[1]                             # probability of survial for endpoint 1 in treated arm
  
  # Endpoint 2
  pi20 <- 1 - p2                                       # probability of survial for endpoint 2 in control arm    
  pi21 <- 1 - p_treated[2]                             # probability of survial for endpoint 2 in treated arm
  
  
  
  ##-- Sample size
  # Endpoint 1
  N_TOT1 <- (2*E1)/(2-pi10-pi11)                       # TOTAL SAMPLE SIZE
  n1 <- N_TOT1/2                                       # PATIENTS PER GROUP NEEDED
  
  # Endpoint 2
  N_TOT2 <- (2*E2)/(2-pi20-pi21)                       # TOTAL SAMPLE SIZE
  n2 <- N_TOT2/2                                       # PATIENTS PER GROUP NEEDED
  
  res <- c(N_TOT1,N_TOT2,E1,E2)
  return(res)
  
}


#######################################################################################
# Function: Different_scenarios 
#
#######################################################################################
# Description: It computes ARE for several combinations of parameters
# Parameters:
# rho0	                Spearman's coefficient between T1 and T2 in control group
# rho1	                Spearman's coefficient between T1 and T2 in treatment group
# beta1	                Shape parameter for a Weibull law for the relevant event
# beta2                 Shape parameter for a Weibull law for the additional event 
# HR1                   Hazard Ratio for a Weibull law for the relevant event
# HR2                   Hazard Ratio for a Weibull law for the additional event
# p1                    Proportion of the relevant event expected in group zero
# p2                    Proportion of the additional event expected in group zero
# case                  Censoring case (1,2,3 or 4)
# copula                Copula used:
#                          Archimedean: "Frank", "Gumbel" or "Clayton"
#                          Elliptical: "Normal"
#                          Other: "FGM"
# rho_thype             Type of correlation (Spearman or Kendall)
#######################################################################################

Different_scenarios <- function(rho0,rho1=rho0, beta1, beta2, HR1, HR2, p1, p2, case = case, copula = copula, rhoType='Spearman'){

  beta1V <- beta1
  beta2V <- beta2
  
  p1V <- p1
  p2V <- p2
  HR1V <- HR1
  
  rhoV <- c(0.01, 0.15, 0.3, 0.5, 0.7, 0.9, 0.95) 

  # IF THERE IS SMALL DISTANCE BETWEEN HR'S, THEN THE PLOT CENTERS THE 5 CURVES IN THE HR_RE VALUE. OTHERWISE, MAXIMUM AND MINIMUM VALUES ARE PLOTTED (HR_RE AND HR_AE), WITH THE CORRESPONDING 3 CURVES INBETWEEN. 
  abs_dif_HR1_HR2 <- abs(HR1-HR2)
  
  if(abs_dif_HR1_HR2 <  0.1) {
    HR2V <- c(HR2-0.1,HR2-0.05,HR2,HR2+0.05, HR2+0.1)
  }
  
  if( (abs_dif_HR1_HR2>= 0.1) & (HR1>HR2)  ) {
    HR2V <- c(HR2,HR2+(abs_dif_HR1_HR2/4),HR2+2*(abs_dif_HR1_HR2/4), HR2+3*(abs_dif_HR1_HR2/4),HR1)
  }
  
  if((abs_dif_HR1_HR2>= 0.1) & (HR1<HR2)){
    HR2V <- c(HR1,HR2-3*(abs_dif_HR1_HR2/4),HR2-2*(abs_dif_HR1_HR2/4), HR2-(abs_dif_HR1_HR2/4),HR2)
  }

  # CREATE THE DATASET
  dataset <- expand.grid(beta1=beta1V,beta2=beta2V,p1=p1V,p2=p2V,HR1=HR1V,HR2=HR2V,rho=rhoV)
  dataset$ARE <- sfApply(dataset[,1:7],1,ARE.ARRAY,case=case,copula=copula,rhoType=rhoType)

  return(dataset)
}

