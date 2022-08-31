

cr.interval <- function(data=NULL,covs=NULL,method="exp",start1=NULL,start2=NULL,k=10,gam=c(10,10),maxt=NULL,maxit=30,abstol=1e-4,pred=FALSE){
  
  library(optimParallel)  
  
  data <- data
  method <- method
  
  library(tibble)
  if(is.data.frame(data)==FALSE & is_tibble(data)==FALSE) stop("Please enter a valid data frame")
  
  if (!(method %in% c("exp","wei","bs"))) stop("Please enter a valid baseline hazards specification. Use 'exp' for the exponential distribution, 'wei' for the Weibull distribution,  or 'bs' for the flexible distribution using B-splines")
  
  gam=gam
  
  startval1 <- start1
  startval2 <- start2
  k <- k
  degree<-3
  covs <- covs
  gam <- gam
  maxit <- maxit
  abstol <- abstol
  
  if(is.null(startval1) | is.null(startval2)) stop("Please enter start values for events 1 and 2 respectively.")
  
  
  if(method=="bs" & length(startval1)!=k) warning("Length of each set of start values should equal k for method 'bs'. Only the first value will be used.")
  if(method=="bs" & length(startval2)!=k) warning("Length of each set of start values should equal k for method 'bs'. Only the first value will be used.")
  if(method=="exp" & length(startval1)!=1) warning("Only one start value must be specified for event 1 for method 'exp'. Only the first value will be used.")    
  if(method=="exp" & length(startval2)!=1) warning("Only one start value must be specified for event 2 for method 'exp'. Only the first value will be used.")
  if(method=="wei" & length(startval1)!=2) warning("Only two start value must be specified for event 1 for method 'wei'. Only the first value will be used.")
  if(method=="wei" & length(startval2)!=2) warning("Only two start value must be specified for event 2 for method 'wei'. Only the first value will be used.")
  
  hes<-TRUE
  covs<-covs
  
  method<-method
  
  mstatdat <- data
  
  
  ###B-spline basis functions for baseline splines
  
  
  if (method=="bs"){
    
    # Splines -----------------------------
    library(splines2)
    library(fda)
    
    #degree<-3
    
    bs.smoother<-function(X){
      #  k=if(is.null(nk)) {length(knots)} else nk
      #  df = k+degree
      #  deriv=2
      
      B <- list()
      #  probs <- (1/k)*(1:(k-1))
      #  if(is.null(knots)){knots=quantile(X,probs=probs,na.rm=T)} else {knots=knots}
      
      #k=10 #degrees of freedom
      
      library(mgcv)
      smoother <-  smooth.construct(s(x,k=k,  bs = "bs", m=c(degree,(degree-1))), data.frame(x=X),NULL)
      
      #B$penmat <- output <- matrix(unlist(smoother$S), ncol = k, byrow = TRUE)
      
      B$knots <- smoother$knots[-(1:3)]
      #B$knots <- smoother$knots[]
      B$knots[1] <- 0
      #B$knots[1] <- max(0,B$knots[1])
      nk=length(B$knots)
      
      B$bs   <- splines2::bSpline(X, knots=B$knots[-c(1,nk)], degree = degree, intercept = TRUE,Boundary.knots = c(B$knots[1],B$knots[nk]))
      
      B$bsint <- splines2::ibs(X, knots=B$knots[-c(1,nk)], degree = degree, intercept = TRUE,Boundary.knots = c(B$knots[1],B$knots[nk]))
      
      
      basisfd <- fda::create.bspline.basis(rangeval=c(B$knots[1],B$knots[nk]),breaks = B$knots,norder=(degree+1))
      B$penmat<- fda::bsplinepen(basisfd, Lfd = 2)
      
      B$nk = nk
      
      return(B)
    }
    
    
    
    
    t1=unique(sort(c(0,mstatdat$L1time[mstatdat$state==1], mstatdat$L1time[mstatdat$state==1]-0.5*mstatdat$R0time[mstatdat$state==1], mstatdat$L1time[mstatdat$state==1]-mstatdat$R0time[mstatdat$state==1],mstatdat$L1time[mstatdat$state==2 & mstatdat$overlap==0], mstatdat$L1time[mstatdat$state==2 & mstatdat$overlap==0]-0.5*mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==0], mstatdat$L1time[mstatdat$state==2 & mstatdat$overlap==0]-mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==0], (mstatdat$L1time[mstatdat$state==2 & mstatdat$overlap==0]+mstatdat$R1time[mstatdat$state==2 & mstatdat$overlap==0])/2, (mstatdat$L1time[mstatdat$state==2 & mstatdat$overlap==0]+mstatdat$R1time[mstatdat$state==2 & mstatdat$overlap==0]-mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==0])/2, (mstatdat$L1time[mstatdat$state==2 & mstatdat$overlap==0]+mstatdat$R1time[mstatdat$state==2 & mstatdat$overlap==0])/2 -mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==0], mstatdat$R1time[mstatdat$state==2 & mstatdat$overlap==0], mstatdat$R1time[mstatdat$state==2 & mstatdat$overlap==0]-0.5*mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==0], mstatdat$R1time[mstatdat$state==2 & mstatdat$overlap==0]-mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==0], mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==1]/4, mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==1]/2, mstatdat$R0time[mstatdat$state==2 & mstatdat$overlap==1])))
    
    t2=unique(sort(c(0,mstatdat$L1time[mstatdat$state==1], mstatdat$L1time[mstatdat$state==1]-0.5*mstatdat$R0time[mstatdat$state==1], mstatdat$L1time[mstatdat$state==1]-mstatdat$R0time[mstatdat$state==1],mstatdat$L1time[mstatdat$state==3], mstatdat$L1time[mstatdat$state==3]-0.5*mstatdat$R0time[mstatdat$state==3], mstatdat$L1time[mstatdat$state==3]-mstatdat$R0time[mstatdat$state==3], (mstatdat$L1time[mstatdat$state==3]+mstatdat$R1time[mstatdat$state==3])/2, (mstatdat$L1time[mstatdat$state==3]+mstatdat$R1time[mstatdat$state==3]-mstatdat$R0time[mstatdat$state==3])/2, (mstatdat$L1time[mstatdat$state==3]+mstatdat$R1time[mstatdat$state==3])/2 -mstatdat$R0time[mstatdat$state==3], mstatdat$R1time[mstatdat$state==3], mstatdat$R1time[mstatdat$state==3]-0.5*mstatdat$R0time[mstatdat$state==3], mstatdat$R1time[mstatdat$state==3]-mstatdat$R0time[mstatdat$state==3])))
    
    #    t1<-unique(sort(c(0,mstatdat$R0time[mstatdat$state==2],mstatdat$L1time[mstatdat$state==2],mstatdat$R1time[mstatdat$state==2])))
    #    t2<-unique(sort(c(0,mstatdat$R0time[mstatdat$state==3],mstatdat$L1time[mstatdat$state==3],mstatdat$R1time[mstatdat$state==3])))
    
    
    t1 = unique(c(t1))
    t2 = unique(c(t2))
    
    #t1 = unique(sort(c(t1,t2)))
    #t2 = unique(sort(c(t1,t2)))
    
    
    bsmooth1<-bs.smoother(X=t1)
    nk<-bsmooth1$nk
    
    smooth1<-bsmooth1$bsint
    npar1<-ncol(smooth1[,1:(nk-1)])
    hsmooth1<-bsmooth1$bs
    #pen1<-bsmooth1$penmat[,1:(nk-1)]
    #pen1<-bsmooth1$penmat[1:(k-1),1:(k-1)]
    pen1<-bsmooth1$penmat[1:(nk-1),1:(nk-1)]
    
    bsmooth2<-bs.smoother(X=t2)
    smooth2<-bsmooth2$bsint
    npar2<-ncol(smooth2[,1:(nk-1)])
    hsmooth2<-bsmooth2$bs
    #pen2<-bsmooth1$penmat[,1:(nk-1)]
    #pen2<-bsmooth2$penmat[1:(k-1),1:(k-1)]
    pen2<-bsmooth2$penmat[1:(nk-1),1:(nk-1)]
    
    
    
    #boundary.knots=c(0,max(t2)
    
    npar<-npar1+npar2
    
    # Now including covariates ------------------------------------------------
    
    ncov1<-length(covs)
    ncov2<-length(covs) #(covariates for each transition)
    ncov<-ncov1+ncov2
    nparcov<-npar+ncov
    
    if (is.null(covs)){
      
      # No covariates  -----------------------------------------------------
      if (length(startval1)!=npar1 | length(startval2)!=npar2){
        p0 = c(rep(startval1[1],npar1),rep(startval2[1],npar2)) #start values
      } else {p0=c(startval1,startval2)}
      
      gam0 = gam   
      
      
      
      # Functions ---------------------------------------------------------------
      
      negloglik.base<-function(par,gam){
        
        Pen = matrix(0,nrow=npar,ncol=npar)
        Pen[1:npar1,1:npar1]=pen1*exp(gam[1])
        Pen[(npar1+1):npar,(npar1+1):npar]=pen2*exp(gam[2])
        
        mstatdat1 = mstatdat[mstatdat$state==1,]
        mstatdat3 =  mstatdat[mstatdat$state==3,]
        mstatdat2 =  mstatdat[mstatdat$state==2 & mstatdat$overlap==0,] 
        mstatdat2o = mstatdat[mstatdat$state==2 & mstatdat$overlap==1,] 
        
        if (nrow(mstatdat1)==0){L1=1} else{   L1= ((1/6)*(exp((-predict(smooth1,mstatdat1[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat1[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))+4*exp((-predict(smooth1,mstatdat1[["L1time"]]-mstatdat1[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat1[["L1time"]]-mstatdat1[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))+exp((-predict(smooth1,mstatdat1[["L1time"]]-mstatdat1[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat1[["L1time"]]-mstatdat1[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))))
        }
        
        L3= ((mstatdat3[["R1time"]]-mstatdat3[["L1time"]])/36)*(
          exp((-predict(smooth1,mstatdat3[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat3[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,mstatdat3[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            4*exp((-predict(smooth1,(mstatdat3[["L1time"]]+mstatdat3[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,(mstatdat3[["L1time"]]+mstatdat3[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,(mstatdat3[["L1time"]]+mstatdat3[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            exp((-predict(smooth1,mstatdat3[["R1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat3[["R1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,mstatdat3[["R1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+ 
            4*exp((-predict(smooth1,mstatdat3[["L1time"]]-mstatdat3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat3[["L1time"]]-mstatdat3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,mstatdat3[["L1time"]]-mstatdat3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            16*exp((-predict(smooth1,(mstatdat3[["L1time"]]+mstatdat3[["R1time"]]-mstatdat3[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,(mstatdat3[["L1time"]]+mstatdat3[["R1time"]]-mstatdat3[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,(mstatdat3[["L1time"]]+mstatdat3[["R1time"]]-mstatdat3[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            4*exp((-predict(smooth1,mstatdat3[["R1time"]]-mstatdat3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat3[["R1time"]]-mstatdat3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,mstatdat3[["R1time"]]-mstatdat3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            exp((-predict(smooth1,mstatdat3[["L1time"]]-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat3[["L1time"]]-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,mstatdat3[["L1time"]]-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            4*exp((-predict(smooth1,((mstatdat3[["L1time"]]+mstatdat3[["R1time"]])/2)-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,((mstatdat3[["L1time"]]+mstatdat3[["R1time"]])/2)-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,((mstatdat3[["L1time"]]+mstatdat3[["R1time"]])/2)-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))+
            exp((-predict(smooth1,mstatdat3[["R1time"]]-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat3[["R1time"]]-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth2,mstatdat3[["R1time"]]-mstatdat3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))
        
        L2 = ((mstatdat2[["R1time"]]-mstatdat2[["L1time"]])/36)*(
          exp((-predict(smooth1,mstatdat2[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            4*exp((-predict(smooth1,(mstatdat2[["L1time"]]+mstatdat2[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,(mstatdat2[["L1time"]]+mstatdat2[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,(mstatdat2[["L1time"]]+mstatdat2[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            exp((-predict(smooth1,mstatdat2[["R1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2[["R1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2[["R1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))+ 
            4*exp((-predict(smooth1,mstatdat2[["L1time"]]-mstatdat2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2[["L1time"]]-mstatdat2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2[["L1time"]]-mstatdat2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            16*exp((-predict(smooth1,(mstatdat2[["L1time"]]+mstatdat2[["R1time"]]-mstatdat2[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,(mstatdat2[["L1time"]]+mstatdat2[["R1time"]]-mstatdat2[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,(mstatdat2[["L1time"]]+mstatdat2[["R1time"]]-mstatdat2[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            4*exp((-predict(smooth1,mstatdat2[["R1time"]]-mstatdat2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2[["R1time"]]-mstatdat2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2[["R1time"]]-mstatdat2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            exp((-predict(smooth1,mstatdat2[["L1time"]]-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2[["L1time"]]-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2[["L1time"]]-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            4*exp((-predict(smooth1,((mstatdat2[["L1time"]]+mstatdat2[["R1time"]])/2)-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,((mstatdat2[["L1time"]]+mstatdat2[["R1time"]])/2)-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,((mstatdat2[["L1time"]]+mstatdat2[["R1time"]])/2)-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            exp((-predict(smooth1,mstatdat2[["R1time"]]-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2[["R1time"]]-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2[["R1time"]]-mstatdat2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))
        
        if (nrow(mstatdat2o)==0){L2o=1} else{         L2o = (mstatdat2o[["R0time"]]/36)*(
          exp((-predict(smooth1,mstatdat2o[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2o[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2o[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))+   
            6*exp((-predict(smooth1,mstatdat2o[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2o[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2o[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))+
            8*exp((-predict(smooth1,mstatdat2o[["R0time"]]/4)[,1:(nk-1)]%*%((exp(par[1:npar1]))))-predict(smooth2,mstatdat2o[["R0time"]]/4)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar]))))*predict(hsmooth1,mstatdat2o[["R0time"]]/4)[,1:(nk-1)]%*%((exp(par[1:npar1])))
        )}
        
        loglik=sum(log(L1))+sum(log(L3))+sum(log(L2))+sum(log(L2o))-0.5*(t(exp(par[(1):npar]))%*%Pen%*%(exp(par[(1):npar])))
        
        return(-loglik)
      }
      
      fbase<-function(x){
        L1= ((1/6)*(exp((-predict(smooth1,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))+4*exp((-predict(smooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))+exp((-predict(smooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))))
        
        L3= ((mstatdat[["R1time"]]-mstatdat[["L1time"]])/36)*(
          exp((-predict(smooth1,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            4*exp((-predict(smooth1,(mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,(mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,(mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            exp((-predict(smooth1,mstatdat[["R1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,mstatdat[["R1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+ 
            4*exp((-predict(smooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            16*exp((-predict(smooth1,(mstatdat[["L1time"]]+mstatdat[["R1time"]]-mstatdat[["R0time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,(mstatdat[["L1time"]]+mstatdat[["R1time"]]-mstatdat[["R0time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,(mstatdat[["L1time"]]+mstatdat[["R1time"]]-mstatdat[["R0time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            4*exp((-predict(smooth1,mstatdat[["R1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,mstatdat[["R1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            exp((-predict(smooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            4*exp((-predict(smooth1,((mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,((mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,((mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))+
            exp((-predict(smooth1,mstatdat[["R1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth2,mstatdat[["R1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))
        
        L2 = ((mstatdat[["R1time"]]-mstatdat[["L1time"]])/36)*(
          exp((-predict(smooth1,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1])))+
            4*exp((-predict(smooth1,(mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,(mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,(mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1])))+
            exp((-predict(smooth1,mstatdat[["R1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["R1time"]])[,1:(nk-1)]%*%(((x[1:npar1])))+ 
            4*exp((-predict(smooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1])))+
            16*exp((-predict(smooth1,(mstatdat[["L1time"]]+mstatdat[["R1time"]]-mstatdat[["R0time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,(mstatdat[["L1time"]]+mstatdat[["R1time"]]-mstatdat[["R0time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,(mstatdat[["L1time"]]+mstatdat[["R1time"]]-mstatdat[["R0time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1])))+
            4*exp((-predict(smooth1,mstatdat[["R1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["R1time"]]-mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1])))+
            exp((-predict(smooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["L1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))+
            4*exp((-predict(smooth1,((mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,((mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,((mstatdat[["L1time"]]+mstatdat[["R1time"]])/2)-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))+
            exp((-predict(smooth1,mstatdat[["R1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["R1time"]]-mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))
        
        L2o = (mstatdat[["R0time"]]/36)*(
          exp((-predict(smooth1,mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))+   
            6*exp((-predict(smooth1,mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1])))+
            8*exp((-predict(smooth1,mstatdat[["R0time"]]/4)[,1:(nk-1)]%*%(((x[1:npar1]))))-predict(smooth2,mstatdat[["R0time"]]/4)[,1:(nk-1)]%*%(((x[(npar1+1):npar]))))*predict(hsmooth1,mstatdat[["R0time"]]/4)[,1:(nk-1)]%*%(((x[1:npar1])))
        )
        
        loglik=sum(log((L1^(I(mstatdat[["state"]]==1)) )*  
                         (L3^(I(mstatdat[["state"]]==3)))*
                         (L2^(I(mstatdat[["state"]]==2)*I(mstatdat[["overlap"]]==0)))*
                         (L2o^(I(mstatdat[["state"]]==2)*I(mstatdat[["overlap"]]==1)))))
        
        return(loglik[[1]])
      }
      
      library(pracma)
      jacob = function(par){
        pracma::jacobian(fbase,x0=par,h=.Machine$double.eps^(1/3))
      }
      
      
      # Prelim algorithm:
      count    <- 0
      Max      <- maxit
      abstol   <- abstol
      
      jj = jacob(c(exp(p0[1:npar])))
      S = colSums(jj)
      I = (t(jj)%*%jj)/nrow(mstatdat)
      
      sqrtm = function (A, kmax = 20, tol = .Machine$double.eps^(1/2)) 
      {
        stopifnot(is.numeric(A), is.matrix(A))
        if (nrow(A) != ncol(A)) 
          stop("Matrix 'A' must be square.")
        P0 <- A
        Q0 <- diag(nrow(A))
        k <- 0
        while (norm(A - P0 %*% P0, "F") > tol && k < kmax) {
          P1 <- 0.5 * (P0 + MASS::ginv(Q0, tol = tol))
          Q1 <- 0.5 * (Q0 + MASS::ginv(P0, tol = tol))
          P0 <- P1
          Q0 <- Q1
          k <- k + 1
        }
        if (k >= kmax)
          k <- -1
        return(list(B = P0, Binv = Q0, k = k, acc = norm(A - P0 %*% 
                                                           P0, "F")))
      }
      
      
      UBRE.g <- function(gam,p0){
        Pen = matrix(0,nrow=npar,ncol=npar)
        Pen[1:npar1,1:npar1]=pen1*exp(gam[1])
        Pen[(npar1+1):npar,(npar1+1):npar]=pen2*exp(gam[2])
        
        p0up = c(exp(p0[1:npar]))
        
        Sp = S-Pen%*%p0up
        Ip = I+Pen
        sq = sqrtm(I,tol=1e-40)
        sqrt.I = sq$B
        inv.sqrt.I = sq$Binv
        #sqrt.I = expm::sqrtm(I)
        #inv.sqrt.I = solve(sqrt.I,tol=1e-40)
        e = inv.sqrt.I%*%S
        z = sqrt.I%*%p0up + e
        inv.pen.I = MASS::ginv(Ip,tol=1e-40)
        
        A    <- sqrt.I%*%inv.pen.I%*%sqrt.I
        ubre <- t(z - A%*%z)%*%(z - A%*%z) + 2*(sum(diag(A))) - npar
        return(ubre[[1]])
      }
      
      
      
      start_time <- Sys.time()
      
      
      while(count < Max){ #from 0 to 100
        
        no_cores <- detectCores() - 2 
        cl <- makeCluster(no_cores)
        clusterEvalQ(cl, library(splines2))
        clusterExport(cl,c("smooth1","smooth2","hsmooth1","hsmooth2","npar1","npar","pen1","pen2","mstatdat","nk","degree"), envir=environment())
        setDefaultCluster(cl=cl)
        
        
        results = optimParallel(
          par = p0,
          negloglik.base,
          gam=gam0,
          control = list(maxit = 10000, trace=1))
        
        
        results$gam=gam0
        
        count <- count + 1
        
        results$countruns = count
        
        print(count)
        
        results$diff = c(abs(exp(p0[1:npar]) - exp(results$par[1:npar])))
        
        max(results$diff)
        
        if((max(results$diff) < abstol) | (count==Max)) {
          results = optimParallel(
            par = p0,
            negloglik.base,
            gam=gam,
            control = list(maxit = 10000, trace=1),hessian=TRUE)
          stopCluster(cl)        
          
          results$gam=gam0
          
          break}
        
        stopCluster(cl)        
        
        p0 = results$par
        
        jj = jacob(c(exp(p0[1:npar])))
        
        S = colSums(jj)
        I = (t(jj)%*%jj)/nrow(mstatdat)
        
        penalty = optim(
          par=gam0,
          UBRE.g,
          p0=p0,
          method = "L-BFGS-B",
          lower=0,
          upper=20,
          control = list(maxit = 1000, trace=1))
        
        gam0 = penalty$par
        
        
      }
      
      
      stop_time <- Sys.time()
      
      results$time = stop_time-start_time
      
      results$covariance = MASS::ginv(results$hessian,tol=1e-40)
      
      if(is.null(maxt)){
        t <- seq(0, round(quantile(mstatdat$R1time,0.95)[[1]],0), by = 1)
      } else {
        t <- seq(0, maxt, by = 1)
      }
      
      newd = data.frame(t)
      par = exp(results$par)
      hazM1 = predict(hsmooth1, t)[,1:(nk-1)]
      cumhazM1 = predict(smooth1, t)[,1:(nk-1)]
      hazM2 = predict(hsmooth2, t)[,1:(nk-1)]
      cumhazM2 = predict(smooth2, t)[,1:(nk-1)]
      
      newd$haz12 = as.vector(hazM1 %*% par[1:npar1])
      newd$haz13 = as.vector(hazM2 %*% par[(npar1 + 1):npar])
      newd$cumhaz12 = as.vector(cumhazM1 %*% par[1:npar1])
      newd$cumhaz13 = as.vector(cumhazM2 %*% par[(npar1 + 1):npar])
      newd$surv = exp(-newd$cumhaz12 - newd$cumhaz13)
      
      for(i in 1:nrow(newd)) {
        
        x=newd[i,]
        
        cif12 <- function(z) {
          cif12 =  (exp(-predict(smooth1, z)[,1:(nk-1)] %*% par[1:npar1] - predict(smooth2, z)[,1:(nk-1)] %*% par[(npar1 + 1):npar])) * (predict(hsmooth1, z)[,1:(nk-1)] %*% par[1:npar1])
          return(cif12)
        }
        newd$cif12[i] = integrate(cif12, 0, x[["t"]])$value
      }
      
      
      for(i in 1:nrow(newd)) {
        
        x=newd[i,]
        
        cif13 <- function(z) {
          cif13 =  (exp(-predict(smooth1, z)[,1:(nk-1)] %*% par[1:npar1] - predict(smooth2, z)[,1:(nk-1)] %*%par[(npar1 + 1):npar])) * (predict(hsmooth2, z)[,1:(nk-1)] %*% par[(npar1 + 1):npar])
          
          return(cif13)}
        
        newd$cif13[i] = integrate(cif13, 0, x[["t"]])$value
      }
      
      
      newd$sum = newd$surv + newd$cif12 + newd$cif13
      
      ####PLOTS#####
      
      library(ggplot2)
      library(ggprism)
      library(ggpubr)
      
      
      trans.plot <- ggplot(data = newd) +
        #geom_line(aes(t,0,col="1->1"),size=1.2)+
        geom_line(aes(t, haz12, col = "Gametocyte initiation"), size = 1.2) +
        geom_line(aes(t, haz13, col = "Malaria resolution without gametocytes"), size = 1.2) +
        ylab("Hazard rate") +
        scale_color_manual(values = c("#E7B800", "#FC4E07")) +
        theme_prism(base_size = 16) +
        theme(legend.position = "bottom") +
        xlab("Time from incident malaria")
      
      #trans.plot+coord_cartesian(ylim=c(0,0.1))
      
      results$trans.plot <-trans.plot
      
      cumprob.plot <- ggplot(data = newd) +
        geom_line(aes(t, cif12, col = "Gametocyte initiation"), size = 1.2) +
        geom_line(aes(t, cif13, col = "Malaria resolution without gametocytes"), size = 1.2) +
        theme_bw() +
        ylab("Cumulative incidence") +
        scale_color_manual(values = c("#E7B800", "#FC4E07")) +
        theme_prism(base_size = 16) +
        theme(legend.position = "bottom") +
        xlab("Time from incident malaria")+
        scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1),labels = scales::percent_format(accuracy = 1))
      
      
      results$cumprob.plot <- cumprob.plot
      
      combiplot <-
        ggarrange(
          cumprob.plot+theme(axis.title.x = element_blank()),
          trans.plot, 
          nrow = 2, 
          labels = c("A","B"),common.legend = TRUE, align="v"
        ) 
      
      combiplot
      
      results$combiplot <- combiplot
      
      
      results$fitdata <- newd
      
      
    } else {
      
      # With covariates ---------------------------------------------------------
      
      
      if (length(startval1)!=npar1 | length(startval2)!=npar2){
        p0 = c(rep(startval1[1],npar1),rep(startval2[1],npar2),rep(0.1,ncov)) #start values
      } else {p0=c(startval1,startval2,rep(0.1,ncov))}
      
      
      
      gam0 = gam
      
      mstatdat_complete<-na.omit(mstatdat %>% dplyr::select("L0time","R0time","L1time","R1time","state","overlap",all_of(covs)))
      
      negloglik.covs<-function(par,gam){
        
        Pen = matrix(0,nrow=nparcov,ncol=nparcov)
        Pen[1:npar1,1:npar1]=pen1*exp(gam[1])
        Pen[(npar1+1):npar,(npar1+1):npar]=pen2*exp(gam[2])
        
        mstatdat_complete[["lp1"]]=as.vector(exp(as.matrix(mstatdat_complete[covs])%*%par[(npar+1):(npar+ncov1)]))
        mstatdat_complete[["lp2"]]=as.vector(exp(as.matrix(mstatdat_complete[covs])%*%par[(npar+ncov1+1):nparcov]))
        
        
        mstatdat_complete1 = mstatdat_complete[mstatdat_complete$state==1,]
        mstatdat_complete3 =  mstatdat_complete[mstatdat_complete$state==3,]
        mstatdat_complete2 =  mstatdat_complete[mstatdat_complete$state==2 & mstatdat_complete$overlap==0,] 
        mstatdat_complete2o = mstatdat_complete[mstatdat_complete$state==2 & mstatdat_complete$overlap==1,] 
        
        if (nrow(mstatdat_complete1)==0){L1=1} else{ 
          L1= ((1/6)*(exp((-predict(smooth1,mstatdat_complete1[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete1[["lp1"]]-predict(smooth2,mstatdat_complete1[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete1[["lp2"]])+
                        
                        4*exp((-predict(smooth1,mstatdat_complete1[["L1time"]]-mstatdat_complete1[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete1[["lp1"]]-predict(smooth2,mstatdat_complete1[["L1time"]]-mstatdat_complete1[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete1[["lp2"]])+exp((-predict(smooth1,mstatdat_complete1[["L1time"]]-mstatdat_complete1[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete1[["lp1"]]-predict(smooth2,mstatdat_complete1[["L1time"]]-mstatdat_complete1[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete1[["lp2"]])))}
        
        if (nrow(mstatdat_complete3)==0){L3=1} else{        
          L3= ((mstatdat_complete3[["R1time"]]-mstatdat_complete3[["L1time"]])/36)*(
            exp((-predict(smooth1,mstatdat_complete3[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,mstatdat_complete3[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,mstatdat_complete3[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              4*exp((-predict(smooth1,(mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,(mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,(mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              exp((-predict(smooth1,mstatdat_complete3[["R1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,mstatdat_complete3[["R1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,mstatdat_complete3[["R1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+ 
              4*exp((-predict(smooth1,mstatdat_complete3[["L1time"]]-mstatdat_complete3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,mstatdat_complete3[["L1time"]]-mstatdat_complete3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,mstatdat_complete3[["L1time"]]-mstatdat_complete3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              16*exp((-predict(smooth1,(mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,(mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,(mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              4*exp((-predict(smooth1,mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              exp((-predict(smooth1,mstatdat_complete3[["L1time"]]-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,mstatdat_complete3[["L1time"]]-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,mstatdat_complete3[["L1time"]]-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              4*exp((-predict(smooth1,((mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]])/2)-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,((mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]])/2)-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,((mstatdat_complete3[["L1time"]]+mstatdat_complete3[["R1time"]])/2)-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]]+
              exp((-predict(smooth1,mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete3[["lp1"]]-predict(smooth2,mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])*predict(hsmooth2,mstatdat_complete3[["R1time"]]-mstatdat_complete3[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete3[["lp2"]])}
        
        if (nrow(mstatdat_complete2)==0){L2=1} else{         
          L2 = ((mstatdat_complete2[["R1time"]]-mstatdat_complete2[["L1time"]])/36)*(
            exp((-predict(smooth1,mstatdat_complete2[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,mstatdat_complete2[["L1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,mstatdat_complete2[["L1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              4*exp((-predict(smooth1,(mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,(mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,(mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              exp((-predict(smooth1,mstatdat_complete2[["R1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,mstatdat_complete2[["R1time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,mstatdat_complete2[["R1time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+ 
              4*exp((-predict(smooth1,mstatdat_complete2[["L1time"]]-mstatdat_complete2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,mstatdat_complete2[["L1time"]]-mstatdat_complete2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,mstatdat_complete2[["L1time"]]-mstatdat_complete2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              16*exp((-predict(smooth1,(mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,(mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,(mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]])/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              4*exp((-predict(smooth1,mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              exp((-predict(smooth1,mstatdat_complete2[["L1time"]]-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,mstatdat_complete2[["L1time"]]-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,mstatdat_complete2[["L1time"]]-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              4*exp((-predict(smooth1,((mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]])/2)-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,((mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]])/2)-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,((mstatdat_complete2[["L1time"]]+mstatdat_complete2[["R1time"]])/2)-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]]+
              exp((-predict(smooth1,mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2[["lp1"]]-predict(smooth2,mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2[["lp2"]])*predict(hsmooth1,mstatdat_complete2[["R1time"]]-mstatdat_complete2[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2[["lp1"]])}
        
        if (nrow(mstatdat_complete2o)==0){L2o=1} else{        
          L2o = (mstatdat_complete2o[["R0time"]]/36)*(
            exp((-predict(smooth1,mstatdat_complete2o[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2o[["lp1"]]-predict(smooth2,mstatdat_complete2o[["R0time"]])[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2o[["lp2"]])*predict(hsmooth1,mstatdat_complete2o[["R0time"]])[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2o[["lp1"]]+   
              6*exp((-predict(smooth1,mstatdat_complete2o[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2o[["lp1"]]-predict(smooth2,mstatdat_complete2o[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2o[["lp2"]])*predict(hsmooth1,mstatdat_complete2o[["R0time"]]/2)[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2o[["lp1"]]+
              8*exp((-predict(smooth1,mstatdat_complete2o[["R0time"]]/4)[,1:(nk-1)]%*%((exp(par[1:npar1]))))*mstatdat_complete2o[["lp1"]]-predict(smooth2,mstatdat_complete2o[["R0time"]]/4)[,1:(nk-1)]%*%((exp(par[(npar1+1):npar])))*mstatdat_complete2o[["lp2"]])*predict(hsmooth1,mstatdat_complete2o[["R0time"]]/4)[,1:(nk-1)]%*%((exp(par[1:npar1])))*mstatdat_complete2o[["lp1"]]
          )}
        
        loglik=sum(log(L1))+sum(log(L3))+sum(log(L2))+sum(log(L2o))-0.5*(t(exp(par[(1):nparcov]))%*%Pen%*%(exp(par[(1):nparcov])))
        
        return(-loglik)
      }
      
      fcovs<-function(x){
        
        Pen = matrix(0,nrow=nparcov,ncol=nparcov)
        Pen[1:npar1,1:npar1]=pen1*exp(gam[1])
        Pen[(npar1+1):npar,(npar1+1):npar]=pen2*exp(gam[2])
        
        mstatdat_complete[["lp1"]]=as.vector(exp(as.matrix(mstatdat_complete[covs])%*%x[(npar+1):(npar+ncov1)]))
        mstatdat_complete[["lp2"]]=as.vector(exp(as.matrix(mstatdat_complete[covs])%*%x[(npar+ncov1+1):nparcov]))
        
        
        
        
        L1= ((1/6)*(exp((-predict(smooth1,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])+4*exp((-predict(smooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])+exp((-predict(smooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])))
        
        
        L3= ((mstatdat_complete[["R1time"]]-mstatdat_complete[["L1time"]])/36)*(
          exp((-predict(smooth1,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            4*exp((-predict(smooth1,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            exp((-predict(smooth1,mstatdat_complete[["R1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,mstatdat_complete[["R1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+ 
            4*exp((-predict(smooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            16*exp((-predict(smooth1,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            4*exp((-predict(smooth1,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            exp((-predict(smooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            4*exp((-predict(smooth1,((mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,((mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,((mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]]+
            exp((-predict(smooth1,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth2,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])
        
        
        L2 = ((mstatdat_complete[["R1time"]]-mstatdat_complete[["L1time"]])/36)*(
          exp((-predict(smooth1,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["L1time"]])[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            4*exp((-predict(smooth1,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            exp((-predict(smooth1,mstatdat_complete[["R1time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R1time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["R1time"]])[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+ 
            4*exp((-predict(smooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            16*exp((-predict(smooth1,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,(mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])/2)[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            4*exp((-predict(smooth1,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            exp((-predict(smooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["L1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            4*exp((-predict(smooth1,((mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,((mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,((mstatdat_complete[["L1time"]]+mstatdat_complete[["R1time"]])/2)-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            exp((-predict(smooth1,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["R1time"]]-mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]])
        
        L2o = (mstatdat_complete[["R0time"]]/36)*(
          exp((-predict(smooth1,mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["R0time"]])[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+   
            6*exp((-predict(smooth1,mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["R0time"]]/2)[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]+
            8*exp((-predict(smooth1,mstatdat_complete[["R0time"]]/4)[,1:(nk-1)]%*%(((x[1:npar1]))))*mstatdat_complete[["lp1"]]-predict(smooth2,mstatdat_complete[["R0time"]]/4)[,1:(nk-1)]%*%(((x[(npar1+1):npar])))*mstatdat_complete[["lp2"]])*predict(hsmooth1,mstatdat_complete[["R0time"]]/4)[,1:(nk-1)]%*%(((x[1:npar1])))*mstatdat_complete[["lp1"]]
        )
        
        loglik=sum(log((L1^(I(mstatdat_complete[["state"]]==1)) )*  
                         (L3^(I(mstatdat_complete[["state"]]==3)))*
                         (L2^(I(mstatdat_complete[["state"]]==2)*I(mstatdat_complete[["overlap"]]==0)))*
                         (L2o^(I(mstatdat_complete[["state"]]==2)*I(mstatdat_complete[["overlap"]]==1)))))
        
        return(loglik[[1]])
      }
      
      
      library(pracma)
      jacob = function(par){
        pracma::jacobian(fcovs,x0=par,h=.Machine$double.eps^(1/3))
      }
      
      
      # Prelim algorithm:
      count    <- 0
      Max      <- maxit
      abstol   <- abstol
      
      jj = jacob(c(exp(p0[1:npar]),p0[(npar+1):nparcov]))
      S = colSums(jj)
      I = (t(jj)%*%jj)/nrow(mstatdat_complete)
      
      sqrtm = function (A, kmax = 20, tol = .Machine$double.eps^(1/2)) 
      {
        stopifnot(is.numeric(A), is.matrix(A))
        if (nrow(A) != ncol(A)) 
          stop("Matrix 'A' must be square.")
        P0 <- A
        Q0 <- diag(nrow(A))
        k <- 0
        while (norm(A - P0 %*% P0, "F") > tol && k < kmax) {
          P1 <- 0.5 * (P0 + MASS::ginv(Q0, tol = tol))
          Q1 <- 0.5 * (Q0 + MASS::ginv(P0, tol = tol))
          P0 <- P1
          Q0 <- Q1
          k <- k + 1
        }
        if (k >= kmax)
          k <- -1
        return(list(B = P0, Binv = Q0, k = k, acc = norm(A - P0 %*% 
                                                           P0, "F")))
      }
      
      
      UBRE.g <- function(gam,p0){
        Pen = matrix(0,nrow=nparcov,ncol=nparcov)
        Pen[1:npar1,1:npar1]=pen1*exp(gam[1])
        Pen[(npar1+1):npar,(npar1+1):npar]=pen2*exp(gam[2])
        
        p0up = c(exp(p0[1:npar]),p0[(npar+1):nparcov])
        
        Sp = S-Pen%*%p0up
        Ip = I+Pen
        sq = sqrtm(I,tol=1e-40)
        sqrt.I = sq$B
        inv.sqrt.I = sq$Binv
        #sqrt.I = expm::sqrtm(I)
        #inv.sqrt.I = solve(sqrt.I,tol=1e-40)
        e = inv.sqrt.I%*%S
        z = sqrt.I%*%p0up + e
        inv.pen.I = MASS::ginv(Ip,tol=1e-40)
        
        A    <- sqrt.I%*%inv.pen.I%*%sqrt.I
        ubre <- t(z - A%*%z)%*%(z - A%*%z) + 2*(sum(diag(A))) - nparcov
        return(ubre[[1]])
      }
      
      
      
      start_time <- Sys.time()
      
      
      while(count < Max){ #from 0 to 100
        
        no_cores <- detectCores() - 2 
        cl <- makeCluster(no_cores)
        clusterEvalQ(cl, library(splines2))
        clusterExport(cl,c("smooth1","smooth2","hsmooth1","hsmooth2","npar1","npar","covs","ncov1","nparcov","pen1","pen2","mstatdat_complete","nk","degree"), envir=environment())
        setDefaultCluster(cl=cl)
        
        
        results = optimParallel(
          par = p0,
          negloglik.covs,
          gam=gam0,
          control = list(maxit = 10000, trace=1))
        
        
        results$gam=gam0
        
        count <- count + 1
        
        results$countruns = count
        
        print(count)
        
        results$diff = c(abs(exp(p0[1:npar]) - exp(results$par[1:npar])),abs(p0[(npar+1):nparcov]-results$par[(npar+1):nparcov]))
        
        max(results$diff)
        
        if((max(results$diff) < abstol) | (count==Max)) {
          results = optimParallel(
            par = p0,
            negloglik.covs,
            gam=gam,
            control = list(maxit = 10000, trace=1),hessian=TRUE)
          stopCluster(cl)        
          
          results$gam=gam0
          
          break}
        
        stopCluster(cl)        
        
        p0 = results$par
        
        jj = jacob(c(exp(p0[1:npar]),p0[(npar+1):nparcov]))
        
        S = colSums(jj)
        I = (t(jj)%*%jj)/nrow(mstatdat_complete)
        
        penalty = optim(
          par=gam0,
          UBRE.g,
          p0=p0,
          method = "L-BFGS-B",
          lower=0,
          upper=20,
          control = list(maxit = 1000, trace=1))
        
        gam0 = penalty$par
        
        
      }
      
      
      stop_time <- Sys.time()
      
      results$time = stop_time-start_time
      
      
      results$covariance<- MASS::ginv(results$hessian,tol=1e-40)
      prop_sigma<-sqrt(diag(results$covariance)[(npar+1):(nparcov)])
      
      estimates=data.frame(trans=rep(c("1->2","1->3"),each=ncov1),covs=rep(covs,2),
                           coef=results$par[(npar+1):nparcov],SE=prop_sigma) %>%
        mutate(coef.lci=coef-1.96*SE,coef.uci=coef+1.96*SE,HR=exp(coef),HR.lci=exp(coef.lci),HR.uci=exp(coef.uci),p.val=2*pnorm(-abs(coef/SE)))
      
      results$estimates=estimates
      
      
      
      results$AIC = 2*(length(results$par)+results$value)
      
      if(pred==TRUE){     
        par = c(exp(results$par[1:npar]),as.vector(na.exclude(results$par[(npar+1):(npar+nparcov)])))
        
        
        if(is.null(maxt)){
          t <- seq(0, round(quantile(mstatdat_complete$R1time,0.95)[[1]],0), by = 1)
        } else {
          t <- seq(0, maxt, by = 1)
        }
        
        expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
        newd= expand.grid.df(data.frame(t),
                             unique(mstatdat_complete[covs]))
        
        newd[["lp1"]]=as.vector(exp(as.matrix(newd[covs])%*%par[(npar+1):(npar+ncov1)]))
        newd[["lp2"]]=as.vector(exp(as.matrix(newd[covs])%*%par[(npar+ncov1+1):nparcov]))
        
        hazM1 = predict(hsmooth1, newd$t)[,1:(nk-1)]
        cumhazM1 = predict(smooth1, newd$t)[,1:(nk-1)]
        hazM2 = predict(hsmooth2, newd$t)[,1:(nk-1)]
        cumhazM2 = predict(smooth2, newd$t)[,1:(nk-1)]
        
        newd$haz12 = as.vector(hazM1 %*% par[1:npar1])*newd[["lp1"]]
        newd$haz13 = as.vector(hazM2 %*% par[(npar1 + 1):npar])*newd[["lp2"]]
        newd$cumhaz12 = as.vector(cumhazM1 %*% par[1:npar1])*newd[["lp1"]]
        newd$cumhaz13 = as.vector(cumhazM2 %*% par[(npar1 + 1):npar])*newd[["lp2"]]
        newd$surv = exp(-newd$cumhaz12 - newd$cumhaz13)
        
        for(i in 1:nrow(newd)) {
          
          x=newd[i,]
          
          cif12 <- function(z) {
            cif12 =  (exp(-(predict(smooth1, z)[,1:(nk-1)] %*% par[1:npar1])*x[["lp1"]] - (predict(smooth2, z)[,1:(nk-1)] %*% par[(npar1 + 1):npar])*x[["lp2"]])) * (predict(hsmooth1, z)[,1:(nk-1)] %*% par[1:npar1])*x[["lp1"]]
            return(cif12)
          }
          newd$cif12[i] = integrate(cif12, 0, x[["t"]])$value
        }
        
        
        for(i in 1:nrow(newd)) {
          
          x=newd[i,]
          
          cif13 <- function(z) {
            cif13 =  (exp(-(predict(smooth1, z)[,1:(nk-1)] %*% par[1:npar1])*x[["lp1"]] - (predict(smooth2, z)[,1:(nk-1)] %*%par[(npar1 + 1):npar])*x[["lp2"]])) * (predict(hsmooth2, z)[,1:(nk-1)] %*% par[(npar1 + 1):npar])*x[["lp2"]]
            
            return(cif13)}
          
          newd$cif13[i] = integrate(cif13, 0, x[["t"]])$value
        }
        
        
        newd$sum = newd$surv + newd$cif12 + newd$cif13
        
        
        
        
        results$fitdata <- newd
        
      }
      
    }
    
    
  } else {
    if (method=="exp"){
      npar1<-1
      npar2<-1
      
      ######Likelihood function#####
      npar<-npar1+npar2
      
      
      negloglik.base<-function(par,data){
        
        epar=par[1:2]^2
        
        L1 = (epar[1]/(((epar[1]+epar[2])^2)*(data[["R0time"]]-data[["L0time"]])))*exp(-(epar[1]+epar[2])*(data[["R0time"]]-data[["L0time"]])) + (epar[1]/(epar[1]+epar[2])) - (epar[1]/(((epar[1]+epar[2])^2)*(data[["R0time"]]-data[["L0time"]])))
        
        L2 = (-epar[1]/(((epar[1]+epar[2])^2)*(data[["R0time"]]-data[["L0time"]])))*(exp(-(epar[1]+epar[2])*(data[["R1time"]]-data[["R0time"]])) - exp(-(epar[1]+epar[2])*(data[["R1time"]]-data[["L0time"]]))-exp(-(epar[1]+epar[2])*(data[["L1time"]]-data[["R0time"]]))+exp(-(epar[1]+epar[2])*(data[["L1time"]]-data[["L0time"]])))
        
        L3 = (-epar[2]/(((epar[1]+epar[2])^2)*(data[["R0time"]]-data[["L0time"]])))*(exp(-(epar[1]+epar[2])*(data[["R1time"]]-data[["R0time"]])) - exp(-(epar[1]+epar[2])*(data[["R1time"]]-data[["L0time"]]))-exp(-(epar[1]+epar[2])*(data[["L1time"]]-data[["R0time"]]))+exp(-(epar[1]+epar[2])*(data[["L1time"]]-data[["L0time"]])))
        
        L4 = (1/(((epar[1]+epar[2]))*(data[["R0time"]]-data[["L0time"]])))*(exp(-(epar[1]+epar[2])*(data[["L1time"]]-data[["R0time"]]))-exp(-(epar[1]+epar[2])*(data[["L1time"]]-data[["L0time"]])))
        
        lik = (L1^(data[["state"]]==2 & data[["overlap"]]==1))*(L2^(data[["state"]]==2 & data[["overlap"]]==0))*(L3^(data[["state"]]==3))*(L4^(data[["state"]]==1))
        
        lik = ifelse(lik==0,0.001,lik)
        
        loglik = sum(log(lik))
        
        return(-loglik)
      }
      
      
      
      # Now including covariates ------------------------------------------------
      
      ncov1<-length(covs)
      ncov2<-length(covs) #(covariates for each transition)
      ncov<-ncov1+ncov2
      nparcov<-npar+ncov
      
      
      negloglik.covs<-function(par,data){
        
        epar=par[1:2]^2
        
        data[["lp1"]]=as.vector(exp(as.matrix(data[covs])%*%par[(npar+1):(npar+ncov1)]))
        data[["lp2"]]=as.vector(exp(as.matrix(data[covs])%*%par[(npar+ncov1+1):nparcov]))
        
        L1 = ((epar[1]*data[["lp1"]])/((((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))^2)*(data[["R0time"]]-data[["L0time"]])))*exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["R0time"]]-data[["L0time"]])) + ((epar[1]*data[["lp1"]])/((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))) - ((epar[1]*data[["lp1"]])/((((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))^2)*(data[["R0time"]]-data[["L0time"]])))
        
        L2 = (-(epar[1]*data[["lp1"]])/((((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))^2)*(data[["R0time"]]-data[["L0time"]])))*(exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["R1time"]]-data[["R0time"]])) - exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["R1time"]]-data[["L0time"]]))-exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["L1time"]]-data[["R0time"]]))+exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["L1time"]]-data[["L0time"]])))
        
        L3 = (-(epar[2]*data[["lp2"]])/((((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))^2)*(data[["R0time"]]-data[["L0time"]])))*(exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["R1time"]]-data[["R0time"]])) - exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["R1time"]]-data[["L0time"]]))-exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["L1time"]]-data[["R0time"]]))+exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["L1time"]]-data[["L0time"]])))
        
        L4 = (1/((((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]])))*(data[["R0time"]]-data[["L0time"]])))*(exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["L1time"]]-data[["R0time"]]))-exp(-((epar[1]*data[["lp1"]])+(epar[2]*data[["lp2"]]))*(data[["L1time"]]-data[["L0time"]])))
        
        lik = (L1^(data[["state"]]==2 & data[["overlap"]]==1))*(L2^(data[["state"]]==2 & data[["overlap"]]==0))*(L3^(data[["state"]]==3))*(L4^(data[["state"]]==1))
        
        lik = ifelse(lik==0,0.001,lik)
        
        loglik = sum(log(lik))
        
        return(-loglik)
        
      }
      
      
      if (is.null(covs)){
        
        # No covariates results -----------------------------------------------------
        p0 = c(startval1,startval2) #start values
        
        library(optimParallel)  
        no_cores <- detectCores() - 2 
        cl <- makeCluster(no_cores)
        #clusterEvalQ(cl, library(splines2))
        clusterExport(cl,c("npar1","npar"), envir=environment())
        setDefaultCluster(cl=cl)
        
        
        time0<-Sys.time()
        results = optimParallel(
          par = p0,
          negloglik.base,
          data = mstatdat,
          #lower=rep(0,4),
          control = list(maxit = 10000, trace=1), hessian = hes)
        time1<-Sys.time()
        
        results$time = time1-time0
        
        stopCluster(cl)
        
        if(is.null(maxt)){
          t <- seq(0, round(quantile(mstatdat$R1time,0.95)[[1]],0), by = 1)
        } else {
          t <- seq(0, maxt, by = 1)
        }
        
        newd = data.frame(t)
        par = results$par
        epar=(par[1:2])^2
        haz1 = function(t){(epar[1])}
        haz2 = function(t){(epar[2])}
        cumhaz1 = function(t){(epar[1]*(t))}
        cumhaz2 = function(t){(epar[2]*(t))}
        
        newd$haz12 = haz1(t)
        newd$haz13 = haz2(t)
        newd$cumhaz12 = cumhaz1(t)
        newd$cumhaz13 = cumhaz2(t)
        newd$surv = exp(-newd$cumhaz12 - newd$cumhaz13)
        
        newd$cif12 = (newd$haz12/(newd$haz12+newd$haz13))-(newd$haz12/(newd$haz12+newd$haz13))*newd$surv
        newd$cif13 = (newd$haz13/(newd$haz12+newd$haz13))-(newd$haz13/(newd$haz12+newd$haz13))*newd$surv        
        
        newd$sum = newd$surv + newd$cif12 + newd$cif13
        
        ####PLOTS#####
        
        library(ggplot2)
        library(ggprism)
        library(ggpubr)
        
        #newd$haz12=ifelse(newd$t>60,NA,newd$haz12)
        
        results$fitdata <- newd
        
        trans.plot <- ggplot(data = newd) +
          #geom_line(aes(t,0,col="1->1"),size=1.2)+
          geom_line(aes(t, haz12, col = "Gametocyte initiation"), size = 1.2) +
          geom_line(aes(t, haz13, col = "Malaria resolution without gametocytes"), size = 1.2) +
          ylab("Hazard rate") +
          scale_color_manual(values = c("#E7B800", "#FC4E07")) +
          theme_prism(base_size = 16) +
          theme(legend.position = "bottom") +
          xlab("Time from incident malaria")
        
        results$trans.plot <-trans.plot
        
        cumprob.plot <- ggplot(data = newd) +
          geom_line(aes(t, cif12, col = "Gametocyte initiation"), size = 1.2) +
          geom_line(aes(t, cif13, col = "Malaria resolution without gametocytes"), size = 1.2) +
          theme_bw() +
          ylab("Cumulative incidence") +
          scale_color_manual(values = c("#E7B800", "#FC4E07")) +
          theme_prism(base_size = 16) +
          theme(legend.position = "bottom") +
          xlab("Time from incident malaria")+
          scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1),labels = scales::percent_format(accuracy = 1))
        
        
        results$cumprob.plot <- cumprob.plot
        
        combiplot <-
          ggarrange(
            cumprob.plot+theme(axis.title.x = element_blank()),
            trans.plot, 
            nrow = 2, 
            labels = c("A","B"),common.legend = TRUE, align="v"
          ) 
        
        combiplot
        
        results$combiplot <- combiplot
        
        
        
        
      } else {
        
        # With covariates ---------------------------------------------------------
        p0 = c(startval1,startval2,rep(0.1,ncov)) #start values
        #start values
        
        mstatdat_complete=na.omit(mstatdat %>% dplyr::select("L0time","R0time","L1time","R1time","state","overlap",all_of(covs)))
        
        no_cores <- detectCores() - 2 
        cl <- makeCluster(no_cores)
        #clusterEvalQ(cl, library(splines2))
        clusterExport(cl,c("npar1","npar"), envir=environment())
        setDefaultCluster(cl=cl)
        
        covs <- covs
        ncov1<-length(covs)
        ncov2<-length(covs) #(covariates for each transition)
        ncov<-ncov1+ncov2
        nparcov<-npar+ncov
        
        time0<-Sys.time()
        results = optimParallel(
          par = p0,
          negloglik.covs,
          data = mstatdat_complete,
          control = list(maxit = 10000, trace=1), hessian = hes)
        time1<-Sys.time()
        
        results$time = time1-time0
        
        stopCluster(cl)
        
        if (hes){
          fisher_info<-solve(results$hessian)
          prop_sigma<-sqrt(diag(fisher_info)[(npar+1):nparcov])
          
          estimates=data.frame(trans=rep(c("1->2","1->3"),each=ncov1),covs=rep(covs,2),
                               coef=results$par[(npar+1):nparcov],SE=prop_sigma) %>%
            mutate(coef.lci=coef-1.96*SE,coef.uci=coef+1.96*SE,HR=exp(coef),HR.lci=exp(coef.lci),HR.uci=exp(coef.uci),p.val=2*pnorm(-abs(coef/SE)))
          
          results$estimates=estimates
        }
        
        if(pred==TRUE){
          
          par = results$par
          
          if(is.null(maxt)){
            t <- seq(0, round(quantile(mstatdat_complete$R1time,0.95)[[1]],0), by = 1)
          } else {
            t <- seq(0, maxt, by = 1)
          }
          
          expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
          newd= expand.grid.df(data.frame(t),
                               unique(mstatdat_complete[covs]))
          
          newd[["lp1"]]=as.vector(exp(as.matrix(newd[covs])%*%par[(npar+1):(npar+ncov1)]))
          newd[["lp2"]]=as.vector(exp(as.matrix(newd[covs])%*%par[(npar+ncov1+1):nparcov]))
          
          epar=(par[1:2])^2
          haz1 = function(t){(epar[1])}
          haz2 = function(t){(epar[2])}
          cumhaz1 = function(t){(epar[1]*(t))}
          cumhaz2 = function(t){(epar[2]*(t))}
          
          newd$haz12 = haz1(newd$t)*newd$lp1
          newd$haz13 = haz2(newd$t)*newd$lp2
          newd$cumhaz12 = cumhaz1(newd$t)*newd$lp1
          newd$cumhaz13 = cumhaz2(newd$t)*newd$lp2
          newd$surv = exp(-newd$cumhaz12 - newd$cumhaz13)
          newd$cif12 = (newd$haz12/(newd$haz12+newd$haz13))-(newd$haz12/(newd$haz12+newd$haz13))*newd$surv
          newd$cif13 = (newd$haz13/(newd$haz12+newd$haz13))-(newd$haz13/(newd$haz12+newd$haz13))*newd$surv
          newd$sum = newd$surv + newd$cif12 + newd$cif13
          
          results$fitdata = newd
          
        }        
      }
      
      results$AIC = 2*(length(results$par)+results$value)
      
      #95% CI for hazard rate parameters
      #since we used squared parameters, we can use a log transform
      
      pars = function(x){log(x^2)}
      library(rootSolve)
      grad = rootSolve::gradient(pars,x=results$par[1:npar])
      SE   = sqrt(diag(grad%*%solve(results$hessian[1:npar,1:npar])%*%t(grad)))
      results$baseline = data.frame(par=exp(pars(results$par[1:npar])),
                                    par.lci = exp(pars(results$par[1:npar])-1.96*SE),
                                    par.lci = exp(pars(results$par[1:npar])+1.96*SE))
      
      
      
    } else {
      if (method=="wei"){
        npar1<-2
        npar2<-2
        
        ######Likelihood function#####
        npar<-npar1+npar2
        
        
        negloglik.base<-function(par,data){
          
          lik=list()
          
          epar=exp(par[1:(npar)])
          
          for(i in 1:nrow(data)) {
            
            z=data[i,]
            
            haz1 = function(x,y){(epar[1]*(epar[2]^epar[1])*(y-x)^(epar[1]-1))}
            haz2 = function(x,y){(epar[3]*(epar[4]^epar[3])*(y-x)^(epar[3]-1))}
            cumhaz1 = function(x,y){(epar[2]*(y-x))^epar[1]}
            cumhaz2 = function(x,y){(epar[4]*(y-x))^epar[3]}
            
            p11=function(x){
              (1/(z[["R0time"]]-z[["L0time"]]))*(exp(-cumhaz1(x,z[["L1time"]])-cumhaz2(x,z[["L1time"]])))
            }
            
            p12<-function(x,y){
              (1/(z[["R0time"]]-z[["L0time"]]))*(exp(-cumhaz1(x,y)-cumhaz2(x,y)))*haz1(x,y)}
            
            p13<-function(x,y){
              (1/(z[["R0time"]]-z[["L0time"]]))*(exp(-cumhaz1(x,y)-cumhaz2(x,y)))*haz2(x,y)}
            
            if (z[["state"]]==1){
              L=integrate(p11,z[["L0time"]],z[["R0time"]])$value
            } else {
              if (z[["state"]]==3){
                L=pracma::integral2(p13,z[["L0time"]],z[["R0time"]],z[["L1time"]],z[["R1time"]],reltol=1e-4)$Q
              } else {
                if (z[["state"]]==2 & z["overlap"]==0){
                  L=pracma::integral2(p12,z[["L0time"]],z[["R0time"]],z[["L1time"]],z[["R1time"]],reltol=1e-4)$Q
                } else
                  L=pracma::integral2(p12,z[["L0time"]],z[["R0time"]],function(x){x},z[["R1time"]],reltol=1e-4)$Q
              }
              
            }
            lik[[i]]=log(L)
          }
          
          loglik=Reduce("+",lik)
          
          return(-loglik)
        }
        
        
        
        # Now including covariates ------------------------------------------------
        
        ncov1<-length(covs)
        ncov2<-length(covs) #(covariates for each transition)
        ncov<-ncov1+ncov2
        nparcov<-npar+ncov
        
        
        negloglik.covs<-function(par,data){
          
          epar=exp(par[1:(npar)])
          #epar=par
          
          lik=list()
          
          data[["lp1"]]=as.vector(exp(as.matrix(data[covs])%*%par[(npar+1):(npar+ncov1)]))
          data[["lp2"]]=as.vector(exp(as.matrix(data[covs])%*%par[(npar+ncov1+1):nparcov]))
          
          for (i in 1:nrow(data)) {
            
            z=data[i,]
            
            haz1 = function(x,y){(epar[1]*(epar[2]^epar[1])*(y-x)^(epar[1]-1))*z[["lp1"]]}
            haz2 = function(x,y){(epar[3]*(epar[4]^epar[3])*(y-x)^(epar[3]-1))*z[["lp2"]]}
            cumhaz1 = function(x,y){z[["lp1"]]*(epar[2]*(y-x))^epar[1]}
            cumhaz2 = function(x,y){z[["lp2"]]*(epar[4]*(y-x))^epar[3]}
            
            
            p11=function(x){
              (1/(z[["R0time"]]-z[["L0time"]]))*(exp(-cumhaz1(x,z[["L1time"]])-cumhaz2(x,z[["L1time"]])))
            }
            
            p12<-function(x,y){
              (1/(z[["R0time"]]-z[["L0time"]]))*(exp(-cumhaz1(x,y)-cumhaz2(x,y)))*haz1(x,y)}
            
            p13<-function(x,y){
              (1/(z[["R0time"]]-z[["L0time"]]))*(exp(-cumhaz1(x,y)-cumhaz2(x,y)))*haz2(x,y)}
            
            
            
            if (z[["state"]]==1){
              L=integrate(p11,z[["L0time"]],z[["R0time"]])$value
            } else {
              if (z[["state"]]==3){
                L=pracma::integral2(p13,z[["L0time"]],z[["R0time"]],z[["L1time"]],z[["R1time"]],reltol=1e-4)$Q
              } else {
                if (z[["state"]]==2 & z["overlap"]==0){
                  L=pracma::integral2(p12,z[["L0time"]],z[["R0time"]],z[["L1time"]],z[["R1time"]],reltol=1e-4)$Q
                } else
                  L=pracma::integral2(p12,z[["L0time"]],z[["R0time"]],function(x){x},z[["R1time"]],reltol=1e-4)$Q
              }
              
            }
            
            lik[[i]]=log(L)
          }
          
          loglik=Reduce("+",lik)
          
          return(-loglik)
        }
        
        
        if (is.null(covs)){
          
          # No covariates results -----------------------------------------------------
          p0 = c(startval1,startval2) #start values
          
          library(optimParallel)  
          no_cores <- detectCores() - 2 
          cl <- makeCluster(no_cores)
          #clusterEvalQ(cl, library(splines2))
          clusterExport(cl,c("npar1","npar"), envir=environment())
          setDefaultCluster(cl=cl)
          
          
          time0<-Sys.time()
          results = optimParallel(
            par = p0,
            negloglik.base,
            data = mstatdat,
            #lower=rep(0,4),
            control = list(maxit = 10000, trace=1), hessian = hes)
          time1<-Sys.time()
          
          results$time = time1-time0
          
          stopCluster(cl)
          
          if(is.null(maxt)){
            t <- seq(0, round(quantile(mstatdat$R1time,0.95)[[1]],0), by = 1)
          } else {
            t <- seq(0, maxt, by = 1)
          }
          
          newd = data.frame(t)
          par = results$par
          epar=exp(par)
          haz1 = function(t){(epar[1]*(epar[2]^epar[1])*(t)^(epar[1]-1))}
          haz2 = function(t){(epar[3]*(epar[4]^epar[3])*(t)^(epar[3]-1))}
          cumhaz1 = function(t){(epar[2]*(t))^epar[1]}
          cumhaz2 = function(t){(epar[4]*(t))^epar[3]}
          
          newd$haz12 = haz1(t)
          newd$haz13 = haz2(t)
          newd$cumhaz12 = cumhaz1(t)
          newd$cumhaz13 = cumhaz2(t)
          newd$surv = exp(-newd$cumhaz12 - newd$cumhaz13)
          
          
          integratey = function(f,a,b){
            if(a==b){0} else{
              integrate(f,a,b)$value
            }
          }
          
          
          for(i in 1:nrow(newd)) {
            x=newd[i,]
            
            cif12 <- function(z) {
              cif12 =  (exp(-cumhaz1(z) - cumhaz2(z)))*haz1(z)
              return(cif12)
            }
            newd$cif12[i] = integratey(cif12, 0, x[["t"]])
          }
          
          for(i in 1:nrow(newd)) {
            
            x=newd[i,]
            
            cif13 <- function(z) {
              cif13 =  (exp(-cumhaz1(z) - cumhaz2(z)))*haz2(z)
              
              return(cif13)}
            
            newd$cif13[i] = integratey(cif13, 0, x[["t"]])
          }
          
          
          newd$sum = newd$surv + newd$cif12 + newd$cif13
          
          ####PLOTS#####
          
          library(ggplot2)
          library(ggprism)
          library(ggpubr)
          
          #newd$haz12=ifelse(newd$t>60,NA,newd$haz12)
          
          results$fitdata <- newd
          
          trans.plot <- ggplot(data = newd) +
            #geom_line(aes(t,0,col="1->1"),size=1.2)+
            geom_line(aes(t, haz12, col = "Gametocyte initiation"), size = 1.2) +
            geom_line(aes(t, haz13, col = "Malaria resolution without gametocytes"), size = 1.2) +
            ylab("Hazard rate") +
            scale_color_manual(values = c("#E7B800", "#FC4E07")) +
            theme_prism(base_size = 16) +
            theme(legend.position = "bottom") +
            xlab("Time from incident malaria")
          
          results$trans.plot <-trans.plot
          
          cumprob.plot <- ggplot(data = newd) +
            geom_line(aes(t, cif12, col = "Gametocyte initiation"), size = 1.2) +
            geom_line(aes(t, cif13, col = "Malaria resolution without gametocytes"), size = 1.2) +
            theme_bw() +
            ylab("Cumulative incidence") +
            scale_color_manual(values = c("#E7B800", "#FC4E07")) +
            theme_prism(base_size = 16) +
            theme(legend.position = "bottom") +
            xlab("Time from incident malaria")+
            scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1),labels = scales::percent_format(accuracy = 1))
          
          
          results$cumprob.plot <- cumprob.plot
          
          combiplot <-
            ggarrange(
              cumprob.plot+theme(axis.title.x = element_blank()),
              trans.plot, 
              nrow = 2, 
              labels = c("A","B"),common.legend = TRUE, align="v"
            ) 
          
          combiplot
          
          results$combiplot <- combiplot
          
          
          
          
        } else {
          
          # With covariates ---------------------------------------------------------
          p0 = c(startval1,startval2,rep(0.1,ncov)) #start values
          #start values
          
          mstatdat_complete=na.omit(mstatdat %>% dplyr::select("L0time","R0time","L1time","R1time","state","overlap",all_of(covs)))
          
          no_cores <- detectCores() - 2 
          cl <- makeCluster(no_cores)
          #clusterEvalQ(cl, library(splines2))
          clusterExport(cl,c("npar1","npar"), envir=environment())
          setDefaultCluster(cl=cl)
          
          ncov1<-length(covs)
          ncov2<-length(covs) #(covariates for each transition)
          ncov<-ncov1+ncov2
          nparcov<-npar+ncov
          
          time0<-Sys.time()
          results = optimParallel(
            par = p0,
            negloglik.covs,
            data = mstatdat_complete,
            control = list(maxit = 10000, trace=1), hessian = hes)
          time1<-Sys.time()
          
          results$time = time1-time0
          
          stopCluster(cl)
          
          if (hes){
            fisher_info<-solve(results$hessian)
            prop_sigma<-sqrt(diag(fisher_info)[(npar+1):nparcov])
            
            estimates=data.frame(trans=rep(c("1->2","1->3"),each=ncov1),covs=rep(covs,2),
                                 coef=results$par[(npar+1):nparcov],SE=prop_sigma) %>%
              mutate(coef.lci=coef-1.96*SE,coef.uci=coef+1.96*SE,HR=exp(coef),HR.lci=exp(coef.lci),HR.uci=exp(coef.uci),p.val=2*pnorm(-abs(coef/SE)))
            
            results$estimates=estimates
          }
          
          
          if(pred==TRUE){
            
            par = results$par
            
            if(is.null(maxt)){
              t <- seq(0, round(quantile(mstatdat_complete$R1time,0.95)[[1]],0), by = 1)
            } else {
              t <- seq(0, maxt, by = 1)
            }
            
            expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
            newd= expand.grid.df(data.frame(t),
                                 unique(mstatdat_complete[covs]))
            
            newd[["lp1"]]=as.vector(exp(as.matrix(newd[covs])%*%par[(npar+1):(npar+ncov1)]))
            newd[["lp2"]]=as.vector(exp(as.matrix(newd[covs])%*%par[(npar+ncov1+1):nparcov]))
            
            epar=exp(par)
            haz1 = function(t){(epar[1]*(epar[2]^epar[1])*(t)^(epar[1]-1))}
            haz2 = function(t){(epar[3]*(epar[4]^epar[3])*(t)^(epar[3]-1))}
            cumhaz1 = function(t){(epar[2]*(t))^epar[1]}
            cumhaz2 = function(t){(epar[4]*(t))^epar[3]}
            
            newd$haz12 = haz1(t)*newd[["lp1"]]
            newd$haz13 = haz2(t)*newd[["lp2"]]
            newd$cumhaz12 = cumhaz1(t)*newd[["lp1"]]
            newd$cumhaz13 = cumhaz2(t)*newd[["lp2"]]
            newd$surv = exp(-newd$cumhaz12 - newd$cumhaz13)
            
            
            integratey = function(f,a,b){
              if(a==b){0} else{
                integrate(f,a,b)$value
              }
            }            
            
            for(i in 1:nrow(newd)) {
              x=newd[i,]
              
              cif12 <- function(z) {
                cif12 =  (exp(-cumhaz1(z)*x[["lp1"]] - cumhaz2(z)*x[["lp2"]]))*haz1(z)*x[["lp1"]]
                return(cif12)
              }
              newd$cif12[i] = integratey(cif12, 0, x[["t"]])
            }
            
            for(i in 1:nrow(newd)) {
              
              x=newd[i,]
              
              cif13 <- function(z) {
                cif13 =  (exp(-cumhaz1(z)*x[["lp1"]] - cumhaz2(z)*x[["lp2"]]))*haz2(z)*x[["lp2"]]
                
                return(cif13)}
              
              newd$cif13[i] = integratey(cif13, 0, x[["t"]])
            }
            
            
            newd$sum = newd$surv + newd$cif12 + newd$cif13
            
            results$fitdata = newd
            
            
          }
          
          
          
        }
        
        
        results$AIC = 2*(length(results$par)+results$value)
        
        #95% CI for hazard rate parameters
        #since we used squared parameters, we can use a log transform
        
        pars = function(x){exp(x)}
        SE   = sqrt(diag(solve(results$hessian[1:npar,1:npar])))
        results$baseline = data.frame(par=pars(results$par[1:npar]),
                                      par.lci = pars(results$par[1:npar]-1.96*SE),
                                      par.lci = pars(results$par[1:npar]+1.96*SE))
        
        
        
        
        
      } 
    }
  }
  
  return(results)
  
}








