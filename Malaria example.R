source("model function.R")

faildt <- readRDS("faildt.RDS")

# Baseline ----------------------------------------------------------------
time0<-Sys.time()
out.exp=cr.interval(data=faildt,method="exp",
                    start1=c(0.1),start2=c(0.1),
                    maxt = 100)
time1<-Sys.time()
ittime=time1-time0

saveRDS(out.exp,"output/analysis_baseline_exp.RDS")

time0<-Sys.time()
out.wei=cr.interval(data=faildt,method="wei",
                    start1=c(0.1,-4),start2=c(0.1,-4),
                    maxt = 100)
time1<-Sys.time()
ittime=time1-time0

saveRDS(out.wei,"output/analysis_baseline_wei.RDS")

set.seed(pi)
start1=-abs(rnorm(10,mean=-5,1))

set.seed(pi*2)
start2=-abs(rnorm(10,mean=-5,1))

time0<-Sys.time()
out.bs=cr.interval(data=faildt,
                   method="bs",
                   start1=start1,
                   start2=start2,
                   abstol = 1e-6,
                   maxit=100,
                   maxt = 100)
out.bs$combiplot
time1<-Sys.time()
ittime=time1-time0

saveRDS(out.bs,"output/analysis_baseline_bs.RDS")


## need to combine all plots

ggarrange(
  out.exp$trans.plot+coord_cartesian(ylim=c(0,0.07)),
  out.wei$trans.plot+coord_cartesian(ylim=c(0,0.07)),
  out.bs$trans.plot+coord_cartesian(ylim=c(0,0.07)),
  out.exp$cumprob.plot,
  out.wei$cumprob.plot,
  out.bs$cumprob.plot,
  labels=c("A.","B.","C.","D.","E.","F."),
  nrow = 2, ncol=3,
  common.legend = TRUE,
  legend="top"
)


ggsave("results/Figure7.pdf",dpi=300,width=12,height=8)


## 15 basis functions to check knots placement.

### we allow randomly generated starting values. If the algorithm fails, a new set of starting values is used.

run=0
out.bs15=NA
maxr=20

while(any(is.na(out.bs15)) & run<=maxr){
  
  out.bs15=tryCatch(cr.interval(data=faildt,covs=NULL,method="bs",start1=-abs(rnorm(15,mean=-7,2)),start2=-abs(rnorm(15,mean=-7,2)),k=15,maxit=500,maxt = 100,abstol=1e-6,pred=TRUE), error=function(err) tryCatch(cr.interval(data=faildt,covs=NULL,method="bs",start1=-abs(rnorm(15,mean=-5,2)),start2=-abs(rnorm(15,mean=-5,2)),k=15,maxit=500,maxt = 100,abstol=1e-6,pred=TRUE), error=function(err) tryCatch(cr.interval(data=faildt,covs=NULL,method="bs",start1=-abs(rnorm(15,mean=-3,2)),start2=-abs(rnorm(15,mean=-3,2)),k=15,maxit=500,maxt = 100,abstol=1e-6,pred=TRUE), error=function(err) NA)))
  
  if(is.na(out.bs15)==T) {
    out.bs15=NA
  } else {
    if((any(is.nan(out.bs15$estimates$SE)) | any(out.bs15$estimates$SE<0.01))==T) {out.bs15=NA} else{print("YAY")}
  }
  print(run)
  run = run+1
}


out.bs$fitdata$t[out.bs$fitdata$haz12==max(out.bs$fitdata$haz12[out.bs$fitdata$t<70])]

out.bs15$fitdata$t[out.bs15$fitdata$haz12==max(out.bs15$fitdata$haz12[out.bs15$fitdata$t<70])]



saveRDS(out.bs15,"output/analysis_baseline_bs15.RDS")

out.bs15$combiplot
ggsave("results/supFigure15.pdf",dpi=300,width=12,height=8)


#AIC

2*length(out.exp$par)+2*out.exp$value
2*length(out.wei$par)+2*out.wei$value
2*length(out.bs$par)+2*out.bs$value

out.exp$value
out.wei$value
out.bs$value


# Technical ---------------------------------------------------------------

## to avoid issues of convergence, we write an algorithm to simulate start times from a normal distribution and repeat until we have no optimization failures.

splinemod<- function(covs=NULL){
  
  run=0
  out.bs2=NA
  maxr=20
  
  while(any(is.na(out.bs2)) & run<=maxr){
    
    out.bs2=tryCatch(cr.interval(data=faildt,covs=covs,method="bs",start1=-abs(rnorm(10,mean=-7,2)),start2=-abs(rnorm(10,mean=-7,2)),maxit=50,abstol=1e-5,pred=FALSE), error=function(err) tryCatch(cr.interval(data=faildt,covs=covs,method="bs",start1=-abs(rnorm(10,mean=-5,2)),start2=-abs(rnorm(10,mean=-5,2)),maxit=50,abstol=1e-5,pred=FALSE), error=function(err) tryCatch(cr.interval(data=faildt,covs=covs,method="bs",start1=-abs(rnorm(10,mean=-3,2)),start2=-abs(rnorm(10,mean=-3,2)),maxit=50,abstol=1e-5,pred=FALSE), error=function(err) NA)))
    
    if(is.na(out.bs2)==T) {
      out.bs2=NA
    } else {
      if((any(is.nan(out.bs2$estimates$SE)) | any(out.bs2$estimates$SE<0.01))==T) {out.bs2=NA} else{print("YAY")}
    }
    print(run)
    run = run+1
  }
  
  return(out.bs2)
  
}

# HBS --------------------------------------------------------------------

covs=c("hbs.homo","hbs.het")

out.exp.hbs = cr.interval(data=faildt,
                          method="exp",
                          covs=covs,
                          start1=c(0.1),
                          start2=c(0.1),
                          pred=FALSE)

out.exp.hbs$estimates


out.wei.hbs = cr.interval(data=faildt,
                          method="wei", 
                          covs=covs,
                          start1=c(0.1,-4),
                          start2=c(0.1,-4),
                          pred=FALSE)

out.wei.hbs$estimates

out.bs.hbs = splinemod(covs=covs)

out.bs.hbs$estimates


out.hbs = list(out.exp.hbs=out.exp.hbs,out.wei.hbs=out.wei.hbs,out.bs.hbs=out.bs.hbs)
saveRDS(out.hbs,"output/analysis_covs_hbs.RDS")


# male -------------------------------------------------------------------

covs=c("male")

out.exp.male = cr.interval(data=faildt,
                           method="exp",
                           covs=covs,
                           start1=c(0.1),
                           start2=c(0.1))
out.exp.male$estimates

out.wei.male = cr.interval(data=faildt,
                           method="wei", 
                           covs=covs,
                           start1=c(0.1,-4),
                           start2=c(0.1,-4))
out.wei.male$estimates



out.bs.male = splinemod(covs=covs)
out.bs.male$estimates

out.male = list(out.exp.male=out.exp.male,out.wei.male=out.wei.male,out.bs.male=out.bs.male)
saveRDS(out.male,"output/analysis_covs_male.RDS")



# qPCR -------------------------------------------------------------------

covs=c("logqpcr")

out.exp.qpcr = cr.interval(data=faildt,
                           method="exp",
                           covs=covs,
                           start1=c(0.1),
                           start2=c(0.1))
out.exp.qpcr$estimates

out.wei.qpcr = cr.interval(data=faildt,
                           method="wei", 
                           covs=covs,
                           start1=c(0.1,-4),
                           start2=c(0.1,-4))
out.wei.qpcr$estimates



out.bs.qpcr = splinemod(covs=covs)
out.bs.qpcr$estimates

out.qpcr = list(out.exp.qpcr=out.exp.qpcr,out.wei.qpcr=out.wei.qpcr,out.bs.qpcr=out.bs.qpcr)
saveRDS(out.qpcr,"output/analysis_covs_qpcr.RDS")


# age --------------------------------------------------------------------

covs=c("age2","age3")

out.exp.age = cr.interval(data=faildt,
                          method="exp",
                          covs=covs,
                          start1=c(0.1),
                          start2=c(0.1))
out.exp.age$estimates

out.wei.age = cr.interval(data=faildt,
                          method="wei", 
                          covs=covs,
                          start1=c(0.1,-4),
                          start2=c(0.1,-4))
out.wei.age$estimates



out.bs.age = splinemod(covs=covs)
out.bs.age$estimates

out.age = list(out.exp.age=out.exp.age,out.wei.age=out.wei.age,out.bs.age=out.bs.age)
saveRDS(out.age,"output/analysis_covs_age.RDS")




# moi --------------------------------------------------------------------

covs=c("moib")

out.exp.moi = cr.interval(data=faildt,
                          method="exp",
                          covs=covs,
                          start1=c(0.1),
                          start2=c(0.1))
out.exp.moi$estimates

out.wei.moi = cr.interval(data=faildt,
                          method="wei", 
                          covs=covs,
                          start1=c(0.1,-4),
                          start2=c(0.1,-4))
out.wei.moi$estimates


out.bs.moi = splinemod(covs=covs)
out.bs.moi$estimates

out.moi = list(out.exp.moi=out.exp.moi,out.wei.moi=out.wei.moi,out.bs.moi=out.bs.moi)
saveRDS(out.moi,"output/analysis_covs_moi.RDS")



# combine results --------------------------------------------------------

out.exp.hbs$estimates$method="exp"
out.wei.hbs$estimates$method="wei"
out.bs.hbs$estimates$method="spline"

out.exp.age$estimates$method="exp"
out.wei.age$estimates$method="wei"
out.bs.age$estimates$method="spline"

out.exp.male$estimates$method="exp"
out.wei.male$estimates$method="wei"
out.bs.male$estimates$method="spline"

out.exp.qpcr$estimates$method="exp"
out.wei.qpcr$estimates$method="wei"
out.bs.qpcr$estimates$method="spline"

out.exp.moi$estimates$method="exp"
out.wei.moi$estimates$method="wei"
out.bs.moi$estimates$method="spline"

out.all = rbind(
  out.exp.age$estimates,out.wei.age$estimates,out.bs.age$estimates,
  out.exp.male$estimates,out.wei.male$estimates,out.bs.male$estimates,
  out.exp.hbs$estimates,out.wei.hbs$estimates,out.bs.hbs$estimates,
  out.exp.qpcr$estimates,out.wei.qpcr$estimates,out.bs.qpcr$estimates,
  out.exp.moi$estimates,out.wei.moi$estimates,out.bs.moi$estimates
)%>%dplyr::select(trans,method,covs,HR,HR.lci,HR.uci,p.val)
out.all = out.all[order(out.all$trans),]
out.all

saveRDS(out.all,"output/analysis_all.RDS")


write.table(out.all, file = "results/allout.csv", append = FALSE, sep = ",", dec = ".",row.names = FALSE, col.names = TRUE)

