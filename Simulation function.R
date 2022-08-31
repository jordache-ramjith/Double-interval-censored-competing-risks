
library(survival)
library(ggplot2)
library(ggprism)

simulation = function(nsim=100 ,n=100, maxt=200, beta1=0.2, beta2=-0.1,width=28){
  
  n = floor(n/2)*2 #make sure it's an even number
  
  library(dplyr)
  
  #use a 4-p log-logistic
  t <- seq(0, maxt, by = 1)
  a1 = 0.014
  p1 = 1.6
  l1=0.026
  h1 = (p1*(a1^p1)*t^(p1-1))/(1+(l1*t)^p1)
  ch1 = ((a1/l1)^p1)*log(1+(l1*t)^p1)
  
  
  #max point of peak (28 days = 4 weeks)
  ((p1-1)/(l1^p1))^(1/p1)
  
  a2 = 0.04
  p2 = 3
  l2 = 0.05
  h2 = ifelse(t<28,0,(p2*(a2^p2)*(t-28)^(p2-1))/(1+(l2*(t-28))^p2))
  ch2 = c(rep(0,29),((a2/l2)^p2)*log(1+(l2*(t[-(1:29)]-28))^p2))

  
  #max point
  ((p2-1)/(l2^p2))^(1/p2) + 28
  
  
  newd = expand.grid(t=t,x=c(0,1))
  newd$h1=h1
  newd$h2=h2
  newd$s=exp(-ch1-ch2)
  newd$haz1 = h1*exp(beta1*newd$x)
  newd$haz2 = h2*exp(beta2*newd$x)
  newd$cumhaz1 = ch1*exp(beta1*newd$x)
  newd$cumhaz2 = ch2*exp(beta2*newd$x)
  newd$surv = exp(-newd$cumhaz1 - newd$cumhaz2)
  
  newd = newd %>% group_by(x) %>% mutate(cif1=cumsum(surv*haz1),cif2=cumsum(surv*haz2))
  
  simplotdata <<- newd
  
  library(ggplot2)
  library(ggpubr)
  ggplot(data=newd, aes(x = t))+
    geom_line(aes(y = h1, col="Event 1"),size=1.2)+
    geom_line(aes(y=h2, col = "Event 2"),size=1.2)+
    theme_prism()+
    theme(legend.position = "top")+
    xlab("Time")+
    ylab("Hazard")+
    scale_color_viridis_d(begin=0.20, end=0.75)+
    geom_vline(aes(xintercept = ((p1-1)/(l1^p1))^(1/p1),col="Event 1"),linetype="dashed")+
    geom_vline(aes(xintercept = ((p2-1)/(l2^p2))^(1/p2) + 28,col="Event 2"),linetype="dashed")
    
    
    
  newd = newd %>% group_by(x) %>% mutate(surv_prob=lag(surv,default = NA)-surv)
  
  newd = newd[newd$t>0,]
  newd0 = newd[newd$x==0,]
  newd1 = newd[newd$x==1,]
  
  
  simdat = list()
  
  for (i in 1:nsim){
    
    simdat[[i]]=data.frame(id=1:n)
    
    simdat[[i]]$x = rep(c(0,1),each=n/2)
    
    
    simdat[[i]]$L0time = 0
    
    fixed = seq(0,728+5*width,width)
    visits=list()
    for(j in 1:n){
      set.seed(paste0(i,j,1))
      visits[[j]] = round(fixed + c(0,rnorm(length(fixed),0,3)[-1]))
      simdat[[i]]$R0time[j] = visits[[j]][2]
    }
    
    for (j in 1:n){
      set.seed(paste0(j,420,i))
      simdat[[i]]$W[j] = round(runif(1, 1, simdat[[i]]$R0time[j]))
    }
    
    #using method by Jan Beyersmann 2009 for the times and events simulation
    
    set.seed(paste0(i,1))
    simdat[[i]]$t[1:(n/2)] = sample(c(newd0$t,maxt+1), size = n/2,replace=TRUE, prob=c(newd0$surv_prob,1-sum(newd0$surv_prob)))
    set.seed(paste0(i,2))
    simdat[[i]]$t[(n/2+1):n] = sample(c(newd1$t,maxt+1), size = n/2,replace=TRUE, prob=c(newd1$surv_prob,1-sum(newd1$surv_prob)))
    
    simdat[[i]]$h1 = ((p1*a1^p1*(simdat[[i]]$t^(p1-1)))/(1+(l1*simdat[[i]]$t)^p1))*exp(simdat[[i]]$x*beta1)
    simdat[[i]]$h2 = ifelse(simdat[[i]]$t<28,0,((p2*a2^p2*((simdat[[i]]$t-28)^(p2-1)))/(1+(l2*(simdat[[i]]$t-28))^p2))*exp(simdat[[i]]$x*beta2))
    
    simdat[[i]]$pr=simdat[[i]]$h1/(simdat[[i]]$h1+simdat[[i]]$h2)
    
    set.seed(paste0(i,3))
    simdat[[i]]$event = 2-rbinom(n,1,simdat[[i]]$pr)
    
    simdat[[i]] = simdat[[i]] %>% mutate(status= ifelse(t>maxt,0,1),
                                         t = ifelse(t > maxt, maxt, t),
                                         event=ifelse(status==0,0,event))
    
    simdat[[i]]$U=simdat[[i]]$W+simdat[[i]]$t
    
    for(j in 1:n){
      val3 = last (which(visits[[j]]<simdat[[i]]$U[j]))
      val4 = first(which(visits[[j]]>=simdat[[i]]$U[j]))
      simdat[[i]]$L1time[j]=visits[[j]][val3]
      simdat[[i]]$R1time[j]=visits[[j]][val4]
      
    }
    
    simdat[[i]] = simdat[[i]] %>% mutate(overlap=ifelse(R1time==R0time & L0time==L1time,1,0),
                                         gam_event=as.numeric(I(event==1)),
                                         res_event=as.numeric(I(event==2)),
                                         R1time=ifelse(status==0,L1time,R1time),
                                         state=ifelse(status==0,1,ifelse(gam_event==1,2,3)))
    
    print(i)
    
  }
  
  return(simdat)
}





