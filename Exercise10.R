rm(list=ls())
setwd("/Users/ta4ha/Documents/Biocomputing/Biocomp-Fall2018-181109-Exercise10")

# Problem 1

# Load packages
library(ggplot2)

# Load data
data <- read.table("data.txt",header = TRUE, sep = ",")

# Create likelihood functions
linearMod<-function(p,x,y){
  B0=p[1]
  B1=p[2]
  sigma=exp(p[3])
  
  pred=B0+B1*x
  nll=-sum(dnorm(x=y,mean=pred,sd=sigma,log=TRUE))
  
  return(nll)
}

quadraticMod<-function(p,x,y,z){
  B0=p[1]
  B1=p[2]
  B2=p[3]
  sigma=exp(p[4])
  
  pred=B0+B1*x+B2*x^2
  nll=-sum(dnorm(x=y,mean=pred,sd=sigma,log=TRUE))
  
  return(nll)
}

# estimate parameters
linearGuess=c(12,12,1)
quadraticGuess=c(12,12,12,1)

fitLinear=optim(par=linearGuess,fn=linearMod,x=data$x,y=data$y)
fitQuadratic=optim(par=quadraticGuess,fn=quadraticMod,x=data$x,y=data$y)

# Run likelihood ratio test

teststat=2*(fitLinear$value-fitQuadratic$value)
# df=length(fitLinear$par)-length(fitQuadratic$par)
df = 1
1-pchisq(teststat,df)

# We would expect the quadratic to be a better fit and return a lower
# nll, but such is not the case. The linear model returns a lower nll.
# If the quadratic had been a better fit, we would have seen p < .05.
# Instead the liklihood ratio test returns 1.


#Problem 2

library(deSolve)

# Custom function that defines model differential equations
LVsim<-function(t,y,p){
  N1=y[1]
  N2=y[2]
  
  R1=p[1]
  K1=p[2]
  a21=p[3]
  R2=p[4]
  K2=p[5]
  a12=p[6]
  
  dN1dt=R1*(1-(N1+a12*N2)/K1)*N1
  dN2dt=R2*(1-(N2+a21*N1)/K2)*N2
  
  return(list(c(dN1dt,dN2dt)))
}

# Simulation 1 - N2 species goes to extinction
times=1:100
y0=c(0.1,0.1)
params1=c(0.5,10,2,0.5,10,0.5)
LVsim1=ode(y=y0,times=times,func=LVsim,parms=params1)
LVout1=data.frame(time=LVsim1[,1],normal=LVsim1[,2],tumor=LVsim1[,3])
ggplot(LVout1,aes(x=time,y=normal))+geom_line()+geom_line(data=LVout1,mapping=aes(x=time,y=tumor),col='red')+theme_classic()

# Simulation 2 - both species stable
times=1:400
y0=c(0.05,0.3)
params2=c(0.5,10,0.5,0.5,10,0.5)
LVsim2=ode(y=y0,times=times,func=LVsim,parms=params2)
LVout2=data.frame(time=LVsim2[,1],normal=LVsim2[,2],tumor=LVsim2[,3])
ggplot(LVout2,aes(x=time,y=normal))+geom_line()+geom_line(data=LVout2,mapping=aes(x=time,y=tumor),col='red')+theme_classic()

# Simulation 3 - N1 species goes to extinction
times=1:100
y0=c(0.1,0.1)
params3=c(0.5,10,0.5,0.5,10,2)
LVsim3=ode(y=y0,times=times,func=LVsim,parms=params3)
LVout3=data.frame(time=LVsim3[,1],normal=LVsim3[,2],tumor=LVsim3[,3])
ggplot(LVout3,aes(x=time,y=normal))+geom_line()+geom_line(data=LVout3,mapping=aes(x=time,y=tumor),col='red')+theme_classic()
