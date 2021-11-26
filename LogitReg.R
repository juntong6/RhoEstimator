####################data simulation################################
install.packages("cmaes",repos=structure(c(CRAN="https://packages.othr.de/cran/")))
install.packages("e1071",repos=structure(c(CRAN="https://packages.othr.de/cran/")))

library(stats4)
library(cmaes)
library(e1071)
k = 5 
n = 500 
vraieloi = t(c(1,1,1,1,1,1))
#thetaconta = t(c(-1,-1,-1,1,-1,1))
SimulationDesCaracteristiques = function(n,k) 
{
  y = matrix(0,nrow=n,ncol=k,byrow=T); 
  a = 0.25
  b = 2
  for (i in 1:n)
  {
    p1 = t(a*(2*runif(k,0,1)-1))
    p2 = t(0.25*(2*runif(k,0,1)-1)+b)
    p3 = t(0.25*(2*runif(k,0,1)-1)-b)
    p = runif(1,0,1)
    if(p<=(1/3))
      y[i,] = p1
    else if((p<=(2/3))&&(p>=(1/3)))
      y[i,] = p2  
    else
      y[i,] = p3 
    
  }
  constcol= matrix(1,nrow=n,ncol=1)
  y=cbind(constcol,y)
  return(y)
}

SimulationDesReponses = function(w,vraieloi)
{
  a = dim(w)
  n=a[1]
  y = matrix(0,nrow=n,ncol=1,byrow=T)   
  for (i in 1:n) 
  {
    y[i]= 2*(runif(1,0,1) < 1/(1 + exp(-w[i,]%*%t(vraieloi)))) - 1
  }
  return(y)
}

#ContaminationResponse = function(w,y,vraieloi,thetaconta,eps)
#{
#  a = dim(w)
#  n = a[1]
#  r = matrix(0,nrow=n,ncol=1,byrow=T)   
#  for (i in 1:n) 
#  { 
#    p = runif(1,0,1)
#    if(p<=eps)
#    r[i]= 2*(runif(1,0,1) < 1/(1 + exp(-w[i,]%*%t(thetaconta)))) - 1
#    else
#    r[i]= y[i]
#  }
#  return(r)
#}
###########################################compute function Rho##########################
chi = function(x) 
{return ((x-1)/(x+1))}

chigrand = function(theta,thetaprime) {
  y=0
  if (thetaprime<0 && theta<0) {y=chi(sqrt((1+exp(thetaprime))/(1+exp(theta))))}
  else {
    if (thetaprime<0 && theta>=0) {
      if (theta<=300) {y=chi(sqrt((1+exp(thetaprime))/(1+exp(theta))))}
      else {y=-1}
    }
    else {
      if (thetaprime>=0 && theta<0) {
        if (thetaprime>300) {y=1}
        else {y=chi(sqrt((1+exp(thetaprime))/(1+exp(theta))))}
      }
      else {
        if (thetaprime<=200) {
          if (theta<=300) {y=chi(sqrt((1+exp(thetaprime))/(1+exp(theta))))}
          else {y=-1}
        }
        else {
          if (theta<=100) {y=1}
          else {
            if (0.5*(thetaprime-theta) >50) {y=1}
            else {y=chi(exp(0.5*(thetaprime-theta)))}
          }
        }
      }
    }
  }
  return(y)
}

MoinsTStatistiqueNew = function(w,y,theta,thetaprime) 
{
  a = dim(w)
  n=a[1]
  r=0
  for (i in 1:n)
  { 
    a1= w[i,]%*%(as.numeric(t(theta)))
    a = -a1*y[i]
    b1= w[i,]%*%(as.numeric(t(thetaprime)))
    b = -b1*y[i]
    r = r + chigrand(a,b)
  }
  return(r)
}

#acbias = function (w,y,theta){
#  s = 0
#  a = dim(w)
#  n = a[1]
#  ytheoretical = matrix(0,nrow=n,ncol=1)
#  for (i in 1:n) {
#    a1=w[i,]%*%theta
#    if(a1>0) ytheoretical[i]=1
#    else ytheoretical[i]=-1
#    if(y[i]!=ytheoretical[i])
#      s = s+1
#  }
#  return(s)
#}

#choosestarter = function(w,y){
#  k=dim(w)[2]
#  start = c(0,0,0,0,0,0)
#  objective = function(p) {
#    q = acbias(w,y,p)
#    return(q)}
#  b=cma_es(start,objective)
#  c=t(b$par)
#  return(c)
#}
choosestarter = function(x,y){
  xp = x
  class = as.numeric(y)
  training_data = data.frame(xp = xp, class = as.factor(class))
  svmfit = svm(class~., 
               data = training_data,
               kernel = "linear", 
               cost = 10,
               scale = FALSE)
  c=t(svmfit$coefs)%*%svmfit$SV
  b=-svmfit$rho
  d=cbind(b,c)
  if((as.numeric(svmfit$labels)[1])==1)
  {d = -d}
  return(d)
}

DescenteDeRhoEstimationNew = function(W,Y,longeurTheta)
{
  count = 0
  start = as.numeric(t(choosestarter(W[,c(-1)],Y)))
  a = start
  for (m in 1:longeurTheta) {
    objective = function(x) {
      y = MoinsTStatistiqueNew(W,Y,t(a),t(x))
      return(y)}
    b = cma_es(start,objective)
    c = t(b$par)
    count = count+1
    # stop condition
    if (MoinsTStatistiqueNew(W,Y,a,c) > -1){
      a=as.numeric(t(c))
      break
    }
    a=as.numeric(t(c))
  }
  newlist=list(first=t(a),second=count)
  return(newlist)
}
###################################compute function MLE####################################
MLE = function(x,y){
  objective = function(theta1,theta2,theta3,theta4,theta5,theta6){
    s = 0
    a = dim(x)
    n = a[1]
    for (i in 1:n) {
      b = x[i,1]*theta1+x[i,2]*theta2+x[i,3]*theta3+x[i,4]*theta4+x[i,5]*theta5+x[i,6]*theta6
      t = log(1+exp(-y[i]*b))
      s = s + t
    }
    return(s)
  }
  mlestimator = stats4::mle(minuslogl=objective,start=list(theta1=0,theta2=0,theta3=0,theta4=0,theta5=0,theta6=0)) 
  mle = t(as.numeric(coef(mlestimator)))
  return(mle)
}
####################################compute function Hellinger-squared####################
DistanceEmpirique = function (w,vraieloi,loiestimer){
  s = 0;
  a = dim(w)
  n = a[1]
  for (i in 1:n) {
    x = 1 - (1/sqrt(((1+exp(-w[i,]%*%t(vraieloi)))*(1+exp(-w[i,]%*%t(loiestimer)))))) - (1/sqrt(((1+exp(w[i,]%*%t(vraieloi)))*(1+exp(w[i,]%*%t(loiestimer))))))
    s = s + x
  }
  s = s/n  
  return(s)
}

#HellingerDistance = function (w,vstar,vprime,u,eps){
#  s = 0;
#  a = dim(w)
#  n = a[1]
#  for (i in 1:n) {
#    p1 = (1-eps)/(1+exp(-w[i,]%*%t(vstar)))+eps/(1+exp(-w[i,]%*%t(vprime)))
#    p0 = (1-eps)/(1+exp(w[i,]%*%t(vstar)))+eps/(1+exp(w[i,]%*%t(vprime)))
#    x = 1 - sqrt(p1)*(1/sqrt((1+exp(-w[i,]%*%t(u))))) - sqrt(p0)*(1/sqrt((1+exp(w[i,]%*%t(u)))))
#    s = s + x
#  }
#  s = s/n  
#  return(s)
#}
###################################100 times experiments################################
longeurTheta = 100
#epsilon = 0.05
l=1
LOL = matrix(0,nrow = 100,ncol = 8,byrow=T)
for (l in 1:100)
{
  W = SimulationDesCaracteristiques(n,k)
  Y = SimulationDesReponses(W,vraieloi) 
  Wmonte = SimulationDesCaracteristiques(1000,k)
  #Yconta = ContaminationResponse(W,Y,vraieloi,thetaconta,epsilon)
  WT=matrix(0,nrow=n+1,ncol=k+1,byrow=T)
  for (i in 1:n)
  {WT[i,]=W[i,]}
  WT[n+1,]=cbind(1,1000*matrix(1,nrow=1,ncol=k,byrow=T))
  YT=matrix(0,nrow=n+1,ncol=1,byrow=T)
  for (i in 1:n)
  {YT[i]=Y[i]}
  YT[n+1]=-1
  loi1 = DescenteDeRhoEstimationNew(W,Y,longeurTheta)
  loiestimer1 = loi1$first
  count1 = loi1$second
  startpoint1 = t(as.numeric(t(choosestarter(W[,c(-1)],Y))))
  mlestimator1 = MLE(W,Y)
  HelDisStart1 = DistanceEmpirique(Wmonte,vraieloi,startpoint1)
  HelDisRho1 = DistanceEmpirique(Wmonte,vraieloi,loiestimer1)
  HelDisMle1 = DistanceEmpirique(Wmonte,vraieloi,mlestimator1)
  loi2 = DescenteDeRhoEstimationNew(WT,YT,longeurTheta)
  loiestimer2 = loi2$first
  count2 = loi2$second
  startpoint2 = t(as.numeric(t(choosestarter(WT[,c(-1)],YT))))
  mlestimator2 = MLE(WT,YT)
  HelDisStart2 = DistanceEmpirique(Wmonte,vraieloi,startpoint2)
  HelDisRho2 = DistanceEmpirique(Wmonte,vraieloi,loiestimer2)
  HelDisMle2 = DistanceEmpirique(Wmonte,vraieloi,mlestimator2)
  #loi3 = DescenteDeRhoEstimationNew(W,Yconta,longeurTheta)
  #loiestimer3 = loi3$first
  #count3 = loi3$second
  #mlestimator3 = MLE(W,Yconta)
  #HelDisRho3 = HellingerDistance(Wmonte,vraieloi,thetaconta,loiestimer3,epsilon)
  #HelDisTTT = DistanceEmpirique(Wmonte,vraieloi,thetaconta)
  #HelTTT = HellingerDistance(Wmonte,vraieloi,thetaconta,thetaconta,epsilon)
  #HelDis = DistanceEmpirique(Wmonte,vraieloi,mlestimator3)
  #HelDisMle3 = HellingerDistance(Wmonte,vraieloi,thetaconta,mlestimator3,epsilon)
  LOL[l,1]=as.numeric(format(HelDisStart1,digits = 3))
  LOL[l,2]=as.numeric(format(HelDisRho1,digits = 3))
  LOL[l,3]=as.numeric(format(count1,digits = 1))
  LOL[l,4]=as.numeric(format(HelDisMle1,digits = 3))
  LOL[l,5]=as.numeric(format(HelDisStart2,digits = 3))
  LOL[l,6]=as.numeric(format(HelDisRho2,digits = 3))
  LOL[l,7]=as.numeric(format(count2,digits = 1))
  LOL[l,8]=as.numeric(format(HelDisMle2,digits = 3))
  #LOL[l,7]=as.numeric(format(HelDisRho3,digits = 3))
  #LOL[l,8]=as.numeric(format(count3,digits = 1))
  #LOL[l,9]=as.numeric(format(HelDisMle3,digits = 3))
  l=l+1
}
logisticresult = data.frame(No.=c(1:100), startpoint1=LOL[,1], loiestimer1=LOL[,2], count1=LOL[,3], mlestimator1=LOL[,4], startpoint2=LOL[,5], loiestimer2=LOL[,6], count2=LOL[,7], mlestimator2=LOL[,8])
write.csv(logisticresult,file = "log.csv")  

