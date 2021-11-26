######################################data simulation#########################################
install.packages("cmaes",repos=structure(c(CRAN="https://packages.othr.de/cran/")))
library(stats4)
library(cmaes)
k = 6
n = 500
realtheta=t(c(0.07,3,4,6,2,1))


generateX=function(n,k){
  x = matrix(0,nrow=n,ncol=(k-1))
  constcol= matrix(1,nrow=n,ncol=1)
  for (i in 1:n)
  {
    x[i,1]=runif(1,min=0,max=0.01)
    x[i,2]=runif(1,min=0,max=0.01)
    x[i,3]=runif(1,min=0,max=0.01)
    x[i,4]=runif(1,min=0,max=0.1)
    x[i,5]=runif(1,min=0,max=0.1)
  }
  x = cbind(constcol,x)
  return(x)
}

generateY=function(x,v){
  n=dim(x)[1]
  y = matrix(0,nrow=n,ncol=1)  
  lambda = matrix(0,nrow=n,ncol=1)
  for (i in 1:n)
  {
    lambda[i]=log(1+exp(x[i,]%*%t(v)))
    y[i]=rexp(1,lambda[i])
  }
  return(y)
}

generateYconta = function(x,y,realtheta,eps)
{
  a = dim(x)
  n = a[1]
  r = matrix(0,nrow=n,ncol=1,byrow=T)   
  for (i in 1:n) 
  { 
    p = runif(1,0,1)
    if(p<=eps){
      r[i] = runif(1,50,60)
    }
    else{
      r[i] = y[i] 
    }
  }
  return(r)
}
###########################################compute function Rho##########################
chi = function(x) 
{return ((x-1)/(x+1))}

chigrand = function(theta,thetaprime,b)
{ 
  y=-chi(sqrt(exp((theta-thetaprime)*b)*(thetaprime/theta)))
  return(y)
}

MoinsTStatistiqueNew = function(x,y,theta,thetaprime) 
{ 
  a = dim(x)
  ndata = a[1]
  sumchi = 0
  for (i in 1:ndata)
  { 
    d1 = x[i,]%*%(as.numeric(t(theta)))
    d2 = x[i,]%*%(as.numeric(t(thetaprime)))
    a1 = log(1+exp(d1))
    b1 = log(1+exp(d2))
    sumchi = sumchi + chigrand(a1,b1,y[i])
  }
  return(sumchi)
}

acbias = function (x,y,theta){
  s = 0
  a = dim(x)
  n = a[1]
  for (i in 1:n) {
    a1 = log(1+exp(x[i,]%*%theta))
    b1 = log(2)/a1
    t = abs(y[i]-b1)
    s = s+t
  }
  return(s/n)
}

choosestarter = function(x,y){
  k=dim(x)[2]
  start = c(0,0,0,0,0,0)
  objective = function(p){
    q = acbias(x,y,p)
    return(q)}
  b = cma_es(start,objective)
  c = t(b$par)
  return(c)
}


DescenteDeRhoEstimationNew = function(W,Y,longeurTheta)
{
  count = 0
  start = as.numeric(t(choosestarter(W,Y)))
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
      d = log(1+exp(b))
      t = log(d)-d*y[i]
      s = s-t
    }
    return(s)
  }
  mlestimator = stats4::mle(minuslogl=objective,start=list(theta1=0,theta2=0,theta3=0,theta4=0,theta5=0,theta6=0)) 
  mle = t(as.numeric(coef(mlestimator)))
  return(mle)
}
####################################compute function Hellinger-squared####################
DistanceEmpirique = function (x,v,u){
  y = 0;
  a = dim(x)
  n = a[1]
  for (i in 1:n) {
    theta1 = log(1+exp(x[i,]%*%t(v)))
    theta2 = log(1+exp(x[i,]%*%t(u)))
    s = 1 - 2*sqrt(theta1*theta2)/(theta1+theta2)
    y = y+s
  }
  y = y/n  
  return(y)
}

HellingerDistance = function (w,vstar,u,eps){
  s = 0;
  a = dim(w)
  n = a[1]
  for (i in 1:n) {
    x = 1
    a1 = log(1+exp(w[i,]%*%t(vstar)))
    c1 = log(1+exp(w[i,]%*%t(u)))
    for (l in 1:10000){
      step = 100/10000
      y = step*l
      p1 = a1*exp(-a1*y)
      p2 = c1*exp(-c1*y)
      if((l>=5001)&&(l<=6000))
      {p3 = 0.1}
      else
      {p3 = 0}
      x = x - sqrt(((1-eps)*p1+eps*p3)*p2)*step
    }
    s = s + x
  }
  s = s/n  
  return(s)
}
###################################100 times experiments################################
longeurTheta = 100
epsilon = 0.05
l=1
LOL = matrix(0,nrow = 100,ncol = 69,byrow=T)
for (l in 1:100)
{
  X = generateX(n,k)
  Y = generateY(X,realtheta) 
  Xmonte = generateX(1000,k)
  Xout=t(c(1,0.005,0.005,0.005,0.05,0.05))
  XT=rbind(X,Xout)
  YT=rbind(Y,1000)
  Yconta = generateYconta(X,Y,realtheta,epsilon)
  loi1 = DescenteDeRhoEstimationNew(X,Y,longeurTheta)
  loiestimer1 = loi1$first
  count1 = loi1$second
  startpoint1 = t(as.numeric(t(choosestarter(X,Y))))
  mlestimator1 = MLE(X,Y)
  HelDisStart1 = DistanceEmpirique(Xmonte,realtheta,startpoint1)
  HelDisRho1 = DistanceEmpirique(Xmonte,realtheta,loiestimer1)
  HelDisMle1 = DistanceEmpirique(Xmonte,realtheta,mlestimator1)
  loi2 = DescenteDeRhoEstimationNew(XT,YT,longeurTheta)
  loiestimer2 = loi2$first
  count2 = loi2$second
  mlestimator2 = MLE(XT,YT) 
  startpoint2 = t(as.numeric(t(choosestarter(XT,YT))))
  HelDisStart2 = DistanceEmpirique(Xmonte,realtheta,startpoint2)
  HelDisRho2 = DistanceEmpirique(Xmonte,realtheta,loiestimer2)
  HelDisMle2 = DistanceEmpirique(Xmonte,realtheta,mlestimator2)
  loi3 = DescenteDeRhoEstimationNew(X,Yconta,longeurTheta)
  loiestimer3 = loi3$first
  count3 = loi3$second
  mlestimator3 = MLE(X,Yconta)
  startpoint3 = t(as.numeric(t(choosestarter(X,Yconta))))
  HelDisStart3star = DistanceEmpirique(Xmonte,realtheta,startpoint3)
  HelDisStart3con = HellingerDistance(Xmonte,realtheta,startpoint3,epsilon)
  HelDisRho3con = HellingerDistance(Xmonte,realtheta,loiestimer3,epsilon)
  HelDisMle3con = HellingerDistance(Xmonte,realtheta,mlestimator3,epsilon)
  HelDisRho3star = DistanceEmpirique(Xmonte,realtheta,loiestimer3)
  HelDisMle3star = DistanceEmpirique(Xmonte,realtheta,mlestimator3)
  LOL[l,1]=as.numeric(format(HelDisStart1,digits = 3))
  LOL[l,2]=as.numeric(format(HelDisRho1,digits = 3))
  LOL[l,3]=as.numeric(format(count1,digits = 1))
  LOL[l,4]=as.numeric(format(HelDisMle1,digits = 3))
  LOL[l,5]=as.numeric(format(HelDisStart2,digits = 3))
  LOL[l,6]=as.numeric(format(HelDisRho2,digits = 3))
  LOL[l,7]=as.numeric(format(count2,digits = 1))
  LOL[l,8]=as.numeric(format(HelDisMle2,digits = 3))
  LOL[l,9]=as.numeric(format(HelDisStart3con,digits = 3))
  LOL[l,10]=as.numeric(format(HelDisRho3con,digits = 3))
  LOL[l,11]=as.numeric(format(count3,digits = 1))
  LOL[l,12]=as.numeric(format(HelDisMle3con,digits = 3))
  LOL[l,13]=as.numeric(format(HelDisStart3star,digits = 3))
  LOL[l,14]=as.numeric(format(HelDisRho3star,digits = 3))
  LOL[l,15]=as.numeric(format(HelDisMle3star,digits = 3))
  LOL[l,16]=as.numeric(format(loiestimer3[1],digits = 5))
  LOL[l,17]=as.numeric(format(loiestimer3[2],digits = 5))
  LOL[l,18]=as.numeric(format(loiestimer3[3],digits = 5))
  LOL[l,19]=as.numeric(format(loiestimer3[4],digits = 5))
  LOL[l,20]=as.numeric(format(loiestimer3[5],digits = 5))
  LOL[l,21]=as.numeric(format(loiestimer3[6],digits = 5))
  LOL[l,22]=as.numeric(format(loiestimer2[1],digits = 5))
  LOL[l,23]=as.numeric(format(loiestimer2[2],digits = 5))
  LOL[l,24]=as.numeric(format(loiestimer2[3],digits = 5))
  LOL[l,25]=as.numeric(format(loiestimer2[4],digits = 5))
  LOL[l,26]=as.numeric(format(loiestimer2[5],digits = 5))
  LOL[l,27]=as.numeric(format(loiestimer2[6],digits = 5))
  LOL[l,28]=as.numeric(format(loiestimer1[1],digits = 5))
  LOL[l,29]=as.numeric(format(loiestimer1[2],digits = 5))
  LOL[l,30]=as.numeric(format(loiestimer1[3],digits = 5))
  LOL[l,31]=as.numeric(format(loiestimer1[4],digits = 5))
  LOL[l,32]=as.numeric(format(loiestimer1[5],digits = 5))
  LOL[l,33]=as.numeric(format(loiestimer1[6],digits = 5))
  LOL[l,34]=as.numeric(format(mlestimator3[1],digits = 5))
  LOL[l,35]=as.numeric(format(mlestimator3[2],digits = 5))
  LOL[l,36]=as.numeric(format(mlestimator3[3],digits = 5))
  LOL[l,37]=as.numeric(format(mlestimator3[4],digits = 5))
  LOL[l,38]=as.numeric(format(mlestimator3[5],digits = 5))
  LOL[l,39]=as.numeric(format(mlestimator3[6],digits = 5))
  LOL[l,40]=as.numeric(format(mlestimator2[1],digits = 5))
  LOL[l,41]=as.numeric(format(mlestimator2[2],digits = 5))
  LOL[l,42]=as.numeric(format(mlestimator2[3],digits = 5))
  LOL[l,43]=as.numeric(format(mlestimator2[4],digits = 5))
  LOL[l,44]=as.numeric(format(mlestimator2[5],digits = 5))
  LOL[l,45]=as.numeric(format(mlestimator2[6],digits = 5))
  LOL[l,46]=as.numeric(format(mlestimator1[1],digits = 5))
  LOL[l,47]=as.numeric(format(mlestimator1[2],digits = 5))
  LOL[l,48]=as.numeric(format(mlestimator1[3],digits = 5))
  LOL[l,49]=as.numeric(format(mlestimator1[4],digits = 5))
  LOL[l,50]=as.numeric(format(mlestimator1[5],digits = 5))
  LOL[l,51]=as.numeric(format(mlestimator1[6],digits = 5))
  LOL[l,52]=as.numeric(format(startpoint3[1],digits = 5))
  LOL[l,53]=as.numeric(format(startpoint3[2],digits = 5))
  LOL[l,54]=as.numeric(format(startpoint3[3],digits = 5))
  LOL[l,55]=as.numeric(format(startpoint3[4],digits = 5))
  LOL[l,56]=as.numeric(format(startpoint3[5],digits = 5))
  LOL[l,57]=as.numeric(format(startpoint3[6],digits = 5))
  LOL[l,58]=as.numeric(format(startpoint2[1],digits = 5))
  LOL[l,59]=as.numeric(format(startpoint2[2],digits = 5))
  LOL[l,60]=as.numeric(format(startpoint2[3],digits = 5))
  LOL[l,61]=as.numeric(format(startpoint2[4],digits = 5))
  LOL[l,62]=as.numeric(format(startpoint2[5],digits = 5))
  LOL[l,63]=as.numeric(format(startpoint2[6],digits = 5))
  LOL[l,64]=as.numeric(format(startpoint1[1],digits = 5))
  LOL[l,65]=as.numeric(format(startpoint1[2],digits = 5))
  LOL[l,66]=as.numeric(format(startpoint1[3],digits = 5))
  LOL[l,67]=as.numeric(format(startpoint1[4],digits = 5))
  LOL[l,68]=as.numeric(format(startpoint1[5],digits = 5))
  LOL[l,69]=as.numeric(format(startpoint1[6],digits = 5))
  l=l+1
}
expresult = data.frame(No.=c(1:100), startpoint1=LOL[,1], loiestimer1=LOL[,2], count1=LOL[,3], mlestimator1=LOL[,4], startpoint2=LOL[,5], loiestimer2=LOL[,6], count2=LOL[,7], mlestimator2=LOL[,8], startpoint3con=LOL[,9], loiestimer3con=LOL[,10],count3=LOL[,11], mlestimator3con=LOL[,12], startpoint3star=LOL[,13], loiestimer3star=LOL[,14],mlestimator3star=LOL[,15]
                       ,loiestimer31=LOL[,16],loiestimer32=LOL[,17],loiestimer33=LOL[,18],loiestimer34=LOL[,19],loiestimer35=LOL[,20],loiestimer36=LOL[,21],loiestimer21=LOL[,22],loiestimer22=LOL[,23],loiestimer23=LOL[,24],loiestimer24=LOL[,25],loiestimer25=LOL[,26],loiestimer26=LOL[,27],loiestimer11=LOL[,28],loiestimer12=LOL[,29],loiestimer13=LOL[,30],loiestimer14=LOL[,31],loiestimer15=LOL[,32],loiestimer16=LOL[,33],mlestimator31=LOL[,34]
                       ,mlestimator32=LOL[,35],mlestimator33=LOL[,36],mlestimator34=LOL[,37],mlestimator35=LOL[,38],mlestimator36=LOL[,39],mlestimator21=LOL[,40],mlestimator22=LOL[,41],mlestimator23=LOL[,42],mlestimator24=LOL[,43],mlestimator25=LOL[,44],mlestimator26=LOL[,45],mlestimator11=LOL[,46],mlestimator12=LOL[,47],mlestimator13=LOL[,48],mlestimator14=LOL[,49],mlestimator15=LOL[,50],mlestimator16=LOL[,51],startpoint31=LOL[,52],startpoint32=LOL[,53],startpoint33=LOL[,54],startpoint34=LOL[,55],startpoint35=LOL[,56],startpoint36=LOL[,57]
                       ,startpoint21=LOL[,58],startpoint22=LOL[,59],startpoint23=LOL[,60],startpoint24=LOL[,61],startpoint25=LOL[,62],startpoint26=LOL[,63],startpoint11=LOL[,64],startpoint12=LOL[,65],startpoint13=LOL[,66],startpoint14=LOL[,67],startpoint15=LOL[,68],startpoint16=LOL[,69])
write.csv(expresult,file = "exp.csv")

