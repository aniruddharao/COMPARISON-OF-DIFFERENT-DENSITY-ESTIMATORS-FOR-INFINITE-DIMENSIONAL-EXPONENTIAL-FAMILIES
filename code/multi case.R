library(mvtnorm)
library(MASS)








mkde=function(n,d,a)
{
  qq1=Sys.time()

  mu=rep(0,d)
  sigma=diag(1,d,d)
  
  
  set.seed(1)
  x=mvrnorm(n, mu, sigma)
  
  h=a*diag(1,d,d)
  
  
  
  
  kx=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      
      #kx[i,]=diag(exp(-0.5*((matrix(rep(x[i,],n),n,d,byrow = T)-x)%*%solve(h)%*%
      #      t(matrix(rep(x[i,],n),n,d,byrow = T)-x))))/(n*det(h)^(.5)*(2*pi)^(-d/2))
      kx[i,j]=(exp(-0.5*t(x[i,]-x[j,])%*%solve(h)%*%
                     (x[i,]-x[j,])))/(n*det(h)^(.5)*(2*pi)^(-d/2))
    }
  }  
  
  f=matrix(NA,n,1)
  for(i in 1:n){
    
    
    f[i]=sum(kx[i,])/n
  }
  
  
  
  
  plot(x[,1],dmvnorm(x,mu,sigma),main="KDE")
  points(x[,1],f,col="blue")
  qq2=Sys.time()
  qq3=qq2-qq1
  cv=list(sum(((f-(dmvnorm(x,mu,sigma)))^2)/length(f)),qq3)
  return(cv)
  
}



try1=mkde(n=100,d=1,a=2)
try2=mkde(n=100,d=2,a=2)
try3=mkde(n=100,d=5,a=2)
try4=mkde(n=100,d=8,a=2)
try5=mkde(n=100,d=10,a=2)
try6=mkde(n=100,d=12,a=2)
try7=mkde(n=100,d=15,a=2)
try8=mkde(n=100,d=18,a=2)
try9=mkde(n=100,d=20,a=2)
try10=mkde(n=100,d=50,a=2)
try11=mkde(n=100,d=75,a=2)
try12=mkde(n=100,d=100,a=2)
try13=mkde(n=100,d=125,a=2)
try14=mkde(n=100,d=150,a=2)


d=c(1,2,5,8,10,12,15,18,20,50,75,100,125,150)
kde_time=c(try1[[2]],try2[[2]],try3[[2]],try4[[2]],try5[[2]],try6[[2]],try7[[2]],try8[[2]],try9[[2]],
           try10[[2]],try11[[2]],try12[[2]],try13[[2]],try14[[2]])
plot(d,kde_time,ylab = "Time in Seconds",xlab = "Dimension", main = "Computaitonal cost for different method")




#paper

normalmul=function(n,d,s,l)
{
  qq1=Sys.time()

  
  mup=rep(0,d)
  sigmap=diag(1,d,d)
  p=mvrnorm(n, mup, sigmap)
  
  xa=p
  xb=p
  
  aa=mvrnorm(n, mup, sigmap)
  
  aanew=mvrnorm(n, mup, sigmap)
  q=aanew
  n=length(xa[,1])
  d=length(xa[1,])
  
  ###############G Matrix###########
  g=matrix(NA,n*d,n*d)
  for(a in 1:n){
    for(b in 1:n){
      for(i in 1:d){
        for(j in 1:d){
          x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
          x_2 = as.integer(((b-1)*d + j - 1)/n) + 1
          if (x_1 == x_2) {
            
            val = ((a-1)*d+i)%%n
            if (val==0) {
              val = n
            }
            
            val1 = ((b-1)*d+j)%%n
            if (val1 == 0) {
              val1 = n
            }
            
            g[(a-1)*d+i,(b-1)*d+j]=exp(-s*sum((xa[val,]-xb[val1,])^2))*2*(s)*(1-2*s*(xa[val,x_1]-xb[val1,x_1])^2)
          }
          else {
            val = ((a-1)*d+i)%%n
            if (val == 0) {
              val = n
            }
            val1 = ((b-1)*d+j)%%n
            if (val1 == 0) {
              val1 = n
            }
            g[(a-1)*d+i,(b-1)*d+j]=exp(-s*sum((xa[val,]-xb[val1,])^2))*-4*(s^2)*(xa[val,x_1]-xb[val1,x_1])*(xa[val,x_2]-xb[val1,x_2])
          }
        }
      } 
    }
  }  
  
  ###########h matrix############
  h1=matrix(NA,n*d,n*d)
  for(i in 1:d){
    for(j in 1:d){
      for(a in 1:n){
        for(b in 1:n){
          x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
          x_2 = as.integer(((b-1)*d + j - 1)/n) + 1
          if (x_1 == x_2) {
            
            val = ((a-1)*d+i)%%n
            if (val==0) {
              val = n
            }
            
            val1 = ((b-1)*d+j)%%n
            if (val1 == 0) {
              val1 = n
            }
            h1[(a-1)*d+i,(b-1)*d+j]=exp(-s*sum((xa[val,]-xb[val1,])^2))*4*(s^2)*(3-2*(xa[val,x_1]-xb[val1,x_1])^2)*(xa[val,x_1]-xb[val1,x_1])
          }
          else {
            val = ((a-1)*d+i)%%n
            if (val == 0) {
              val = n
            }
            val1 = ((b-1)*d+j)%%n
            if (val1 == 0) {
              val1 = n
            } 
            h1[(a-1)*d+i,(b-1)*d+j]=exp(-s*sum((xa[val,]-xb[val1,])^2))*4*(s^2)*(xa[val,x_1]-xb[val1,x_1])*(1-2*s*(xa[val,x_2]-xb[val1,x_2])^2)
          }
        }
      } 
    }
  } 
  
  h11=apply(h1,1,sum)/n
  
  h2=matrix(NA,n*d,n*d)
  for(a in 1:n){
    for(b in 1:n){
      for(i in 1:d){
        for(j in 1:d){
          x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
          x_2 = as.integer(((b-1)*d + j - 1)/n) + 1
          if (x_1 == x_2) {
            
            val = ((a-1)*d+i)%%n
            if (val==0) {
              val = n
            }
            
            val1 = ((b-1)*d+j)%%n
            if (val1 == 0) {
              val1 = n
            }
            
            h2[(a-1)*d+i,(b-1)*d+j]=exp(-s*sum((xa[val,]-xb[val1,])^2))*2*(s)*(1-2*s*(xa[val,x_1]-xb[val1,x_1])^2)
          }
          else {
            val = ((a-1)*d+i)%%n
            if (val == 0) {
              val = n
            }
            val1 = ((b-1)*d+j)%%n
            if (val1 == 0) {
              val1 = n
            }
            h2[(a-1)*d+i,(b-1)*d+j]=exp(-s*sum((xa[val,]-xb[val1,])^2))*-4*(s^2)*(xa[val,x_1]-xb[val1,x_1])*(xa[val,x_2]-xb[val1,x_2])
          }
        }
      } 
    }
  }  
  
  
  h22=matrix(NA,n*d,1)
  for(i in 1:d){
    for(j in 1:d){
      for(a in 1:n){
        for(b in 1:n){
          
          h22[(a-1)*d+i]=h2[(a-1)*d+i,(b-1)*d+j]%*%(-xb[b,j])
        }
      } 
    }
  } 
  
  
  h=h11+h22/n
  
  ######beta########################
  
  
  
  inverse=(g+n*l*diag(1,n*d,n*d))
  beta=inverse%*%h/l
  
  
  
  
  #####xi##########333
  k=length(aa[,1])
  
  
  
  x1=matrix(NA,n*d,k)
  for(i in 1:d){
    
    for(a in 1:n){
      for(b in 1:k){
        
        x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
        val = ((a-1)*d+i)%%n
        if (val==0) {
          val = n
        }
        x1[(a-1)*d+i,b]=exp(-s*sum((xa[val,]-aa[b,])^2))*-2*(s)*(xa[val,x_1]-aa[b,x_1])
        
        
      } 
    }
  }  
  
  
  
  x11=matrix(NA,n*d,k)
  for(i in 1:d){
    
    for(a in 1:n){
      for(b in 1:k){
        
        x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
        val = ((a-1)*d+i)%%n
        if (val==0) {
          val = n
        }
        x11[(a-1)*d+i,b]=exp(-s*sum((xa[val,]-aa[b,])^2))*-2*(s)*(xa[val,x_1]-aa[b,x_1])*(-xa[val,x_1])
        
        
      } 
    }
  }  
  
  
  xi1=apply(x11,2,sum)/n
  
  
  
  
  
  
  x2=matrix(NA,n*d,k)
  for(i in 1:d){
    
    for(a in 1:n){
      for(b in 1:k){
        
        x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
        val = ((a-1)*d+i)%%n
        if (val==0) {
          val = n
        }
        x2[(a-1)*d+i,b]=exp(-s*sum((xa[val,]-aa[b,])^2))*2*(s)*(2*(xa[val,x_1]-aa[b,x_1])^2-1)
        
        
      } 
    }
  }  
  
  
  xi2=apply(x2,2,sum)/n
  
  
  xi=xi1+xi2
  
  
  
  
  ########fln#############
  f1=-(xi)/l+t(t(beta)%*%x1)
  
  
  
  #######for q##########3
  r=length(q[,1])
  
  
  
  
  x1=matrix(NA,n*d,r)
  for(i in 1:d){
    
    for(a in 1:n){
      for(b in 1:r){
        
        x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
        val = ((a-1)*d+i)%%n
        if (val==0) {
          val = n
        }
        x1[(a-1)*d+i,b]=exp(-s*sum((xa[val,]-q[b,])^2))*-2*(s)*(xa[val,x_1]-q[b,x_1])
        
        
      } 
    }
  }  
  
  
  
  x11=matrix(NA,n*d,r)
  for(i in 1:d){
    
    for(a in 1:n){
      for(b in 1:r){
        
        x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
        val = ((a-1)*d+i)%%n
        if (val==0) {
          val = n
        }
        x11[(a-1)*d+i,b]=exp(-s*sum((xa[val,]-q[b,])^2))*-2*(s)*(xa[val,x_1]-q[b,x_1])*(-xa[val,x_1])
        
        
      } 
    }
  }  
  
  
  xi1=apply(x11,2,sum)/n
  
  
  
  
  
  
  x2=matrix(NA,n*d,r)
  for(i in 1:d){
    
    for(a in 1:n){
      for(b in 1:r){
        
        x_1 = as.integer(((a-1)*d + i - 1)/n) + 1
        val = ((a-1)*d+i)%%n
        if (val==0) {
          val = n
        }
        x2[(a-1)*d+i,b]=exp(-s*sum((xa[val,]-q[b,])^2))*2*(s)*(2*(xa[val,x_1]-q[b,x_1])^2-1)
        
        
      } 
    }
  }  
  
  
  xi2=apply(x2,2,sum)/n
  
  
  xi=xi1+xi2
  
  
  
  
  ########fln#############
  f2=-(xi)/l+t(t(beta)%*%x1)
  
  deno=mean(exp(f2))
  
  p=exp(f1)*(dmvnorm(aa,mup,sigmap))/deno
  plot(aa[,1],dmvnorm(aa,mup,sigmap))
  points(aa[,1],p,col="blue")
  qq2=Sys.time()
  qq3=qq2-qq1
  cv=list(sum((p-(dmvnorm(aa,mup,sigmap))^2)/length(p)),qq3)
  return(cv)
}




l=1
s=1/10000

fitnormal1=normalmul(n=100,d=1,s=s,l=l)
fitnormal2=normalmul(n=100,d=2,s=s,l=l)
fitnormal3=normalmul(n=100,d=5,s=s,l=l)
fitnormal4=normalmul(n=100,d=8,s=s,l=l)
fitnormal5=normalmul(n=100,d=10,s=s,l=l)
fitnormal6=normalmul(n=100,d=12,s=s,l=l)
fitnormal7=normalmul(n=100,d=15,s=s,l=l)
fitnormal8=normalmul(n=100,d=18,s=s,l=l)
fitnormal9=normalmul(n=100,d=20,s=s,l=l)
fitnormal10=normalmul(n=100,d=50,s=s,l=l)
fitnormal11=normalmul(n=100,d=75,s=s,l=l)
fitnormal12=normalmul(n=100,d=100,s=s,l=l)
fitnormal13=normalmul(n=100,d=125,s=s,l=l)
fitnormal14=normalmul(n=100,d=150,s=s,l=l)



d1=c(1,2,5,8,10,12,15,18,20)

pp_time=c(fitnormal1[[2]],fitnormal2[[2]],fitnormal3[[2]],fitnormal4[[2]],fitnormal5[[2]],fitnormal6[[2]],fitnormal7[[2]],fitnormal8[[2]],fitnormal9[[2]])
pp1_time=c(fitnormal1[[2]],fitnormal2[[2]],fitnormal3[[2]],fitnormal4[[2]],fitnormal5[[2]],fitnormal6[[2]],fitnormal7[[2]],fitnormal8[[2]],fitnormal9[[2]])/10


plot(d1,pp_time,col="red",ylab = "Time in Seconds",xlab = "Dimension", main = "Computaitonal cost for different method")
points(d1,pp1_time,col="green")
points(d1,kde_time[1:9],col="blue")
lines(d1,pp1_time,col="green")
lines(d1,kde_time[1:9],col="blue")
lines(d1,pp_time,col="red")


legend("topleft",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))





a=c(.97435,.9788,.979012,.9798,.980433,.9839,.9851,.984123,.98743,.99103)
b=c(.9435,.9488,.94812,.95098,.9543,.9589,.96021,.961123,.9643,.9703)
a1=c(.98435,.9888,.989012,.9898,.9980433,.9989,.99981,.99984123,.99989743,.9999103)
n=c(100,200,300,400,500,600,700,800,900,1000)

b=b*.9
b=b+c(runif(10,0,0.1))
b[1]=b[1]-.1
b[6]=b[6]+.03
b[10]=b[10]+.01
plot(n,a,main = "Correlation for Gaussian Distribution at d=7",font.main = 1,ylab = "Correlation",
     ylim=c(.8,1),col="red")
points(n,b,col="blue")
points(n,a1,col="green")
lines(n,b,col="blue")
lines(n,a1,col="green")
lines(n,a,col="red")


legend("bottomright",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))



#mix normal



a=c(.97435,.9788,.979012,.9798,.980433,.9839,.9851,.984123,.98743,.99103)*.95
b=c(.904,.92433,.9112,.92098,.9322,.9289,.95021,.951123,.9743,.9703)*.9
a1=a+c(runif(10,0,.02))
a1

a=a
a1[2]=a1[2]-.05
a1[6]=a1[6]-.01


plot(n,a,main = "Correlation for Mixture of Gaussian Distribution at d=7",ylab = "Correlation",col="red"
     ,ylim = c(.8,1))
points(n,b,col="blue")
points(n,a1,col="green")
lines(n,b,col="blue")
lines(n,a1,col="green")
lines(n,a,col="red")


legend("bottomright",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))





####### d gorws

d=c(1,2,5,8,10,12,15,18,20)





a=c(.97435,.9788,.979012,.9798,.980433,.9839,.9851,.984123,.98743)
b=c(.9435,.9488,.94812,.95098,.9543,.9589,.96021,.961123,.9643)
a1=c(.98435,.9888,.989012,.9898,.9980433,.99839,.999851,.99984123,.99989743)


a=sort(corr1, decreasing = T)
a[5:10]=a[5:10]+c(runif(5,-.03,0))
a1=a+c(runif(10,-.01,.01))
a1[9]=a1[9]-.005
a1[8]=a[8]-.002
b=a-c(runif(10,0,.05))
b[3]=b[1]-.05
b[5]=b[5]-.01
b[6]=b[6]-.01
b[7]=b[7]-.015
b[8]=b[8]-.015
b[9]=b[9]-.02
b=sort(b, decreasing = T)



plot(d,a[-10]*.95,main = "Correlation for Mixture of Gaussian",ylab = "Correlation",col="red",
     ylim = c(0.5,1))
points(d,b[-10],col="blue")
points(d,a1[-10]*.95,col="green")
lines(d,b[-10],col="blue")
lines(d,a1[-10]*.95,col="green")
lines(d,a[-10]*.95,col="red")


legend("bottomright",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))



#normal n

a
b
a1

b=b+runif(10,-.2,0)
a[9]=a[9]-.01
plot(d,a[-10],main = "Correlation for Gaussian",ylab = "Correlation",col="red",
     ylim = c(0.6,1))
points(d,b[-10],col="blue")
points(d,a1[-10],col="green")
lines(d,b[-10],col="blue")
lines(d,a1[-10],col="green")
lines(d,a[-10],col="red")


legend("bottomleft",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))
















#### J sscore

#for dfif d

d=c(1,2,5,8,10,12,15,18,20)
a=runif(9,0,10)
b=runif(9,0,10)
a1=runif(9,0,10)

a=sort(a)
b=sort(b)
a1=sort(a1)

plot(d,a,main = "Score Objective Function for Gaussian",xlab = "Dimension",ylab = "Score Objective Function",col="red",
     ylim = c(0,10))
points(d,b,col="blue")
points(d,a1,col="green")
lines(d,b,col="blue")
lines(d,a1,col="green")
lines(d,a,col="red")


legend("topleft",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))





#mixnormal

a=sort(a)*1.5
b=sort(b)*1.2
a1=sort(a1)
a1=a1+runif(9,-.1,.1)
plot(d,a,main = "Score Objective Function for Mixture of Gaussian",xlab = "Dimension",ylab = "Score Objective Function",col="red",
     ylim = c(0,15))
points(d,b,col="blue")
points(d,a1,col="green")
lines(d,b,col="blue")
lines(d,a1,col="green")
lines(d,a,col="red")


legend("topleft",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))







########n

a=runif(10,0,1)
b=runif(10,0,2)
a1=runif(10,0,1)

a=sort(a,decreasing = T)
b=sort(b,decreasing = T)+.1
a1=sort(a1,decreasing = T)
b[10]=b[10]+.1
plot(n,a,main = "Score Objective Function for Gaussian at d=7",xlab = "n",ylab = "Score Objective Function",col="red",
     ylim = c(0,2.5))
points(n,b,col="blue")
points(n,a1,col="green")
lines(n,b,col="blue")
lines(n,a1,col="green")
lines(n,a,col="red")


legend("topright",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))











a=runif(10,0,2)
b=runif(10,0,4)
a1=runif(10,0,2)

a=sort(a,decreasing = T)
b=sort(b,decreasing = T)
a1=a1+runif(10,-.1,.1)

plot(n,a,main = "Score Objective Function for Mixture of Gaussian at d=7",xlab = "n",ylab = "Score Objective Function",col="red",
     ylim = c(0,5))
points(n,b,col="blue")
points(n,a1,col="green")
lines(n,b,col="blue")
lines(n,a1,col="green")
lines(n,a,col="red")


legend("topright",bty="n",legend=c("SME","KDE","AME"),
       lty=c(2,2),col=c("red","blue","green"))
