###########exp

n=500
set.seed(1)
xa=xb=rexp(n, 2)
aa=xa
set.seed(1)
anew=rnorm(n, 0, 100)
q=anew

sme_den_exp=function(n,s,l,xa,aa,q)
{
  qw1=Sys.time()
  set.seed(1)
  xa=xb=rexp(n, 2)
  aa=xa
  set.seed(1)
  anew=rnorm(n, 0, 100)
  q=anew

  n=length(xa)
  g=matrix(NA,n,n)
  
  for(a in 1:n){
    for(b in 1:n){
      
      
      g[a,b]=exp(-s*(xa[a]-xb[b])^2)*2*(s)*(1-2*s*(xa[a]-xb[b])^2)
      
    }
  }
  
  h1=matrix(NA,n,n)
  for(a in 1:n){
    for(b in 1:n){
      
      h1[a,b]=exp(-s*(xa[a]-xb[b])^2)*4*(s^2)*(3-2*((xa[a]-xb[b])^2))*(xa[a]-xb[b])
    }
  } 
  
  h11=apply(h1,1,sum)/n
  
  h2=matrix(NA,n,n)
  for(a in 1:n){
    for(b in 1:n){
      
      
      h2[a,b]=exp(-s*(xa[a]-xb[b])^2)*2*(s)*(1-2*s*(xa[a]-xb[b])^2)
      
      
    }
  }
  
  
  
  h22=matrix(NA,n,n)
  for(a in 1:n){
    for(b in 1:n){
      
      h22[a,b]=h2[a,b]*(-xb[b]/100)
      
    }
  }  
  
  h222=apply(h22,1,sum)/n
  
  h=h11+h222
  
  
  
  inverse=solve(g+n*l*diag(1,n,n))
  beta=inverse%*%h/l
  
  
  
  k=length(aa)
  
  x1=matrix(NA,n,k)
  for(a in 1:n){
    for(b in 1:k){
      x1[a,b]=exp(-s*(xa[a]-aa[b])^2)*-2*(s)*(xa[a]-aa[b])
      
    }
  }  
  
  
  
  x11=matrix(NA,n,k)
  for(a in 1:n){
    for(b in 1:k){
      
      x11[a,b]=x1[a,b]%*%(-xa[a]/100)
    }
  }  
  
  
  xi1=apply(x11,2,sum)/n
  
  
  x22=matrix(NA,n,k)
  
  for(a in 1:n){
    for(b in 1:k){
      
      x22[a,b]=exp(-s*(xa[a]-aa[b])^2)*2*(s)*(2*s*(xa[a]-aa[b])^2-1)
      
    }
  }  
  
  
  x2=apply(x22,2,sum)/n
  
  xi2=x2
  
  xi=xi1+xi2
  
  
  bx1=matrix(NA,n,k)
  for(i in 1:n){
    for(j in 1:k){
      bx1[i,j]=beta[i]*x1[i,j]
    }
  }
  
  
  f1=-xi/l+ apply(bx1,2,sum)
  
  
  ##################for q
  k=length(q)
  x1=matrix(NA,n,k)
  for(a in 1:n){
    for(b in 1:k){
      x1[a,b]=exp(-s*(xa[a]-q[b])^2)*-2*(s)*(xa[a]-q[b])
      
    }
  }  
  
  
  
  x11=matrix(NA,n,k)
  for(a in 1:n){
    for(b in 1:k){
      
      x11[a,b]=x1[a,b]%*%(-xa[a]/100)
    }
  }  
  
  
  xi1=apply(x11,2,sum)/n
  
  
  x22=matrix(NA,n,k)
  
  for(a in 1:n){
    for(b in 1:k){
      
      x22[a,b]=exp(-s*(xa[a]-q[b])^2)*2*(s)*(2*s*(xa[a]-q[b])^2-1)
      
    }
  }  
  
  
  x2=apply(x22,2,sum)/n
  
  xi2=x2
  
  xi=xi1+xi2
  
  
  bx1=matrix(NA,n,k)
  for(i in 1:n){
    for(j in 1:k){
      bx1[i,j]=beta[i]*x1[i,j]
    }
  }
  
  
  f2=-xi/l+ apply(bx1,2,sum)
  
  deno=mean(exp(f2))
  
  p=exp(f1)*(dnorm(aa,0,100))/deno
  plot(aa,dexp(aa,2))
  points(aa,p,col="pink")
  qw2=Sys.time()
  qw3=qw2-qw1
  cv=list(sum(((p-(dexp(aa,2)))^2)/length(p)),p,qw3)
  return(cv)
  
}
fit_exp1=sme_den_exp(n=100,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp2=sme_den_exp(n=200,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp3=sme_den_exp(n=300,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp4=sme_den_exp(n=400,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp5=sme_den_exp(n=500,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp6=sme_den_exp(n=600,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp7=sme_den_exp(n=700,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp8=sme_den_exp(n=800,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp9=sme_den_exp(n=900,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp10=sme_den_exp(n=1000,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 
fit_exp11=sme_den_exp(n=2000,s=1,l=1/1000,xa=xa,aa=xa,q=anew) 


n1=c(100,200,300,400,500,600,700,800,900,1000)
tt=c(fit_exp1[[3]],fit_exp2[[3]],fit_exp3[[3]],fit_exp4[[3]],fit_exp5[[3]],fit_exp6[[3]],
     fit_exp7[[3]],fit_exp8[[3]],fit_exp9[[3]],fit_exp10[[3]],fit_exp11[[3]],fit_exp12[[3]])
plot(n,tt,ylab = "time in seconds",xlab = "n",main = "plot for Time")


kde=function(n,h)
{
  qw1=Sys.time()

  set.seed(1)
  x=rexp(n,2)
  
  
  h1=1.06*((var(x))^.5)*n^(-1/5)
  
  a=x
  
  
  
  f=matrix(NA,n,1)
  for(i in 1:n){
    
    f[i]=sum(exp(-0.5*((a[i]-x)/h)^2))/(n*h*sqrt(2*pi))
    
  }
  
  f1=matrix(NA,n,1)
  for(i in 1:n){
    
    f1[i]=sum(exp(-0.5*((a[i]-x)/h1)^2))/(n*h1*sqrt(2*pi))
    
  }
  
  
  
  
  plot(a,dexp(a,2),main="kde")
  points(a,f,col="blue")
  points(a,f1,col="red")
  cv1=sum(((f-(dexp(a,2)))^2)/length(f))
  cv2=sum(((f1-(dexp(a,2)))^2)/length(f))
  qw2=Sys.time()
  qw3=qw2-qw1
  cv=list(cv1,cv2,qw3,f,f1)
  return(cv)
  
}


fit_exp1=kde(n=100,h=.035)
fit_exp2=kde(n=200,h=.035)
fit_exp3=kde(n=300,h=.035)
fit_exp4=kde(n=400,h=.035)
fit_exp5=kde(n=500,h=.035)
fit_exp6=kde(n=600,h=.035)
fit_exp7=kde(n=700,h=.035)
fit_exp8=kde(n=800,h=.035)
fit_exp9=kde(n=900,h=.035)
fit_exp10=kde(n=1000,h=.035)
fit_exp11=kde(n=2000,h=.035)
fit_exp12=kde(n=5000,h=.035)

tt1=c(fit_exp1[[3]],fit_exp2[[3]],fit_exp3[[3]],fit_exp4[[3]],fit_exp5[[3]],fit_exp6[[3]],
     fit_exp7[[3]],fit_exp8[[3]],fit_exp9[[3]],fit_exp10[[3]],fit_exp11[[3]],fit_exp12[[3]])

aexp=function(n,s,l)
{
  qw1=Sys.time()
  
  set.seed(1)
  xa=rexp(n,2)
  aa=xa
  k1=matrix(NA,n,n)
  for (i in 1:n){
    for (j in 1:n){
      k1[i,j]=-2*s*(xa[i]-xa[j])*exp(-s*(xa[i]-xa[j])^2)
    }
  }
  
  a=matrix(NA,n,n)
  for (i in 1:n){
    for (j in 1:n){
      a[i,j]=sum(k1[,i]*k1[,j])/n
    }
  }
  
  k2=matrix(NA,n,n)
  for (i in 1:n){
    for (j in 1:n){
      k2[i,j]=2*s*(2*s*(xa[i]-xa[j])^2-1)*exp(-s*(xa[i]-xa[j])^2)
    }
  }
  
  
  q1=-xa/100
  
  b=matrix(NA,n,1)
  for (i in 1:n){
    b[i]=sum(k2[,i]+q1[i]*k1[,i])/n
  }
  
  
  al=-solve(a+l*diag(1,n,n))%*%b
  
  sa=length(aa)
  
  kx=matrix(NA,sa,n)
  for (i in 1:sa){
    for (j in 1:n){
      kx[i,j]=exp(-s*(aa[i]-xa[j])^2)
    }
  }
  
  
  
  ff=matrix(NA,sa,n)
  for (i in 1:sa){
    for (j in 1:n){
      ff[i,j]=al[j]*kx[i,j]
    }
  }
  
  f1=apply(ff,1,sum)
  
  
  set.seed(10)
  daa=rnorm(1000,0,100)
  ee=length(daa)
  dkx=matrix(NA,ee,n)
  for (i in 1:ee){
    for (j in 1:n){
      dkx[i,j]=exp(-s*(daa[i]-xa[j])^2)
    }
  }
  
  
  
  dff=matrix(NA,ee,n)
  for (i in 1:ee){
    for (j in 1:n){
      dff[i,j]=al[j]*dkx[i,j]
    }
  }
  
  df1=apply(dff,1,sum)
  
  
  
  deno=mean(exp(df1))
  #infinity
  
  p=exp(f1)*dnorm(aa,0,100)/deno
  
  plot(aa,dexp(aa,2),xlim=c(0,3))
  points(aa,p,col="green")
  qw2=Sys.time()
  qw3=qw2-qw1
  cv3=list(sum(((p-(dexp(aa,2)))^2)/length(p)),p,qw3)
  return(cv3)
}

fit_exp1=aexp(n=100,s=30,l=6)
fit_exp2=aexp(n=200,s=30,l=6)
fit_exp3=aexp(n=300,s=30,l=6)
fit_exp4=aexp(n=400,s=30,l=6)
fit_exp5=aexp(n=500,s=30,l=6)
fit_exp6=aexp(n=600,s=30,l=6)
fit_exp7=aexp(n=700,s=30,l=6)
fit_exp8=aexp(n=800,s=30,l=6)
fit_exp9=aexp(n=900,s=30,l=6)
fit_exp10=aexp(n=1000,s=30,l=6)
fit_exp11=aexp(n=2000,s=30,l=6)
fit_exp12=aexp(n=5000,s=30,l=6)


tt2=c(fit_exp1[[3]],fit_exp2[[3]],fit_exp3[[3]],fit_exp4[[3]],fit_exp5[[3]],fit_exp6[[3]],
      fit_exp7[[3]],fit_exp8[[3]],fit_exp9[[3]],fit_exp10[[3]],fit_exp11[[3]])


plot(n1[1:11],tt[1:11],ylab = "time in seconds",xlab = "n",main = "Computational Cost for each method",col="red",ylim=c(0,18))
lines(n1[1:11],tt[1:11],col="red")
lines(n1[1:11],tt1[1:11],col="blue")
lines(n1[1:11],tt2,col="green")
points(n1[1:11],tt1[1:11],col="blue")
points(n1[1:11],tt2,col="green")


legend("topleft",bty="n",legend=c("KDE","SME","AME"),
       lty=c(2,2),col=c("blue","red","green"))

fit_exp[[1]]
try1[[1]]
try1[[2]]
try[[1]]
