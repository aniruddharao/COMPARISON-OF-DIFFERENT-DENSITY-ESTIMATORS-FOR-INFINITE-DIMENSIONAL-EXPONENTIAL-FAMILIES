###########exp

n=500
set.seed(1)
xa=xb=rexp(n, 2)
aa=xa
set.seed(1)
anew=rnorm(n, 0, 100)
q=anew

sme_den_exp=function(s,l,xa,aa,q)
{
  
  
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
  
  cv=list(sum(((p-(dexp(aa,2)))^2)/length(p)),p)
  return(cv)
  
}
fit_exp=sme_den_exp(s=1,l=1/1000,xa=xa,aa=xa,q=anew) 


kde=function(n,h)
{
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
  cv=list(cv1,cv2,f,f1)
  return(cv)
  
}


try1=kde(n=500,h=.035)


aexp=function(n,s,l)
{
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
  
  cv3=list(sum(((p-(dexp(aa,2)))^2)/length(p)),p)
  return(cv3)
}

try=aexp(n=500,s=30,l=6)


plot(aa,dexp(aa,2),ylab = "p(x)",xlab = "x",main = "Correlation for Gaussian at n=500")
points(aa,fit_exp[[2]],col="red")
points(aa,try1[[3]],col="blue")
points(aa,try[[2]],col="green")



legend("topright",bty="n",legend=c("Original","SME","KDE","AME"),
       lty=c(2,2),col=c("black","red","blue","green"))

fit_exp[[1]]
try1[[1]]
try1[[2]]
try[[1]]


###########nor



n=500
set.seed(1)
xa=xb=rnorm(n,0,1)
aa=xa
set.seed(1)
anew=rnorm(n, 0, 100)
q=anew

den_normal=function(s,l,xa,aa,q)
{
  
  
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
  plot(aa,dnorm(aa,0,1))
  points(aa,p,col="pink")
  
  cv=list(sum(((p-(dnorm(aa,0,1)))^2)/length(p)),p)
  return(cv)
  
}
fit_normal=den_normal(s=1/50,l=1/220,xa=xa,aa=xa,q=anew) 

kde=function(n,h)
{
  set.seed(1)
  x=rnorm(n,0,1)
  
  
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
  
  
  
  
  plot(a,dnorm(a,0,1),ylim=c(0,.5),main="KDE")
  points(a,f,col="blue")
  points(a,f1,col="red")
  cv1=sum(((f-(dnorm(a,0,1)))^2)/length(f))
  cv2=sum(((f1-(dnorm(a,0,1)))^2)/length(f))
  cv=list(cv1,cv2,f,f1)
  return(cv)
  
}


try1=kde(n=500,h=.35)




anorm=function(n,s,l)
{
  set.seed(1)
  xa=rnorm(n,0,1)
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
  
  plot(aa,dnorm(aa,0,1),main = "approx")
  points(aa,p,col="green")
  
  cv=list(sum(((p-(dnorm(aa,0,1)))^2)/length(p)),p)
  return(cv)
}


try2=anorm(n=500,s=1/4,l=1/3.2)



plot(aa,dnorm(aa,0,1),ylab = "p(x)",xlab = "x",main = "Normal")
points(aa,fit_normal[[2]],col="red")
points(aa,try1[[3]],col="blue")
points(aa,try2[[2]],col="green")



legend("topright",bty="n",legend=c("Original","SME","KDE","AME"),
       lty=c(2,2),col=c("black","red","blue","green"))

fit_exp[[1]]
try1[[1]]
try1[[2]]
try2[[1]]




#####cacuhy


kde=function(n,h)
{
  set.seed(1)
  x=rcauchy(n,0,2)
  
  
  h1=1.06*((var(x))^.5)*n^(-1/5)
  a=x
  aa=x
  
  
  
  f=matrix(NA,n,1)
  for(i in 1:n){
    
    f[i]=sum(exp(-0.5*((a[i]-x)/h)^2))/(n*h*sqrt(2*pi))
    
  }
  
  f1=matrix(NA,n,1)
  for(i in 1:n){
    
    f1[i]=sum(exp(-0.5*((a[i]-x)/h1)^2))/(n*h1*sqrt(2*pi))
    
  }
  
  
  
  
  plot(a,dcauchy(a,0,2),ylim=c(0,.25),xlim=c(-10,10),main="kde")
  points(aa,f,col="blue")
  points(a,f1,col="red")
  cv1=sum(((f-(dcauchy(a,0,2)))^2)/length(f))
  cv2=sum(((f1-(dcauchy(a,0,2)))^2)/length(f))
  cv=list(cv1,cv2,f,f1)
  return(cv)
  
}


try1=kde(n=500,h=.46)





anorm=function(n,s,l)
{
  set.seed(1)
  xa=rcauchy(n,0,2)
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
  
  plot(aa,dcauchy(aa,0,2),xlim=c(-10,10),main = "approx")
  points(aa,p,col="green")
  
  cv=list(sum(((p-(dcauchy(aa,0,2)))^2)/length(p)),p)
  
}


try=anorm(n=500,s=1/10,l=1/25)




den_cauchy=function(s,l,xa,aa,q)
{
  
  
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
  plot(aa,dcauchy(aa,0,2),xlim=c(-10,10))
  points(aa,p,col="pink")
  
  cv=list(sum(((p-(dcauchy(aa,0,2)))^2)/length(p)),p)
  return(cv)
  
}

fit_cauchy3=den_cauchy(l=1/30,s=1/52,xa=xa,aa=xa,q=anew) 



plot(aa,dcauchy(aa,0,2),ylab = "p(x)",xlab = "x",main = "Cauchy",xlim = c(-10,10))
points(aa,fit_cauchy3[[2]],col="red")
points(aa,try1[[3]],col="blue")
points(aa,try[[2]],col="green")



legend("topright",bty="n",legend=c("Original","SME","KDE","AME"),
       lty=c(2,2),col=c("black","red","blue","green"))

fit_exp[[1]]
try1[[1]]
try1[[2]]
try2[[1]]






######mix normal


kde=function(n,h)
{
  set.seed(1)
  
  x1=rnorm(n/2,0,1)
  set.seed(1)
  x2=rnorm(n/2,4,2)
  
  x=c(x1,x2)
  
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
  
  
  
  
  plot(a,.5*dnorm(a,0,1) +.5*dnorm(a,4,2),main="kde")
  points(a,f,col="blue")
  points(a,f1,col="red")
  cv1=sum(((f-(.5*dnorm(a,0,1) +.5*dnorm(a,4,2)))^2)/length(f))
  cv2=sum(((f1-(.5*dnorm(a,0,1) +.5*dnorm(a,4,2)))^2)/length(f))
  cv=list(cv1,cv2,f,f1)
  return(cv)
  
}


try4=kde(n=500,h=.46)


anorm=function(n,s,l)
{
  set.seed(1)
  x1=rnorm(n/2,0,1)
  set.seed(1)
  x2=rnorm(n/2,4,2)
  
  xa=c(x1,x2)
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
  
  p=exp(f1)*(dnorm(aa,0,100))/deno
  
  plot(aa,(.5*dnorm(aa,0,1) +.5*dnorm(aa,4,2)),ylim=c(0,.25))
  points(aa,p,col="green")
  
  cv=list(sum(((p-(.5*dnorm(aa,0,1) +.5*dnorm(aa,4,2)))^2)/length(p)),p)
  return(cv)
}


try=anorm(n=500,s=1/4,l=1/10)


den_mix=function(s,l,xa,aa,q)
{
  
  
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
  
  
  f2=-xi/l + apply(bx1,2,sum)
  
  deno=mean(exp(f2))
  
  p=exp(f1)*(dnorm(aa,0,100))/deno
  plot(aa,.5*dnorm(aa,0,1)+.5*dnorm(aa,4,2))
  points(aa,p,col="pink")
  
  cv=list(sum(((p-(.5*dnorm(aa,0,1)+.5*dnorm(aa,4,2)))^2))/length(p),p)
  return(cv)
  
}

fit_mix=den_mix(s=1,l=1/100000,xa=xa,aa=xa,q=anew) 

plot(aa,(.5*dnorm(aa,0,1) +.5*dnorm(aa,4,2)),ylab = "PDF",xlab = "Range",main = "Mix Normal")
points(aa,fit_mix[[2]],col="red")
points(aa,try4[[3]],col="blue")
points(aa,try[[2]],col="green")

legend("topright",bty="n",legend=c("Original","SME","KDE","AME"),
       lty=c(2,2),col=c("black","red","blue","green"))

fit_mix[[1]]
try4[[1]]
try4[[2]]
try[[1]]
