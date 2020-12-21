library(plot3D)
library(plotrix)
library(plotly)
library(rgl)
library(scatterplot3d)


x1=runif(n=300, min = -5, max = 5)
x2=runif(n=300, min = -5, max = 5)
y=(((x1)^2)/4+((x2)^2)/4)*2

scatter3D(x1, x2, y, colvar = y, col = NULL, add = FALSE)
points3D(x1, x2, y)

lines3D(x1, x2, y)

x=data.frame(x1,x2)
a=data.frame(runif(n=50, min = -5, max = 5),runif(n=50, min = -5, max = 5))
b=(((a[,1])^2)/4+((a[,2])^2)/4)*2
points3D(a[,1],a[,2], b)


x11 = seq(-5,5,len=1000)
x22 = seq(-5,5,len=1000)


###checking

l=2^(-14)
s=2^(-4)
n=length(x[,1])
k=matrix(NA,n,n)
for(i in 1:n){
  for(j in 1:n){
    k[i,j]=exp(-s*(((x[i,1]-x[j,1])^2)+((x[i,2]-x[j,2])^2)))
  }
}  

#m=10 # no. of vectors
#for(i in 1:n){
#  for(j in 1:n){
#    sum=0
#    for(h in 1:m){
#      sum= sum+((x[i,k]-x[j,k])^2)
#    }
#    k[i,j]=exp(-s*sum)
#  }
#}  

I=diag(n)


lx3=length(x[,1])
x3=x
k1=matrix(NA,n,lx3)
for(i in 1:n){
  for(j in 1:lx3){
    k1[i,j]=exp(-s*(((x[i,1]-x3[j,1])^2)+((x[i,2]-x3[j,2])^2)))
  }
}  


inverse=solve(k+l*I)

ff=y%*%inverse%*%k1
cv=sum(((ff-y)^2)/length(ff))
cv
points3D(x1, x2, ff,main='estimated')
points3D(x1, x2, y,main='original')

fun2(x =f1, y =ys1, a=a, b=b, s=2^(-4), l = 2^(-14))
fun2(x =x, y =y, a=x, b=y, s=2^(-4), l = 2^(-14))
points3D(f1[,1],f1[,2],ys1)


######function
fun2=function(x,y,a,b, s,l)
{
  
  
  
  n=length(x[,1])
  k=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      k[i,j]=exp(-s*(((x[i,1]-x[j,1])^2)+((x[i,2]-x[j,2])^2)))
    }
  }  
  
 
  I=diag(n)
  
  
  lx3=length(a[,1])
  x3=a
  k1=matrix(NA,n,lx3)
  for(i in 1:n){
    for(j in 1:lx3){
      k1[i,j]=exp(-s*(((x[i,1]-x3[j,1])^2)+((x[i,2]-x3[j,2])^2)))
    }
  }  
  
  
  inverse=solve(k+l*I)
  
  ff=y%*%inverse%*%k1
  cv=sum(((ff-b)^2)/length(ff))
  return(cv)
}
l=c(2^(-14),2^(-13),2^(-12),2^(-11),2^(-10),2^(-9))
s=c(2^(-6),2^(-5),2^(-4),2^(-3),2^(-2),1/2)
v=length(l)
res2 = matrix(NA, v,v)


for(i in 1:v){
  for(j in 1:v){
    
    res2[i,j] = fun2(x = x , y =y, a=a, b=b, s=s[i], l = l[j])
  }
}
res2



########5 cv
d=data.frame(x,y)

xq1=matrix(d$x1,200,5)
xq2=matrix(d$x2,200,5)
ys=matrix(d$y,200,5)



xs1=xq1[,-1]
xs1=as.vector((xs1))
xs2=xq1[,-2]
xs2=as.vector(xs2)
xs3=xq1[,-3]
xs3=as.vector((xs3))
xs4=xq1[,-4]
xs4=as.vector((xs4))
xs5=xq1[,-5]
xs5=as.vector((xs5))


xss1=xq2[,-1]
xss1=as.vector((xss1))
xss2=xq2[,-2]
xss2=as.vector(xss2)
xss3=xq2[,-3]
xss3=as.vector((xss3))
xss4=xq2[,-4]
xss4=as.vector((xss4))
xss5=xq2[,-5]
xss5=as.vector((xss5))

ax1=xq1[,1]
ax1=as.vector((ax1))
ax2=xq1[,2]
ax2=as.vector((ax2))
ax3=xq1[,3]
ax3=as.vector((ax3))
ax4=xq1[,4]
ax4=as.vector((ax4))
ax5=xq1[,5]
ax5=as.vector((ax5))

axx1=xq2[,1]
axx1=as.vector((axx1))
axx2=xq2[,2]
axx2=as.vector((axx2))
axx3=xq2[,3]
axx3=as.vector((axx3))
axx4=xq2[,4]
axx4=as.vector((axx4))
axx5=xq2[,5]
axx5=as.vector((axx5))



f1=data.frame( xs1,xss1)
f2=data.frame( xs2,xss2)
f3=data.frame( xs3,xss3)
f4=data.frame( xs4,xss4)
f5=data.frame( xs5,xss5)

e1=data.frame(ax1,axx1)
e2=data.frame(ax2,axx2)
e3=data.frame(ax3,axx3)
e4=data.frame(ax4,axx4)
e5=data.frame(ax5,axx5)


ys1=(((f1[,1])^2)/4+((f1[,2])^2)/4)*2
ys1=as.vector((ys1))
ys2=(((f2[,1])^2)/4+((f2[,2])^2)/4)*2
ys2=as.vector((ys2))
ys3=(((f3[,1])^2)/4+((f3[,2])^2)/4)*2
ys3=as.vector((ys3))
ys4=(((f4[,1])^2)/4+((f4[,2])^2)/4)*2
ys4=as.vector((ys4))
ys5=(((f5[,1])^2)/4+((f5[,2])^2)/4)*2
ys5=as.vector((ys5))

ay1=(((e1[,1])^2)/4+((e1[,2])^2)/4)*2
ay1=as.vector((ay1))
ay2=(((e2[,1])^2)/4+((e2[,2])^2)/4)*2
ay2=as.vector((ay2))
ay3=(((e3[,1])^2)/4+((e3[,2])^2)/4)*2
ay3=as.vector((ay3))
ay4=(((e4[,1])^2)/4+((e4[,2])^2)/4)*2
ay4=as.vector((ay4))
ay5=(((e5[,1])^2)/4+((e5[,2])^2)/4)*2
ay5=as.vector((ay5))




l=c(2^(-14),2^(-13),2^(-12),2^(-11),2^(-10),2^(-9))
s=c(2^(-6),2^(-5),2^(-4),2^(-3),2^(-2),1/2)

v=length(s)

points3D(f1[,1],f1[,2],ys1)
points3D(f2[,1],f2[,2],ys2)
points3D(f3[,1],f3[,2],ys3)
points3D(f4[,1],f4[,2],ys4)
points3D(f5[,1],f5[,2],ys5)

res1 = matrix(NA, v,v)
for(i in 1:v){
  for(j in 1:v){
    
    res1[i,j] = fun2(x =f1 , y =ys1, a=e1, b=ay1, s=s[i], l = l[j])
  }
}
res1


res2 = matrix(NA, v,v)


for(i in 1:v){
  for(j in 1:v){
    
    res2[i,j] = fun2(x = f2 , y =ys2, a=e2, b=ay2, s=s[i], l = l[j])
  }
}
res2

res3 = matrix(NA, v,v)


for(i in 1:v){
  for(j in 1:v){
    
    res3[i,j] = fun2(x =f3 , y =ys3, a=e3, b=ay3, s=s[i], l = l[j])
  }
}
res3

res4 = matrix(NA, v,v)


for(i in 1:v){
  for(j in 1:v){
    
    res4[i,j] = fun2(x = f4, y =ys4, a=e4, b=ay4, s=s[i], l = l[j])
  }
}
res4

res5 = matrix(NA, v,v)


for(i in 1:v){
  for(j in 1:v){
    
    res5[i,j] = fun2(x =f5 , y =ys5, a=e5, b=ay5, s=s[i], l = l[j])
  }
}
res5

######combine

rv1=as.vector(res1)
rv2=as.vector(res2)
rv3=as.vector(res3)
rv4=as.vector(res4)
rv5=as.vector(res5)

rrv1=rank(rv1)
rrv2=rank(rv2)
rrv3=rank(rv3)
rrv4=rank(rv4)
rrv5=rank(rv5)

lf=rep(l,each=v)
sf=rep(s,times=v)

bd=data.frame(sf,lf,rv1,rrv1,rv2,rrv2,rv3,rrv3,rv4,rrv4,rv5,rrv5)

meancv=bd$rv1*bd$rrv1+bd$rv2*bd$rrv2+bd$rv3*bd$rrv3+bd$rv4*bd$rrv4+bd$rv5*bd$rrv5
r=rank(meancv)
bdf=data.frame(sf,lf,meancv,r)
best=which(bdf$r==1)
best1=bdf[best,]
best1


##plotting

l=2^(-14)
s=2^(-4)


n=length(x[,1])
k=matrix(NA,n,n)
for(i in 1:n){
  for(j in 1:n){
    k[i,j]=exp(-s*(((x[i,1]-x[j,1])^2)+((x[i,2]-x[j,2])^2)))
  }
}  


I=diag(n)


lx3=length(x[,1])
x3=x
k1=matrix(NA,n,lx3)
for(i in 1:n){
  for(j in 1:lx3){
    k1[i,j]=exp(-s*(((x[i,1]-x3[j,1])^2)+((x[i,2]-x3[j,2])^2)))
  }
}  


inverse=solve(k+l*I)

ff=y%*%inverse%*%k1
cv=sum(((ff-y)^2)/length(ff))
cv
par(mfrow=(1:2))
points3D(x1, x2, ff)
title("Esitmate")
points3D(x1,x2,y,pch=33)
title("Original")
scatterplot3d(x1,x2,ff)
title("Esitmate")
scatterplot3d(x1,x2,y,color  ="springgreen")
title("Original")


#cool=matrix(NA,n,n)
#for(i in 1:n){
#  for(j in 1:n){
#    cool[i,j]=(((x11[i])^2)/4+((x22[j])^2)/4)*2
#  }
#}  
#cool0=function(x,y)
#{
#  (((x)^2)/4+((y)^2)/4)*2
#}
#cool1=outer(x11,x22,cool0)




########multivariate x
l=2^(-14)
s=2^(-4)
n=300
start_r <- Sys.time()
k=matrix(NA,n,n)
for(i in 1:n){
  for(j in 1:n){
     sum=0
      for(h in 1:2){
        sum= sum+((x[i,h]-x[j,h])^2)
      }
      k[i,j]=exp(-s*sum)
    }
  }  
  
end_r <- Sys.time()
time_r=end_r-start_r


start_t <- Sys.time()
n=length(x[,1])
m=length(x)
k1=matrix(NA,n,n)
for(i in 1:n){
  for(j in 1:n){
   
    k1[i,j]=exp(-s*sum((x[i,]-x[j,])^2))
  }
}    

end_t <- Sys.time()
time_t=end_t-start_t


I=diag(n)


lx3=length(x[,1])
x3=x
k1=matrix(NA,n,lx3)
for(i in 1:n){
  for(j in 1:lx3){
    k1[i,j]=exp(-s*(((x[i,1]-x3[j,1])^2)+((x[i,2]-x3[j,2])^2)))
  }
}  


inverse=solve(k+l*I)

ff=y%*%inverse%*%k1
cv=sum(((ff-y)^2)/length(ff))
cv

