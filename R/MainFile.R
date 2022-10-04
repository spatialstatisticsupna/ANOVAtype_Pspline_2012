


########### THIS FILE CONTAING THE CODE TO FIT A P-spline ANOVA type model.
#### Ugarte, M. D., Goicoa, T., Etxeberria, J., & Militino, A. F. (2012). 
#### A P-spline ANOVA type model in space-time disease mapping. 
#### Stochastic Environmental Research and Risk Assessment, 26(6), 835-845.


#### Author: Jaione Etxeberria (jaione.etxeberria@unavarra.es)

rm(list=ls(all=TRUE))

library(splines)
library(mgcv)
library(nlme)
library(MASS)
library(mvtnorm)
library(maptools)
library(RColorBrewer)
library(splancs)
library(spdep)

setwd("C:/......./ANOVA") #Directorio donde se ubica la carpeta ANOVA


dat<-read.table(file="Data/ProstH.txt", header = TRUE)

################################################################################
################## CALCULO DE LOS CASOS ESPERADOS ##############################



bspline<-function(x, xl, xr, ndx, bdeg){
 dx <- (xr-xl)/ndx
 knots <- seq(xl-bdeg*dx, xr+bdeg*dx, by=dx)
 B <- spline.des(knots,x,bdeg+1,0*x)$design
 B
 }

bdeg<-3
pord<-2



x1<-unique(dat$long)
#x1<-((x1-mean(x1))/sqrt(var(x1)))
x1<-(x1-min(x1))/1000
#x1<-(x1-min(x1))/max(x1)


x2<-unique(dat$lat)
#x2<-((x2-mean(x2))/sqrt(var(x2)))
x2<-(x2-min(x2))/1000
#x2<-(x2-min(x2))/max(x2)


x3<-unique(dat$agno)
#x3<-((x3-mean(x3))/sqrt(var(x3))) 
x3<-(x3-min(x3))/1000    
#x3<-1:10
#x3<-(x3-mean(x3))/max(x3)


ndx1<-10 #9#min(ceiling(length(x1)/4),40)#8
dis1 <- (max(x1)-min(x1))/ndx1
x1l<-min(x1)-dis1*0.05
x1r<-max(x1)+dis1*0.05
dx1 <- (x1r-x1l)/ndx1
knots1<-seq(x1l-bdeg*dx1, x1r+bdeg*dx1, by=dx1)


ndx2<-10 #9#min(ceiling(length(x2)/4),40) #x
dis2 <- (max(x2)-min(x2))/ndx2
x2l<-min(x2)-dis2*0.05
x2r<-max(x2)+dis2*0.05
dx2 <- (x2r-x2l)/ndx2
knots2<-seq(x2l-bdeg*dx2, x2r+bdeg*dx2, by=dx2)


ndx3<-8 #min(ceiling(length(x3)/4),40)#2#
dis3 <- (max(x3)-min(x3))/ndx3
x3l<-min(x3)-dis3*0.05
x3r<-max(x3)+dis3*0.05
dx3 <- (x3r-x3l)/ndx3
knots3<-seq(x3l-bdeg*dx3, x3r+bdeg*dx3, by=dx3)



B1<-spline.des(knots1,x1,bdeg+1,0*x1)$design
B2<-spline.des(knots2,x2,bdeg+1,0*x2)$design
B3<-spline.des(knots3,x3,bdeg+1,0*x3)$design


####### Producto tensorial por filas

Rten<-function(X1,X2){
one1<-matrix(1,1,ncol(X1))
one2<-matrix(1,1,ncol(X2))
kronecker(X1,one2)*kronecker(one1,X2)
}

########  H-transform of an array A by a matrix X (Currie et al 2006)

H = function(X, A){ d = dim(A)
M = matrix(A, nrow = d[1])
XM=X %*% M
array(XM, c(nrow(XM), d[-1]))
}
########## Rotation of an array A (Currie et al 2006)

Rotate = function(A){ d = 1:length(dim(A))
d1 = c(d[-1], d[1])
aperm(A, d1)
}

######### Rotated H-transform of an array A by a matrix X (Currie et al 2006)
RH = function(X, A){ Rotate(H(X, A))}



###### Matrices penalizacion

###### Penalizacion 1

m1=ncol(B1)
D1=diff(diag(m1),differences=pord)
P1.svd=svd(t(D1)%*%D1)
U1s=(P1.svd$u)[,1:(m1-pord)] # le quito dos columnas por los dos autovalores que son 0
U1n=(P1.svd$u)[,(m1-pord+1):m1] # dos columnas de los dos autovalores que son 0
d1=(P1.svd$d)[1:(m1-pord)]
Delta1=diag(d1)


###### Penalizacion 2

m2=ncol(B2)
D2=diff(diag(m2),differences=pord)
P2.svd=svd(t(D2)%*%D2)
U2s=(P2.svd$u)[,1:(m2-pord)] # le quito dos columnas por los dos autovalores que son 0
U2n=(P2.svd$u)[,(m2-pord+1):m2] # dos columnas de los dos autovalores que son 0
d2=(P2.svd$d)[1:(m2-pord)]
Delta2=diag(d2)


###### Penalizacion 3

m3=ncol(B3)
D3=diff(diag(m3),differences=pord)
P3.svd=svd(t(D3)%*%D3)
U3s=(P3.svd$u)[,1:(m3-pord)] # le quito dos columnas por los dos autovalores que son 0
U3n=(P3.svd$u)[,(m3-pord+1):m3] # dos columnas de los dos autovalores que son 0
d3=(P3.svd$d)[1:(m3-pord)]
Delta3=diag(d3)


############ MATRIZ X


X1<-cbind(1,x1)
X2<-cbind(1,x2)
Xt<-cbind(1,x3)

Xs<-Rten(X2,X1)

One.t<-Xt[,1]
One.s<-Xs[,1]

xs<-Xs[,4]
xt<-x3
xhat<-cbind(x1,x2,xs) 

#dim(kronecker(One.s,xt))
#dim(kronecker(One.t,xt))

X<-cbind(kronecker(Xs,One.t),kronecker(One.s,xt),kronecker(xhat,xt))

############ MATRIZ Z


Z1<-B1%*%U1s
Z2<-B2%*%U2s
Z3<-B3%*%U3s

Z.s<-cbind(Rten(Z2,X1),Rten(X2,Z1),Rten(Z2,Z1))
Z.t<-Z3

One.t<-Xt[,1]
One.s<-Xs[,1]


Z.b1<-kronecker(Z.s,One.t)
Z.b2<-kronecker(One.s,Z.t)  
Z.b3<-kronecker(Z.s,xt)
Z.b4<-kronecker(xhat,Z.t)
Z.b5<-kronecker(Z.s,Z.t)

Z<-cbind(Z.b1,Z.b2,Z.b3,Z.b4,Z.b5)

k1<-m1
k2<-m2
k3<-m3

dump("k1",file="dumpsdata/k1")
dump("k2",file="dumpsdata/k2")
dump("k3",file="dumpsdata/k3")


dump("X",file="dumpsdata/X")
dump("Z",file="dumpsdata/Z")

dump("Delta1",file="dumpsdata/Delta1")
dump("Delta2",file="dumpsdata/Delta2")
dump("Delta3",file="dumpsdata/Delta3")

dump("B1",file="dumpsdata/B1")
dump("B2",file="dumpsdata/B2")
dump("B3",file="dumpsdata/B3")



