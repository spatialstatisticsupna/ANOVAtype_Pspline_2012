
################################################################################
########### La matriz de covarianza del Spline



I2<-as.matrix(diag(1,nrow=2, ncol=2))
I3<-as.matrix(diag(1,nrow=3, ncol=3))

I12<- diag(1,nrow=k1-2, ncol=k1-2)
I22<- diag(1,nrow=k2-2, ncol=k2-2)
I32<- diag(1,nrow=k3-2, ncol=k3-2)




###### Fs 

M1s<-diag(kronecker(Delta2,I2))
F1<-lambda2*M1s

##1
M2s<-diag(kronecker(I2,Delta1))
F2<-lambda1*M2s

##2 

M3s1<-diag(kronecker(I22,Delta1))
M3s2<-diag(kronecker(Delta2,I12))
##3
F3<-lambda1*M3s1 + lambda2*M3s2

###### Ft
##4
F4<-lambda3*diag(Delta3)

###### Fst
##5

M1st<-diag(kronecker(Delta2,I2))
F5<-tau2*M1st

##6

M2st<-diag(kronecker(I2,Delta1))

F6<-tau1*M2st

##7
M3st1<-diag(kronecker(I22,Delta1)) 
M3st2<-diag(kronecker(Delta2,I12))

F7<-tau1*M3st1 + tau2*M3st2 

##8
M4st<-diag(kronecker(I3,Delta3))

F8<-tau3*M4st

##9
M5st1<-diag(kronecker(kronecker(Delta2,I2),I32))
M5st2<-diag(kronecker(kronecker(I22,I2),Delta3))

F9<-tau2*M5st1+tau3*M5st2

##10
M6st1<-diag(kronecker(kronecker(I2,Delta1),I32))
M6st2<-diag(kronecker(kronecker(I2,I12),Delta3))

F10<-tau1*M6st1+tau3*M6st2

##11
M7st1<-diag(kronecker(kronecker(I22,Delta1),I32))
M7st2<-diag(kronecker(kronecker(Delta2,I12),I32))
M7st3<-diag(kronecker(kronecker(I22,I12),Delta3))

F11<-tau1*M7st1+ tau2*M7st2 + tau3*M7st3


diagonalG<-c(1/F1,1/F2,1/F3,1/F4,1/F5,1/F6,1/F7,1/F8,1/F9,1/F10,1/F11)


dim1<-length(diagonalG)

G<-diag(diagonalG,ncol=dim1,nrow=dim1)

####### DERIVADA  DE  G.b

cero1<-c(rep(0,length(F1))) 
cero2<-c(rep(0,length(F2)))
cero3<-c(rep(0,length(F3)))
cero4<-c(rep(0,length(F4)))
cero5<-c(rep(0,length(F5)))
cero6<-c(rep(0,length(F6)))
cero7<-c(rep(0,length(F7)))
cero8<-c(rep(0,length(F8)))
cero9<-c(rep(0,length(F9)))
cero10<-c(rep(0,length(F10)))
cero11<-c(rep(0,length(F11)))


#### Lambda1
dF1<-diag(c(cero1,M2s,M3s1,cero4,cero5,cero6,cero7,cero8,cero9,cero10,cero11))

#### Lambda2
dF2<-diag(c(M1s,cero2,M3s2,cero4,cero5,cero6,cero7,cero8,cero9,cero10,cero11))

#### Lambda3
dF3<-diag(c(cero1,cero2,cero3,diag(Delta3),cero5,cero6,cero7,cero8,cero9,cero10,cero11))

#### Tau1
dF4<-diag(c(cero1,cero2,cero3,cero4,cero5,M2st,M3st1,cero8,cero9,M6st1,M7st1))

#### Tau2
dF5<-diag(c(cero1,cero2,cero3,cero4,M1st,cero6,M3st2,cero8,M5st1,cero10,M7st2))

#### Tau3
dF6<-diag(c(cero1,cero2,cero3,cero4,cero5,cero6,cero7,M4st,M5st2,M6st2,M7st3))


dG1<--G%*%dF1%*%G
dG2<--G%*%dF2%*%G
dG3<--G%*%dF3%*%G
dG4<--G%*%dF4%*%G
dG5<--G%*%dF5%*%G
dG6<--G%*%dF6%*%G



#dGb<-c(unlist(dG1),unlist(dG2),unlist(dG3),unlist(dG4),unlist(dG5),unlist(dG6))
#dG<- array(dGb,c(dim(dG1)[1], dim(dG1)[2],6))



dV1<-Z%*%dG1%*%t(Z)
dV2<-Z%*%dG2%*%t(Z)
dV3<-Z%*%dG3%*%t(Z)
dV4<-Z%*%dG4%*%t(Z)
dV5<-Z%*%dG5%*%t(Z)
dV6<-Z%*%dG6%*%t(Z)


