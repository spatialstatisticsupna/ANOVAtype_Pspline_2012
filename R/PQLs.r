################################################################################
############## PQL PARA EL MODELO DE SPLINES ###################################



repeat {


param<-param


lambda1<-param[1]
lambda2<-param[2]
lambda3<-param[3]

tau1<-param[4]
tau2<-param[5]
tau3<-param[6]

source("COVs.r")

## MODELO DE SPLINES

J <- length(b)

M <- 1

    repeat {

      n.eta<-offset1 + X%*%a +Z%*%b

      mu.b<-exp(n.eta)

      log.mu.b <- log(mu.b)

      W <- diag(as.vector(mu.b)) 
      
      if(sum(as.vector(mu.b))==Inf|| prod(as.vector(mu.b))==0){break}

      W.inv <-  diag(1/as.vector(mu.b))  # ginv (W)  #solve(W)

      #if(sum(W.inv)==Inf){break}

      y.work <- n.eta + W.inv%*%as.matrix(y - as.vector(mu.b)) - offset1 ## Working vector

      y.work <- as.vector(y.work)

      V <- W.inv + Z %*% G %*% t(Z)

      V.inv <- solve (V)

      XVX <- t(X) %*% V.inv %*% X

      XVX.inv <- solve (XVX)

      XVY <- t(X) %*% V.inv %*% y.work

      P <- c(a,b)


      a.hat <- XVX.inv %*% XVY   #a.hat for the fixed efects

      b.hat <- G %*% t(Z) %*% V.inv %*% (y.work - X %*% a.hat) #b.hat for the randome efects


      P.hat <- c(as.vector(a.hat), as.vector(b.hat))

      dP <- sum(abs(P.hat - P))/length(P)
      
      
      a<-a.hat

      b<-b.hat

print(dP)

      if(dP < 0.0001){break}

      if(is.na(dP)==TRUE || dP > 10000 ){break}

      }


if(is.na(dP)==TRUE || dP > 10000 || sum(as.vector(mu.b))==Inf|| prod(as.vector(mu.b))==0){counter<-1000}else{


m <- exp(as.vector(a.hat))

lambda.ij <- exp(as.vector(b.hat[1:J, 1]))

mu.b <- exp(log(esp) + X %*% a.hat + Z %*% b.hat)

YXa <- y.work - X %*% a

y.hat <- mu.b/esp

res <- (y.hat - y)/sqrt(y.hat)



N.v<- length(param)

R<- length(y)

##dV <- array(0, c(R, R, N.v))

#    for(i in 1:N.v)
#    {
#    dV[, , i] <- Z %*% dG[, , i] %*% t(Z)
#    }
    
    
P <- V.inv - V.inv %*% X %*% XVX.inv %*% t(X) %*% V.inv

#tr.pdV <- 1:N.v

#U <- matrix(1, N.v, 1)

#    for(i in 1:N.v)
#    {
#    tr.pdV[i] <- sum((diag(P %*% dV[, , i])))

#    U[i, 1] <-  (1/2) * (t(YXa) %*% V.inv %*% dV[, , i]%*% V.inv %*% YXa - tr.pdV[i])
#    }

U <- matrix(1, N.v, 1)

U.part1<-t(YXa) %*% V.inv 
U.part2<-V.inv %*% YXa

PdV1<-P %*% dV1
PdV2<-P %*% dV2
PdV3<-P %*% dV3
PdV4<-P %*% dV4
PdV5<-P %*% dV5
PdV6<-P %*% dV6


U[1, 1] <-  (1/2) * (U.part1 %*% dV1%*% U.part2 - (sum((diag(PdV1)))))
U[2, 1] <-  (1/2) * (U.part1 %*% dV2%*% U.part2 - (sum((diag(PdV2)))))
U[3, 1] <-  (1/2) * (U.part1 %*% dV3%*% U.part2 - (sum((diag(PdV3)))))
U[4, 1] <-  (1/2) * (U.part1 %*% dV4%*% U.part2 - (sum((diag(PdV4)))))
U[5, 1] <-  (1/2) * (U.part1 %*% dV5%*% U.part2 - (sum((diag(PdV5)))))
U[6, 1] <-  (1/2) * (U.part1 %*% dV6%*% U.part2 - (sum((diag(PdV6)))))

G.s<-matrix(0, N.v, N.v)


G.s[1, 1]<- (1/2) * (sum(diag(PdV1 %*% PdV1)))
G.s[1, 2]<- (1/2) * (sum(diag(PdV1 %*% PdV2)))
G.s[1, 3]<- (1/2) * (sum(diag(PdV1 %*% PdV3)))
G.s[1, 4]<- (1/2) * (sum(diag(PdV1 %*% PdV4)))
G.s[1, 5]<- (1/2) * (sum(diag(PdV1 %*% PdV5)))
G.s[1, 6]<- (1/2) * (sum(diag(PdV1 %*% PdV6)))


G.s[2, 1]<- (1/2) * (sum(diag(PdV2 %*% PdV1)))
G.s[2, 2]<- (1/2) * (sum(diag(PdV2 %*% PdV2)))
G.s[2, 3]<- (1/2) * (sum(diag(PdV2 %*% PdV3)))
G.s[2, 4]<- (1/2) * (sum(diag(PdV2 %*% PdV4)))
G.s[2, 5]<- (1/2) * (sum(diag(PdV2 %*% PdV5)))
G.s[2, 6]<- (1/2) * (sum(diag(PdV2 %*% PdV6)))


G.s[3, 1]<- (1/2) * (sum(diag(PdV3 %*% PdV1)))
G.s[3, 2]<- (1/2) * (sum(diag(PdV3 %*% PdV2)))
G.s[3, 3]<- (1/2) * (sum(diag(PdV3 %*% PdV3)))
G.s[3, 4]<- (1/2) * (sum(diag(PdV3 %*% PdV4)))
G.s[3, 5]<- (1/2) * (sum(diag(PdV3 %*% PdV5)))
G.s[3, 6]<- (1/2) * (sum(diag(PdV3 %*% PdV6)))


G.s[4, 1]<- (1/2) * (sum(diag(PdV4 %*% PdV1)))
G.s[4, 2]<- (1/2) * (sum(diag(PdV4 %*% PdV2)))
G.s[4, 3]<- (1/2) * (sum(diag(PdV4 %*% PdV3)))
G.s[4, 4]<- (1/2) * (sum(diag(PdV4 %*% PdV4)))
G.s[4, 5]<- (1/2) * (sum(diag(PdV4 %*% PdV5)))
G.s[4, 6]<- (1/2) * (sum(diag(PdV4 %*% PdV6)))


G.s[5, 1]<- (1/2) * (sum(diag(PdV5 %*% PdV1)))
G.s[5, 2]<- (1/2) * (sum(diag(PdV5 %*% PdV2)))
G.s[5, 3]<- (1/2) * (sum(diag(PdV5 %*% PdV3)))
G.s[5, 4]<- (1/2) * (sum(diag(PdV5 %*% PdV4)))
G.s[5, 5]<- (1/2) * (sum(diag(PdV5 %*% PdV5)))
G.s[5, 6]<- (1/2) * (sum(diag(PdV5 %*% PdV6)))



G.s[6, 1]<- (1/2) * (sum(diag(PdV6 %*% PdV1)))
G.s[6, 2]<- (1/2) * (sum(diag(PdV6 %*% PdV2)))
G.s[6, 3]<- (1/2) * (sum(diag(PdV6 %*% PdV3)))
G.s[6, 4]<- (1/2) * (sum(diag(PdV6 %*% PdV4)))
G.s[6, 5]<- (1/2) * (sum(diag(PdV6 %*% PdV5)))
G.s[6, 6]<- (1/2) * (sum(diag(PdV6 %*% PdV6)))

#G.s<-matrix(0, N.v, N.v)
#    for(i in 1:N.v){
#        for(j in 1:N.v){
#         G.s[i, j]<- (1/2) * (sum(diag(P %*% dV[, ,i] %*% P %*% dV[, , j])))
#              }}

G.inv<-ginv(G.s)

Newparam<-param + (G.inv %*% U)/2

ds<-max(c(abs(Newparam - param), as.vector(abs(U))))

param<-Newparam

paramet<-c(a,b,Newparam)

#conv<-sum(abs(oldparam-paramet))
conv<-sum(abs(oldparam-paramet))/length(oldparam)

oldparam<-paramet

counter<-counter+1

print(conv)

print(counter)

}
if(conv < 0.0001 | counter == 1000 ){break}
}

















