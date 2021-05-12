
require("MASS")
#require("parallel")
#########################################   INPUT   #########################################
# data   : original time-series data matrix(row : time / col : variables)
# lag    : p in VAR(p)
# beta   : cointegrating vector(known)
# qn     : quantile of abs of error correction term Z (0~1)
# kn     : minimum number of observations in each regime
# m1     : number of grid for iterating when finding gamma hat
# nb     : number of bootstrapping

# g1,g2  : threshold parameters(gamma)
# Z      : error correction term matrix
# rhs    : deltax_t-1 term matrix(design matrix)
# lhs    : deltax_t term matrix

# gm     : gamma value(symmetric gamma : (-theta,theta) )
# ap     : alpha1 | alpha2 matrix
# mu     : level of the model
# t      : sample size

# rho    : VAR coefficients matrix(phi)

# ns     : number of simulation
# mod.type : "size" or "power" (specify the type of test)


###############################################################################################
###############################   Data simulation procedure    ################################

#  simulation code for all the cases
#  OUTPUT : Rejection ratio of three tests according to the input information 
simulation<-function(t0,ap0,gm0,rho0,mu,t,beta,nb,mod.type){ #mod.type =="size" or "power"
  if(t==100){
    m1 <- 0.1*t+16
  }else{
    m1 <- 0.06*t
  }
  
  if(mod.type=="size"){
    res10<-res5<-matrix(0,3,3)
    rownames(res10)<-c("supW","HW","ADF");colnames(res10)<-c("phi0_0.1","phi1_0.1","phi2_0.1")
    rownames(res5)<-c("supW","HW","ADF");colnames(res5)<-c("phi0_0.05","phi1_0.05","phi2_0.05")
    
    for(it in 1:length(t0)){
      ns<-t0[it]
      for(ir in 1:3){
        rho <- rho0[(2*(ir-1)+1):(2*ir),]
        
        sim.res<-simu.size(ns,mu,t,rho,beta,nb,m1)
        for(i in 1:3){
          res10[i,ir]<-sim.res[[i]][1,1]
          res5[i,ir]<- sim.res[[i]][1,2]
        }
      }
    }
    res.tot<-cbind(res10,res5)
    write.csv(res.tot,"Test of size_table3 result.csv")
    return(res.tot)
  }else{
    res10_1<-res5_1<-matrix(0,3,3);    res10_2<-res5_2<-matrix(0,3,3)
    rownames(res10_1)<-rownames(res5_1)<-rownames(res10_2)<-rownames(res5_2)<-c("supW","HW","ADF")
    colnames(res10_1)<-colnames(res10_2)<-c("gam0_0.1","gam1_0.1","gam2_0.1")
    colnames(res5_1)<-colnames(res5_2)<-c("gam0_0.05","gam1_0.05","gam2_0.05")
    res.list2<-list();res.list2<-list(Case1_alpha2=cbind(res10_1,res5_1),Case2_alpha2=cbind(res10_2,res5_2))
    
    for(it in 1:length(t0)){
      ns<-t0[it]
      for(ia in 1:(nrow(ap0)/2)){
        ap <- ap0[(2*(ia-1)+1):(2*ia),]
        for(ig in 1:length(gm0)){
          gm<-gm0[ig]
          
          sim.res<-simu.pow(ns,gm,ap,mu,t,beta,nb,m1)
          for(i in 1:3){
            res.list2[[ia]][i,ig]<-sim.res[[i]][1,1]
            res.list2[[ia]][i,(3+ig)]<-sim.res[[i]][1,2]
          }
          #mat.num<-c(it,ia,ig)
          #print(mat.num);print(sim.res)
        }
      }
    }
    res.tot<- rbind(matrix(c("Case1","-","-","-","-","-"),1,6),res.list2[[1]]
                    ,matrix(c("Case2","-","-","-","-","-"),1,6),res.list2[[2]]) 
    write.csv(res.tot,"Test of power_table4 result.csv")
    return(res.tot)
  }
}

#  simulation code for each case
#  OUTPUT : Rejection ratio of boot.sup-wald/HW/ADF tests under NULL
simu.size<-function(ns,mu,t,rho,beta,nb,m1){  
  wn<-adf1<-matrix(0,ns,1)
  bp1<-matrix(0,ns,3)
  lapply(1:ns, function(i){
    list<-dgp_linear(rho,t);  xs<-as.matrix(list$x)
    bp1[i,]<<- boot.w(xs,1,beta,qn,kn,m1,nb)
    wn[i,] <<- hwcoint(xs,1,beta)
    z<-xs[,1]-xs[,2]
    adf1[i,]<<-adft(z,1)
  })
  
  rej.bp<-matrix(0,3,2);rej.hw<-rej.adft<-matrix(0,1,2)
  colnames(rej.bp)<-colnames(rej.hw)<-colnames(rej.adft)<-c("0.1","0.05")
  rownames(rej.bp)<-c("joint","supw1","supw2");rownames(rej.hw)<-"HW";rownames(rej.adft)<-"ADF"
  
  for(j in 1:ncol(bp1)){
    rej.bp[j,]<-c(sum(bp1[,j]< 0.1)/ns,sum(bp1[,j]<0.05)/ns)
  } 
  rej.hw[1,]<-c((sum(wn>8.30)/ns),(sum(wn>10.18)/ns))
  rej.adft[1,]<-c((sum(adf1<(-2.59))/ns),(sum(adf1<(-2.89))/ns))
  
  list.res<-list(Boot.sup = rej.bp,H.W=rej.hw,ADF=rej.adft)
  return(list.res)
}
#  OUTPUT : Rejection ratio of boot.sup-wald/HW/ADF tests under Alternative
simu.pow<-function(ns,gm,ap,mu,t,beta,nb,m1){
  wn<-adf1<-matrix(0,ns,1)
  bp1<-matrix(0,ns,3)
  
  lapply(1:ns, function(i){
    xs <- dgp_band(gm,ap,mu,t)
    bp1[i,]<<- boot.w(xs,1,beta,qn,kn,m1,nb)
    wn[i,] <<- hwcoint(xs,1,beta)
    z<-xs[,1]-xs[,2]
    adf1[i,]<<-adft(z,1)
  })
  
  rej.bp<-matrix(0,3,2);rej.hw<-rej.adft<-matrix(0,1,2)
  colnames(rej.bp)<-colnames(rej.hw)<-colnames(rej.adft)<-c("0.1","0.05")
  rownames(rej.bp)<-c("joint","supw1","supw2");rownames(rej.hw)<-"HW";rownames(rej.adft)<-"ADF"
  
  for(j in 1:ncol(bp1)){
    rej.bp[j,]<-c(sum(bp1[,j]< 0.1)/ns,sum(bp1[,j]<0.05)/ns)
  } 
  rej.hw[1,]<-c((sum(wn>8.30)/ns),(sum(wn>10.18)/ns))  
  rej.adft[1,]<-c((sum(adf1<(-2.59))/ns),(sum(adf1<(-2.89))/ns))  
  list.res<-list(Boot.sup = rej.bp,H.W=rej.hw,ADF=rej.adft)
  
  return(list.res)
}  





###############################################################################################
#############################   Data construction for simulation   ############################

# 1. dxt=mu+alpha%*%ZU ZL + e
# OUTPUT : Band-TVECM based random generated data matrix
dgp_band<-function(gm,ap,mu,t){
  t1<-t
  u<-mvrnorm(t1,c(0,0),diag(1,2))
  dxi<-matrix(0,2,1)
  xi<-matrix(0,2,1)
  dx<-matrix(0,t1,2)
  x<-matrix(0,t1,2)
  for(i in 1:t1){
    zi<-xi[1,1]-xi[2,1]
    zl <- ifelse(zi <= (-gm), 1, 0) *zi
    zu <- ifelse(zi > gm, 1, 0) * zi
    dxi<-mu+zl*t(ap[1,,drop=F])+zu*t(ap[2,,drop=F])+t(u[i,,drop=F]) #2x1
    xi<-xi+dxi
    dx[i,]<-dxi
    x[i,]<-xi
  }
  return(x)
}

# 2. dxt=phi%*%dxt_1+e
# OUTPUT : VAR based random generated data matrix
dgp_linear<-function(rho,t){
  u <- mvrnorm(t,c(0,0),diag(1,2))
  dxi<-matrix(0,2,1)
  xi<-matrix(0,2,1)
  dx<-matrix(0,t,2)
  x<-matrix(0,t,2)
  
  for(i in 1:t){
    dxi<-rho%*%dxi+t(u[i,,drop=F]) #2x1
    xi<-xi+dxi              #2x1
    dx[i,]<-t(dxi)
    x[i,]<-t(xi)
  }
  
  list(x=x,dx=dx)
}




###############################################################################################
###############################   Derive test_statistics   ####################################

### 1. Bootstraped Sup-W p-value function  
# 1) Residual-based Bootstrap
# OUTPUT : p-value of bootstrap-based sup-wald test for joint and marginal sup-wald stat.
boot.w<-function(data, lag, beta ,qn,kn,m1,nb){
  #Construct data matrix and major numeric
  x<-as.matrix(data)
  n<-nrow(x)
  p<<-ncol(x)
  k<-lag
  
  #Construct VECM data matrix
  dx_lhs<-embed(diff(x),k+1)[,1:p]    # deltax t(1,...,p)
  dx_lag<-embed(diff(x),k+1)[,-(1:p)] # deltax t-1(1,...,p) ~ deltax t-k(1,...,p)
  dx_rhs<-cbind(1,dx_lag)
  
  if(k>0){
    Z<-(x%*%beta)[-c(1:k,n)]
  }else{
    Z<-(x%*%beta)[-n]
  }
  
  #Estimation from gamma.hat : argmin det(Sigma)   
  list.est<-est.w(Z,dx_rhs,dx_lhs, beta ,qn,kn,m1)
  sup<-as.numeric(list.est[[1]])
  gamma.hat<-as.numeric(list.est[[2]][2:3]) #2x1
  Z.hat <- cbind(ifelse(Z <= gamma.hat[1], 1, 0) * Z,ifelse(Z > gamma.hat[2], 1, 0) * Z) #n-k-1 x 2
  X.hat <- cbind(Z.hat, dx_rhs) #n-k-1 x (2+1+pk)
  para.hat <- solve(crossprod(X.hat))%*%t(X.hat)%*%dx_lhs #(2+1+pk)xp
  res.hat <- dx_lhs - X.hat%*%para.hat #n-k-1 x p
  #result.hat <- find.w(gamma.hat[1], gamma.hat[2], Z, dx_rhs, dx_lhs,kn)
  
  #Derive bootstrapped supw distribution 
  supw.b<-matrix(0,nb,3)
  
  list.boot<-lapply(1:nb, function(xxx){
    #Define bootstrap based deltax, x, and Z
    if(k>0){
      dx_lhs.b <- matrix(0, (n-1), p) #n-1 x p
      dx_lhs.b[1:k, ] <- diff(x)[1:k, ]
      X.b <- matrix(0,n,p)            #n x p
      X.b[1:(k+1), ] <- x[1:(k+1), ]  #knowing deltax1 == knowing x1,x2
    }else{
      dx_lhs.b <- matrix(0, (n-1), p) #n-1 x p
      X.b <- matrix(0,n,p)            #n x p
      X.b[1:(k+1), ] <- x[1:(k+1), ]  #knowing deltax1 == knowing x1,x2
    }
    
    #generate residual resampling data
    sample<-sample(1:nrow(res.hat), nrow(res.hat), replace=TRUE)
    res.boot.sam<-res.hat[sample,]  #n-k-1 x p
    #Construct bootstrap based deltax, x
    if(k>0){
      for(i in 1:(n-k-1)){
        tdx<-t(dx_lhs.b[i:(k+i-1),,drop=F]) #tdx : [deltax1,...,deltaxk] : pxk
        tdx.tr<-tdx
        for(j in 1:k){
          tdx.tr[,j]<-tdx[,(k+j-1)] #tdx.tr : [deltaxk, ... , deltax1] :pxk
        }
        dx_lhs.b[(k+i),]<-matrix(tdx.tr,nrow=1)%*%para.hat[-c(1:3),]+res.boot.sam[i,]    #1xp
        X.b[(k+i+1),]<-X.b[(k+i),]+dx_lhs.b[(k+i),] # 1xp
      }
    }else{
      for(i in 1:(n-k-1)){
        dx_lhs.b[(k+i),]<-res.boot.sam[i,]    #1xp
        X.b[(k+i+1),]<-X.b[(k+i),]+dx_lhs.b[(k+i),] # 1xp
      }
    }
    
    #Construct bootstrapped VECM data matrixes
    dx_lhs.b<-embed(diff(X.b),k+1)[,1:p]      # deltax t(1,...,p)
    dx_lag.b<-embed(diff(X.b),k+1)[,-(1:p)] # deltax t-1(1,...,p) ~ deltax t-k(1,...,p)
    dx_rhs.b<-cbind(1,dx_lag.b)
    
    if(k>0){
      Z.b<-(X.b%*%beta)[-c(1:k,n)]
    }else{
      Z.b<-(X.b%*%beta)[-n]
    }
    
    #Construct bootstrapped Threshold parameter gamma
    sprt.b <- sort(unique(abs(Z.b)))
    sprt.b <- sprt.b[floor(qn*length(Z.b))]
    sprt.b <- seq(from=(-sprt.b),by=((2*sprt.b)/(m1-1)),length.out=m1)
    n.gam.b<- length(sprt.b)
    
    #find bootstrapped Sup-Wald among gamma para.space / argmin sigma.hat 
    sw.b<-matrix(0,n.gam.b,n.gam.b)
    sw1.b<-matrix(0,n.gam.b,n.gam.b)
    sw2.b<-matrix(0,n.gam.b,n.gam.b)
    
    for(i in 1:n.gam.b){
      g1.b<-sprt.b[i]
      for(j in i:n.gam.b){
        g2.b<-sprt.b[j]
        result.b <- try(find.w(g1.b,g2.b,Z.b,dx_rhs.b,dx_lhs.b,kn), silent = TRUE)
        if (inherits(result.b, "try-error")) {
          sw.b[i, j] <- sw1.b[i, j] <- sw2.b[i,j]<-NA
        }else {
          sw.b[i, j] <- result.b$Wald
          sw1.b[i, j] <- result.b$W1
          sw2.b[i, j] <- result.b$W2
        }
      }
    }
    return(c(max(sw.b,na.rm = T),max(sw1.b,na.rm = T),max(sw2.b,na.rm = T)))
  })
  for( i in 1:nb ){
    supw.b[i,]<-list.boot[[i]]
  }
  
  
  #Calculate p-value and critical value from bootstrapped asymptotic distribution of supw
  pval.b <- c(mean(ifelse(supw.b[,1] > sup[1], 1, 0)),
              mean(ifelse(supw.b[,2] > sup[2], 1, 0)),
              mean(ifelse(supw.b[,3] > sup[3], 1, 0))  )
  
  #generate results
  #conseq <- list(supWald=sup[1],gamma1=gam.sup1,gamma2=gam.sup2,Parameter.hat=para.hat,Wald.boots = supw.b,
  #               Pvalue.Boot=pval.b)
  return(pval.b)
  
}

# 2) Derive Sup-Wald statistics
# OUTPUT - joint/marginal sup-wald statistics 
#        - min. of det. of variance matrix of residual
#        - gamma which minimizes determinant of sigma
est.w<-function(Z,dx_rhs,dx_lhs, beta ,qn,kn,m1){
  
  #Construct Threshold parameter gamma
  sprt<-sort(unique(abs(Z)))
  sprt <- sprt[floor(qn*length(sprt))]
  sprt <- seq(from=-sprt,by=2*sprt/(m1-1),length.out=m1)
  n.gam<-length(sprt)
  
  #find Sup-Wald among gamma para.space / argmin sigma.hat 
  sw<-matrix(0,n.gam,n.gam)
  sw1<-matrix(0,n.gam,n.gam)
  sw2<-matrix(0,n.gam,n.gam)
  ssig<-matrix(NA,n.gam,n.gam)
  
  for(i in 1:n.gam){
    g1<-sprt[i]
    for(j in i:n.gam){
      g2<-sprt[j]
      result <- try(find.w(g1,g2,Z,dx_rhs,dx_lhs,kn), silent = TRUE)
      if (inherits(result, "try-error")) {
        sw[i, j] <- sw1[i, j] <- sw2[i,j] <- ssig[i,j]<-NA
      } else {
        sw[i, j] <- result$Wald
        sw1[i, j] <- result$W1
        sw2[i, j] <- result$W2
        ssig[i, j] <- result$det
      }
      
    }
  }
  sup <- c(max(sw,na.rm = T),max(sw1,na.rm = T),max(sw2,na.rm = T))
  argsup<- which(sw==sup[1],arr.ind = T)
  gam.sup1<-sprt[argsup[1]]
  gam.sup2<-sprt[argsup[2]]
  
  min.sigma<- min(ssig,na.rm = T)
  arginf<- which(ssig==min.sigma,arr.ind = T)
  gam.inf1<-sprt[arginf[1]]
  gam.inf2<-sprt[arginf[2]]
  
  res.list1<-list("Sup.Wald"=sup,"Min.sigma&gamma"=c(min.sigma,gam.inf1,gam.inf2))
}

# 3) Find Wald statistics
# OUTPUT : list of joint/marginal wald statistics and det. of variance matrix of residual
find.w<-function(g1,g2,Z,rhs,lhs,kn){
  if(sum(ifelse(Z <= g1, 1, 0))>kn & sum(ifelse(Z > g2, 1, 0))>kn ){
    Z.gam<-cbind(ifelse(Z <= g1, 1, 0) * Z,ifelse(Z > g2, 1, 0) * Z) #n-k-1 x2
    X.tvecm<-cbind(Z.gam,rhs) #n-k-1 x (3+pk)
    bhat<- solve(crossprod(X.tvecm))%*%t(X.tvecm)%*%lhs #3+pk x p
    resi <-lhs - X.tvecm%*%bhat  #n-k-1 x p
    sig <-crossprod(resi)/nrow(resi) #pxp
    invrs <- solve(crossprod(X.tvecm))%x%sig
    
    sseij<-det(sig)
    w<- matrix(t(bhat[1:2,]), nrow = 1)%*%solve(invrs[1:(2*p),1:(2*p)])%*% matrix(t(bhat[1:2,]), ncol = 1)
    w1<- matrix(t(bhat[1,]), nrow = 1)%*%solve(invrs[1:p,1:p])%*% matrix(t(bhat[1,]), ncol = 1)
    w2<- matrix(t(bhat[2,]), nrow = 1)%*%solve(invrs[(p+1):(2*p),(p+1):(2*p)])%*% matrix(t(bhat[2,]), ncol = 1)
  }
  list(Wald = w, W1 = w1, W2 = w2, det = sseij)
}



### 2. HW test statistics
hwcoint<-function(data,lag,beta){
  #Construct data matrix and major numerics
  x<-as.matrix(data)
  n<-nrow(x)
  m<-ncol(x)
  p<-lag
  
  #Construct VECM data matrix
  dx_lhs<-embed(diff(x),p+1)[,1:m]    # deltax t(1,...,m)
  dx_lag<-embed(diff(x),p+1)[,-(1:m)] # deltax t-1(1,...,m) ~ deltax t-p(1,...,m)
  dx_rhs<-cbind(1,dx_lag)
  
  Z<-(x%*%beta)[-c(1:p,n)] #n-p-1 x 1
  
  M<-diag(1,n-p-1)-dx_rhs%*%solve(t(dx_rhs)%*%dx_rhs)%*%t(dx_rhs) #n-p-1 x n-p-1 / I - pi under H0(alpha=0)
  My<-M%*%dx_lhs #n-p-1 x m  / (I-pi)Y = e  under H0(alpha=0)
  Mz<-M%*%Z  #n-p-1 x 1  
  
  #estimate alpha.hat, residual, sigma and Wald stat.
  ahat<-solve(crossprod(Mz))%*%crossprod(Mz,My) #1xm     #####why solve(z'M2z)z'M2y? not z'mz z'my?
  res<-My-Mz%*%ahat #n-p-1 x m  /residual
  sig<-crossprod(res)/nrow(res) #mxm / var(res)
  Wn<- t(matrix(ahat, ncol = 1))%*%solve(solve( crossprod(Mz) ) %x% sig)%*% matrix(ahat, ncol = 1)
  
  return(Wn)
}



### 3. ADF test statistics
adft<-function(data,lag){
  #Construct data matrix and major numerics
  x<-as.matrix(data)
  n<-nrow(x)
  m<-ncol(x)         ########always 1?
  p<-lag
  
  #Construct VECM data matrix
  dx_lhs<-embed(diff(x),p+1)[,1:m]    # deltax t(1,...,m) (n-p-1 x m)
  dx_lag<-embed(diff(x),p+1)[,-(1:m)] # deltax t-1(1,...,m) ~ deltax t-p(1,...,m) (n-p-1 x mp)
  dx_rhs<-cbind(x[-c(1:p,n),],1,dx_lag) #(n-p-1 x m(p+1)+1)
  
  irhs<-solve(crossprod(dx_rhs))  #m(p+1)+1 x m(p+1)+1
  btdf<-irhs%*%t(dx_rhs)%*%dx_lhs #m(p+1)+1 x m   /inv(X'X)%*%X'%*%y = beta.hat
  ep<-dx_lhs-dx_rhs%*%btdf        #n-p-1 x m  /y-X%*%beta.hat = residual
  sg<-crossprod(ep)/nrow(ep)      #mxm  /var(res).hat
  adfbi<-btdf[1,]/sqrt(sg*(irhs[1,1])) # var(beta.hat)
  
  return(adfbi)
}




