#install.packages("MASS")
#library(MASS)
#help("mvrnorm")
require("MASS")
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
# rhs/dx_rhs    : deltax_t-1 term matrix(design matrix)
# lhs/dx_lhs    : deltax_t term matrix

## Data import
data<-get.data("se46.csv")

# Setting parameters
lag<-0;beta<-c(1,-1);qn<-0.9;kn<-10;m1<-100;nb<-10

# Main Function
boot.w.detail(data,lag,beta,qn,kn,m1,nb)



# Residual-based Bootstrap
# OUTPUT : Sup-Wald, Gamma.hat, Coefficient.hat, P-value.Boot,
#          Minimum number of obs. in each regime(kn), grid number(m1)
boot.w.detail<-function(data, lag, beta ,qn,kn,m1,nb){
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
  gamma.res<-c( paste0(round(gamma.hat[1],4),"(",list.est[[2]][4],")"),paste0(round(gamma.hat[2],4),"(",list.est[[2]][5],")") )
  gamma.mat<-data.frame();gamma.mat[1,1]<-gamma.res[1];gamma.mat[2,1]<-gamma.res[2]
  rownames(gamma.mat)<-c("gamma1","gamma2");colnames(gamma.mat)<-"value(qunatile)"
  if(k>0){
    rownames(para.hat)<-c("Coefficient of ECT1","Coefficient of ECT2","Constants",
                          paste( rep("Coefficients of lag",(p*k)) ,rep(1:k,each=p),rep("for",(p*k)),rep(1:p,k),rep("-th variable",(p*k)) ) )
  }else{
    rownames(para.hat)<-c("Coefficient of ECT1","Coefficient of ECT2","Constants")
  }
  
  colnames(para.hat)<-paste0( rep("X",p),1:p  )
  pval.b<-matrix(pval.b,3,1);rownames(pval.b)<-c("Joint sup","1st-marginal sup","2nd-marginal sup")
  
  conseq <- list("Sup-Wald stat."=sup[1],
                 "Gamma hat and its quantile in ECT values"=gamma.mat,
                 "Estimation of Coefficients"=para.hat,
                 "P-value of Bootstrap-based Sup-wald stat."=pval.b,
                 "Minimum number of observation in each regime"=kn,
                 "Number of iterating grid"=m1)
  #lapply(conseq, function(x) write.table( data.frame(x), 'Result of bootstrapping_se46.csv'  , append= T, sep=',' ))
  print(conseq)
  
}

# Derive Sup-Wald statistics
# OUTPUT - joint/marginal sup-wald statistics 
#        - min. of det. of variance matrix of residual
#        - gamma which minimizes determinant of sigma
#        - empirical cdf value of gamma1 and gamma2 in the ECT values
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
      }
      else {
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
  quan1<-round(sum(ifelse(Z <= gam.inf1, 1, 0))/length(Z),3) 
  quan2<-round(sum(ifelse(Z <= gam.inf2, 1, 0))/length(Z),3)
  
  res.list1<-list("Sup.Wald"=sup,"Min.sigma&gamma"=c(min.sigma,gam.inf1,gam.inf2,quan1,quan2))
}

# Find Wald statistics
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

# Data importing procedure
# OUTPUT : Warning message if there is no se46.csv file in working dir. or se46 data matrix
get.data<-function(file.name){
  if(file.access(file.name,0)==-1){
    print(paste("warning - there is no '",file.name,
                "'file in your working directory. please get the '",
                file.name,"'file in your working directory :",getwd()) )
  }else{
    data.o<-read.csv(file.name,header = F)
    data<-data.o[,c(1,14)]
    print("Data importing has been done successfully")
  }
  return(data)
}
