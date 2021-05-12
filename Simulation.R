
#########################################   INPUT   #########################################
# t0   : simulation number
# t    : sample size (100 or 250)
# nb   : bootstrapping number 
# mu   : level of the model
# beta : cointegrating vector(known)
# ap0(alpha) : threshold coefficients
# gm0(gamma) : threshold parameters
# rho0(phi)  : VAR coefficients matrix
# mod.type   : "size" or "power" (Size : table3 / Power : table4)

#########################################   OUTPUT   ########################################
# p-value of each cases according to significant level / alpha / gamma / rho / mod.type



###############################    PARAMETER VALUE SETTING   ################################
t0<-1              
t<-100             
nb<-1        
qn<-1;kn<-10
mu<-matrix(c(0,0),2,1) 
beta<-c(1,-1)          
ap0<-rbind(matrix(c(-0.1,0,0,0.1),2,2),matrix(c(-0.1,-0.3,0,0),2,2))  
# 1) when alpha2 is -0.3, what is alpha1? -0.1 or ?
# 2) In the paper, alpha1 is -0.1, alpha2 is 0.1 but in the gauss code, alpha1 is 0.1, alpha2 is -0.1
#    What is the right value?
gm0<-c(5,8,10)  
rho0 <- rbind(matrix(0,2,2),matrix(c(-0.2,-0.1,0,-0.2),2,2),matrix(c(-0.2,-0.1,-0.1,-0.2),2,2))
mod.type<-"size" #"power"


###################################      SIMULATION     #####################################
getwd()
source("Functions_lapply.R")
simulation(t0,ap0,gm0,rho0,mu,t,beta,nb,mod.type)
system.time(simulation(t0,ap0,gm0,rho0,mu,t,beta,nb,mod.type))




