#kasus 1
#fungsi log likelihood
loglik=expression(log(lambda*exp(-lambda*x)))
#turunan pertama
dloglik=D(loglik,"lambda")
dloglik
#turunan kedua
ddloglik=D(dloglik,"lambda")
ddloglik

#contoh 1
#fungsi mle
mle.exp=function(x,lambda0){
  ##fungsi turunan pertama
  dloglik=function(x,lambda){
    return(sum((exp(-lambda * x) - lambda * (exp(-lambda * x) * x))/(lambda * 
                                                                       exp(-lambda * x))))
  }
  ##fungsi turunan kedua
  ddloglik=function(x,lambda){
    return(sum(-((exp(-lambda * x) * x + ((exp(-lambda * x) * x) - lambda * 
                                            (exp(-lambda * x) * x * x)))/(lambda * exp(-lambda * x)) + 
                   (exp(-lambda * x) - lambda * (exp(-lambda * x) * x)) * (exp(-lambda * 
                                                                                 x) - lambda * (exp(-lambda * x) * x))/(lambda * exp(-lambda * 
                                                                                                                                       x))^2)))
  }
  #batas toleransi
  tol=10^(-6)
  #inisiasi
  iter=1
  lambda=lambda0
  #lambda hat
  lambda.hat=lambda-dloglik(x,lambda)/ddloglik(x,lambda)
  #proses iteratif
  while (abs(lambda.hat-lambda)>tol) {
    lambda=lambda.hat
    lambda.hat=lambda-dloglik(x,lambda)/ddloglik(x,lambda)
    iter=iter+1
  }
  #output
  cat('Nilai estimasi lambda =',lambda.hat,'\n')
  cat('Nilai diperoleh pada iterasi ke',iter,'\n')
}

#data
data=c(0.038,0.177,0.131,0.044,0.048,0.034,0.088,0.086,0.298,0.383,0.283, 0.092,0.044,0.009,0.075,0.132,0.013,0.585,0.180,0.101)
#nilai awal = 7
mle.exp(data,7)
#nilai awal = rata - rata sampel
mle.exp(data,mean(data))
#nilai awal = 16
mle.exp(data,16)

#cek estimasi parameter
library(MASS)
fitdistr(data,'exponential')


#contoh 2
#fungsi log likelihood
loglik=expression(log(exp(-mu)*(mu^(x))/factorial(x)))
#turunan pertama
dloglik=D(loglik,"mu")
dloglik
#turunan kedua
ddloglik=D(dloglik,"mu")
ddloglik
#fungsi mle
mle.poi=function(x,mu0){
  ##fungsi turunan pertama
  dloglik=function(x,mu){
    return(sum((exp(-mu) * (mu^((x) - 1) * (x)) - exp(-mu) * (mu^(x)))/factorial(x)/(exp(-mu) * 
                                                                                       (mu^(x))/factorial(x))))
  }
  
  ##fungsi turunan kedua
  ddloglik=function(x,mu){
    return(sum((exp(-mu) * (mu^(((x) - 1) - 1) * ((x) - 1) * (x)) - exp(-mu) * 
                  (mu^((x) - 1) * (x)) - (exp(-mu) * (mu^((x) - 1) * (x)) - 
                                            exp(-mu) * (mu^(x))))/factorial(x)/(exp(-mu) * (mu^(x))/factorial(x)) - 
                 (exp(-mu) * (mu^((x) - 1) * (x)) - exp(-mu) * (mu^(x)))/factorial(x) * 
                 ((exp(-mu) * (mu^((x) - 1) * (x)) - exp(-mu) * (mu^(x)))/factorial(x))/(exp(-mu) * 
                                                                                           (mu^(x))/factorial(x))^2))
  }
  #batas toleransi
  tol=10^(-6)
  #inisiasi
  iter=1
  mu=mu0
  #mu hat
  mu.hat=mu-dloglik(x,mu)/ddloglik(x,mu)
  #proses iteratif
  while(abs(mu.hat-mu)>tol){
    mu=mu.hat
    mu.hat=mu-dloglik(x,mu)/ddloglik(x,mu)
    iter=iter+1
  }
  #output
  cat('Nilai estimasi mu =',mu.hat,'\n')
  cat('Nilai diperoleh pada iterasi ke',iter,'\n')
}

#data
set.seed(77)
x=rpois(25,5)
x

#nilai awal = 4
mle.poi(x,4)
#nilai awal = rata - rata sampel
mle.poi(x,mean(x))
#nilai awal = 12
mle.poi(x,12)

#cek estimasi parameter
library(MASS)
fitdistr(x,'poisson')


################################
# kasus 2
#cara komputasi turunannya
loglik=expression(log(1/(gamma(a)*b^a)*x^(a-1)*exp(-x/b)))
dla=D(loglik,"a")
dla
ddla=D(dla,"a")
ddla
dlb=D(loglik,"b")
dlb
ddlb=D(dlb,"b")
ddlb
ddlab=D(dla,"b")
ddlab

#cara membuat fungsi MLE
mle.gamma=function(x,a0,b0){
  dla=function(x,a,b){
    return(sum((1/(gamma(a) * b^a) * (x^(a - 1) * log(x)) - (gamma(a) *
                                                               digamma(a)*b^a + gamma(a) * (b^a * log(b)))/(gamma(a) * b^a)^2 * x^(a -
                                                                                                                                     1)) * exp(-x/b)/(1/(gamma(a) * b^a) * x^(a - 1) * exp(-x/b))))
  }
  dlb=function(x,a,b){
    return(sum((1/(gamma(a) * b^a) * x^(a - 1) * (exp(-x/b) * (x/b^2)) -
                  gamma(a)*(b^(a - 1) * a)/(gamma(a) * b^a)^2 * x^(a - 1) * exp(-
                                                                                  x/b))/(1/(gamma(a) * b^a) * x^(a - 1) * exp(-x/b))))
  }
  ddla=function(x,a,b){
    return(sum((1/(gamma(a) * b^a) * (x^(a - 1) * log(x) * log(x)) -
                  (gamma(a) * digamma(a) * b^a + gamma(a) * (b^a * log(b)))/(gamma(a) *
                                                                               b^a)^2 * (x^(a - 1) * log(x)) - ((((gamma(a) * digamma(a) * digamma(a)
                                                                                                                   + gamma(a) * trigamma(a)) * b^a + gamma(a) * digamma(a) * (b^a * log(b))
                                                                                                                  + (gamma(a) * digamma(a) * (b^a * log(b)) + gamma(a) * (b^a * log(b) *
                                                                                                                                                                            log(b))))/(gamma(a) * b^a)^2 - (gamma(a) * digamma(a) * b^a + gamma(a)
                                                                                                                                                                                                            * (b^a * log(b))) * (2 * ((gamma(a) * digamma(a) * b^a + gamma(a) * (b^a
                                                                                                                                                                                                                                                                                 * log(b))) * (gamma(a) * b^a)))/((gamma(a) * b^a)^2)^2) * x^(a - 1) +
                                                                                                                  (gamma(a) * digamma(a) * b^a + gamma(a) * (b^a * log(b)))/(gamma(a) *
                                                                                                                                                                               b^a)^2 * (x^(a - 1) * log(x)))) * exp(- x/b)/(1/(gamma(a) * b^a) * x^(a
                                                                                                                                                                                                                                                     - 1) * exp(-x/b)) - (1/(gamma(a) * b^a) * (x^(a - 1) * log(x)) - (gamma(a)
                                                                                                                                                                                                                                                                                                                       * digamma(a) * b^a + gamma(a) * (b^a * log(b)))/(gamma(a) * b^a)^2 *
                                                                                                                                                                                                                                                                            x^(a - 1)) * exp(-x/b) * ((1/(gamma(a) * b^a) * (x^(a - 1) * log(x)) -
                                                                                                                                                                                                                                                                                                         (gamma(a) * digamma(a) * b^a + gamma(a) * (b^a * log(b)))/(gamma(a) *
                                                                                                                                                                                                                                                                                                                                                                      b^a)^2 * x^(a - 1)) * exp(-x/b))/(1/(gamma(a) * b^a) *
                                                                                                                                                                                                                                                                                                                                                                                                          x^(a - 1) * exp(-x/b))^2))
  }
  ddlb=function(x,a,b){
    return(sum((1/(gamma(a) * b^a) * x^(a - 1) * (exp(-x/b) * (x/b^2) *
                                                    (x/b^2) - exp(-x/b) * (x * (2 * b)/(b^2)^2)) - gamma(a) * (b^(a - 1) *
                                                                                                                 a)/(gamma(a) * b^a)^2 * x^(a - 1) * (exp(-x/b) * (x/b^2)) - ((gamma(a)
                                                                                                                                                                               * (b^((a - 1) - 1) * (a - 1) * a)/(gamma(a) * b^a)^2 - gamma(a) * (b^(a
                                                                                                                                                                                                                                                     - 1) * a) * (2 * (gamma(a) * (b^(a - 1) * a) * (gamma(a) *
                                                                                                                                                                                                                                                                                                       b^a)))/((gamma(a) * b^a)^2)^2) * x^(a - 1) * exp(-x/b) + gamma(a) * (b^(a
                                                                                                                                                                                                                                                                                                                                                                               - 1) * a)/(gamma(a) * b^a)^2 * x^(a - 1) * (exp(-x/b) *
                                                                                                                                                                                                                                                                                                                                                                                                                             (x/b^2))))/(1/(gamma(a) * b^a) * x^(a - 1) * exp(-x/b)) - (1/(gamma(a)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           * b^a) * x^(a - 1) * (exp(-x/b) * (x/b^2)) - gamma(a) * (b^(a - 1) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      a)/(gamma(a) * b^a)^2 * x^(a - 1) * exp(-x/b)) * (1/(gamma(a) * b^a) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          x^(a - 1) * (exp(-x/b) * (x/b^2)) - gamma(a) * (b^(a - 1) * a)/(gamma(a)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          * b^a)^2 * x^(a - 1) * exp(-x/b))/(1/(gamma(a) * b^a) * x^(a - 1) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               exp(-x/b))^2))
  }
  ddlab=ddlba=function(x,a,b){
    return(sum(((1/(gamma(a) * b^a) * (x^(a - 1) * log(x)) - (gamma(a) *
                                                                digamma(a)*b^a + gamma(a) * (b^a * log(b)))/(gamma(a) * b^a)^2 * x^(a -
                                                                                                                                      1)) * (exp(-x/b) * (x/b^2)) - (gamma(a) * (b^(a - 1) * a)/(gamma(a) *
                                                                                                                                                                                                   b^a)^2 * (x^(a - 1) * log(x)) + ((gamma(a) * digamma(a) *(b^(a - 1) *
                                                                                                                                                                                                                                                               a) + gamma(a) * (b^(a - 1) * a * log(b) + b^a * (1/b)))/(gamma(a) *
                                                                                                                                                                                                                                                                                                                          b^a)^2 - (gamma(a) * digamma(a) * b^a + gamma(a) * (b^a * log(b))) * (2
                                                                                                                                                                                                                                                                                                                                                                                                * (gamma(a) * (b^(a - 1) * a) * (gamma(a) * b^a)))/((gamma(a) * b^a)^2)^2)
                                                                                                                                                                     * x^(a - 1)) * exp(-x/b))/(1/(gamma(a) * b^a) * x^(a - 1) * exp(-x/b))
               - (1/(gamma(a) * b^a) * (x^(a - 1) * log(x)) - (gamma(a) * digamma(a) *
                                                                 b^a + gamma(a) * (b^a * log(b)))/(gamma(a) * b^a)^2 * x^(a - 1)) * exp(-
                                                                                                                                          x/b) * (1/(gamma(a) * b^a) * x^(a - 1) * (exp(-x/b) * (x/b^2)) - gamma(a)
                                                                                                                                                  * (b^(a - 1) * a)/(gamma(a) * b^a)^2 * x^(a - 1) * exp(-
                                                                                                                                                                                                           x/b))/(1/(gamma(a) * b^a) * x^(a - 1) * exp(-x/b))^2))
  }
  tol=10^(-6)
  iter=1
  teta=matrix(c(a0,b0),nrow=2)
  a=teta[1,1]
  b=teta[2,1]
  Gm=matrix(c(dla(x,a,b),dlb(x,a,b)),nrow=2)
  Hm=matrix(c(ddla(x,a,b),ddlba(x,a,b),ddlab(x,a,b),ddlb(x,a,b)),nrow=2,
            ncol=2)
  teta.hat=teta-solve(Hm)%*%Gm
  galat=abs(teta.hat-teta)
  while (any(galat>tol)) {
    teta=teta.hat
    a=teta[1,1]
    b=teta[2,1]
    Gm=matrix(c(dla(x,a,b),dlb(x,a,b)),nrow=2)
    Hm=matrix(c(ddla(x,a,b),ddlba(x,a,b),ddlab(x,a,b),ddlb(x,a,b)),nrow=2,
              ncol=2)
    teta.hat=teta-solve(Hm)%*%Gm
    galat=abs(teta.hat-teta)
    iter=iter+1
  }
  
  cat('Nilai estimasi alpha =',teta.hat[1,1],'\n')
  cat('Nilai estimasi bera =',teta.hat[2,1],'\n')
  cat('Nilai diperoleh pada iterasi ke',iter,'\n')
}

#data
data=c(1.043,0.235,0.708,1.426,1.055,0.561,0.693,1.113,0.317,0.542,0.328,0.826,0.347,0.604,0.976,0.533,0.284,0.387,0.506,0.556,0.331,0.366,0.166,0.347,0.426)
mle.gamma(data,0.25,0.5) #konvergen
mle.gamma(data,2,4) #tidak konvergen

library(MASS)

#cek estimasi parameter
fitdistr (data, 'gamma')
1/6.456548 


mle.mult=function(fungsi,x,p,nP){
  tol=10^(-6)
  iter=1
  dfungsi=NULL
  for(i in 1:length(p)){
    dfungsi[[i]]=D(parse(text=fungsi),nP[i])
  }
  ddfungsi=NULL
  for(i in 1:length(p)){
    for(j in 1:length(p)){
      ddfungsi[[(i*(i-1)+j)]]=D(D(parse(text=fungsi),nP[i]),nP[j])
    }
  }
  for(i in 1:length(p)){
    assign(nP[i],p[i])
  }
  Hm=matrix(nrow=length(p),ncol=length(p))
  for(i in 1:length(p)){
    for(j in 1:length(p)){
      Hm[i,j]=sum(eval(ddfungsi[[(i*(i-1)+j)]]))
    }
  }
  Gm=matrix(nrow=length(p))
  for(i in 1:length(p)){
    Gm[i]=sum(eval(dfungsi[[i]]))
  }
  pHat=p-solve(Hm)%*%Gm
  galat=abs(pHat-p)
  while (any(galat>tol)) {
    for(i in 1:length(p)){
      assign(nP[i],pHat[i,1])
    }
    p=pHat
    Hm=matrix(nrow=length(p),ncol=length(p))
    for(i in 1:length(p)){
      for(j in 1:length(p)){
        Hm[i,j]=sum(eval(ddfungsi[[(i*(i-1)+j)]]))
      }
    }
    Gm=matrix(nrow=length(p))
    for(i in 1:length(p)){
      Gm[i]=sum(eval(dfungsi[[i]]))
    }
    pHat=p-solve(Hm)%*%Gm
    galat=abs(pHat-p)
    iter=iter+1
  }
  for(i in 1:length(p)){
    cat('Nilai estimasi',nP[i],'=',pHat[i,1],'\n')
  }
  cat('Nilai diperoleh pada iterasi ke',iter,'\n')
}

mle.mult("log(1/(gamma(alpha)*beta^alpha)*x^(alpha-1)*exp(-x/beta))",c(1.043,0.235,0.708,1.426,1.055,0.561,0.693,1.113,0.317,0.542,0.328,0.826,0.347,0.604,0.976,0.533,0.284,0.387,0.506,0.556,0.331,0.366,0.166,0.347,0.426),c(0.25,0.5),c("alpha","beta"))




simulate_queue <- function(arrival_rate, service_rate, num_servers, num_customers) {
  # Initialize vectors to store results
  Cust <- R1 <- interarrival <- arrival <- Serv_Start <- R2 <- Server <- Serv_Time <- Time_Serv_End <- numeric(num_customers)
  Cust_in_Queue <- Cust_in_System <- Time_in_Queue <- Time_in_System <- numeric(num_customers)
  Serv_End <- numeric(num_servers)
  
  for (i in 1:num_customers) {
    # Generate Cust
    Cust[i] <- i
    
    # Generate R1 and R2
    set.seed(i)
    R1[i] <- runif(1)
    set.seed(i + 3)
    R2[i] <- runif(1)
    
    # Generate Interarrival
    interarrival[i] <- -1 / arrival_rate * log(1 - R1[i])
    
    # Generate Server Time
    Serv_Time[i] <- -1 / service_rate * log(1 - R2[i])
    
    # Handle the rest of your logic here
    if (i == 1) {
      arrival[i] <- interarrival[i]
      Server[i] <- 1
      Serv_Start[i] <- arrival[i]
      Time_Serv_End[i] <- Serv_Start[i] + Serv_Time[i]
      Cust_in_Queue[i] <- 0
      Cust_in_System[i] <- 1
      Time_in_Queue[i] <- 0
    } else {
      arrival[i] <- arrival[i - 1] + interarrival[i]
      for (j in 1:num_servers) {
        if (arrival[i] > Serv_End[j]) {
          Server[i] <- j
          Serv_Start[i] <- max(arrival[i], Serv_End[j])
          Time_Serv_End[i] <- Serv_Start[i] + Serv_Time[i]
          Cust_in_Queue[i] <- max(0, i - num_servers)
          Cust_in_System[i] <- Cust_in_Queue[i] + num_servers
          Time_in_Queue[i] <- ifelse(Cust_in_Queue[i] > 0, Time_Serv_End[i - num_servers] - arrival[i], 0)
          break
        }
      }
      if (arrival[i] < min(Serv_End)) {
        k <- which.min(Serv_End)
        Server[i] <- k
        Serv_Start[i] <- arrival[i]
        Time_Serv_End[i] <- Serv_Start[i] + Serv_Time[i]
        Cust_in_Queue[i] <- Cust_in_Queue[i - 1] + 1
        Cust_in_System[i] <- Cust_in_System[i - 1] + 1
        Time_in_Queue[i] <- Serv_Start[i] - Time_Serv_End[i - 1]
      }
    }
    Time_in_System[i] <- Time_in_Queue[i] + Serv_Time[i]
    Serv_End[Server[i]] <- Time_Serv_End[i]
  }
  
  result <- data.frame(
    Cust = Cust,
    R1 = R1,
    Interarrival = interarrival,
    Arrival = arrival,
    Server = Server,
    Serv_Start = Serv_Start,
    R2 = R2,
    Serv_Time = Serv_Time,
    Serv_End = Serv_End,
    Cust_in_Queue = Cust_in_Queue,
    Cust_in_System = Cust_in_System,
    Time_in_Queue = Time_in_Queue,
    Time_in_System = Time_in_System
  )
  
  return(result)
}
#this is new line


# Call the function and store the result
result <- simulate_queue(arrival_rate = 30, service_rate = 10, num_servers = 3, num_customers = 30000)
result
