library('e1071')
Fmeasure <- function(x){
  pr <- x[1,1]/sum(x[1,])
  rc <- x[1,1]/sum(x[,1])
  return(2*pr*rc/(pr+rc))
}
Gmean <- function(x){
  TPR <- x[1,1]/sum(x[,1])
  TNR <- x[2,2]/sum(x[,2])
  return(sqrt(TPR*TNR))
}
AUC <- function(x){
  FPR <- x[1,2]/sum(x[,2])
  TPR <- x[1,1]/sum(x[,1])
  
  return((1+(1-FPR))*TPR/2+(1-FPR)*(1-TPR)/2)
}
STRPart <- function(dt, r, p){  
  dt <- as.data.frame(dt)
  dt[,r] <- as.factor(dt[,r]) 
  n.class <- length(levels(dt[,r])) 
  dt.lst <- vector("list",n.class) 
  ndt.lst <- vector("list",n.class) 
  dt.trn.lst <-vector("list",n.class)  
  dt.tst.lst <- vector("list",n.class)  
  
  i <- 1     
  while(i <= n.class){   
    dt.lst[[i]] <- dt[dt[,r]==levels(dt[,r])[i],] 
    ndt.lst[[i]] <- nrow(dt.lst[[i]])             
    dt.lst[[i]] <- dt.lst[[i]][sample(ndt.lst[[i]]),]  
    dt.trn.lst[[i]] <- dt.lst[[i]][1:round(ndt.lst[[i]]*p),] 
    dt.tst.lst[[i]] <- dt.lst[[i]][(round(ndt.lst[[i]]*p)+1):ndt.lst[[i]],]    
    i <- i+1 
  }    
  j <- 2 
  d.train <- dt.trn.lst[[1]] 
  d.test <- dt.tst.lst[[1]]     
  while(j <= n.class){ 
    b <- dt.trn.lst[[j]]    
    d.train <- rbind(d.train,b)  
    c <- dt.tst.lst[[j]]     
    d.test <- rbind(d.test,c)
    j <- j+1                 
  }  
  return(list(d.train,d.test)) 
}
####
ACO <- function(training,testing,ITA=50,ant_n=50,rho=0.8){
  
  # sigma <- 5 #RBF kernal parameter
  # svm(x = tmp_trn[,-3],y = tmp_trn[,3],kernel = exp(-(abs(u-v)^2)/2*sigma^2),cost = 500)
  S_major <- training[training[,3]==0,]
  S_minor <- training[training[,3]==1,]
  pheromone <- matrix(1,nrow(S_major),2)
  OPS <- 0
  for(i in 1:ITA){
    fitness_old <- 0
    for(j in 1:ant_n){
      sel <-runif(nrow(S_major)) > pheromone[,1]/rowSums(pheromone)
      SS <- S_major[sel,]
      tmp_trn <- rbind(SS,S_minor)
      classifier <- svm(x = tmp_trn[,-3],y = tmp_trn[,3],cost = 500)
      mat <- table(predict(classifier,testing[,-3]),testing[,3])
      fitness_new <- 1/3*Fmeasure(mat)+1/3*Gmean(mat)+1/3*AUC(mat)
      
      if(fitness_new > fitness_old){
        best_sel <- sel
        best_set <- SS
        best_fitness <- fitness_new
        fitness_old <- fitness_new
      }
    }
    if(best_fitness>OPS){
      OPS <- best_fitness
      result_set <- best_set
      result_sel <- best_sel
    }
    
    delta <- matrix(0,nrow(S_major),2)
    delta[best_sel,2] <-1
    delta[!best_sel,1] <-1
    delta[delta==1] <- 1/0.1/ant_n*best_fitness
    pheromone <- rho*pheromone+delta
    print(paste('iteration',as.character(i),'finished.'))
  }
  output <- list()
  output$set <- rbind(result_set,S_minor)
  output$sel <- result_sel
  return(output)
}
####

data <- read.csv("dataset3.csv",h=T)
plot(data[,-3])
points(data[data[,3]==1,-3],col=2)

part <- STRPart(data,3,0.66)
training <- part[[1]]
test <- part[[2]]


maj <- training[training[,3]==0,]
min <- training[training[,3]==1,]
REC  <- vector("numeric",nrow(maj))
names(REC) <- rownames(maj)

for(rep in 1:100){
select <- STRPart(training,3,0.66)
S <- select[[1]]
V <- select[[2]]
name <-rownames(S[S[,3]==0,])
tmp <- ACO(training = S,testing = V)$sel
REC[name[tmp]] <- REC[name[tmp]]+1
print(rep)
}

under_result <- rbind(maj[names(REC[order(REC,decreasing = T)][1:30]),],min)

plot(training[,-3],col=training[,3])
plot(under_result[,-3],col=under_result[,3])


apply(maj[1:2],1,mean)
apply(min[1:2],1,mean)
