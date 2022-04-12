

#===============================================================================
#Simulation: independent hypothesis
#===============================================================================

sfdr <- function(m=20){
  p.true <- 1/4
  mu <- rep(0,m)
  mu[1:round(p.true*m)] <- 2.8
  y <- rnorm(rep(1,m),mu)
  y
  p.val <- 2*(1-pnorm(abs(y)))
  v.fwer <- p.val<0.05/m
  v.fdr  <- p.val[order(p.val)] < 1:m/m*.05
  pw.fwer <- mean(v.fwer[mu>0])
  pw.fdr  <- mean(v.fdr[mu[order(p.val)]>0])
  cbind(pw.fwer,pw.fdr) 
 }

#number of hypothesis 
numh <- sapply(1:10,function(x) 2^x)
#number of simulations
B <- 10^5


sres <- replicate(B,sapply(numh , sfdr))
res <- apply(sres,1:2,mean)


par(mfrow=c(1,2))
plot(c(1,max(numh)),c(0,1),type='n',axes=F ,xlim=c(0,64),   
      xlab = 'number of hypothesis', ylab = 'power')
axis(1,at=numh)
axis(2)
points(numh ,res[1,],typ='l',lwd=2,lty=2)
points(numh ,res[2,],typ='l',lwd=2)

plot(c(1,max(numh)),c(0,1),type='n',axes=F , 
      xlab = 'number of hypothesis', ylab = 'power')
axis(1,at=numh)
axis(2)
points(numh ,res[1,],typ='l',lwd=2,lty=2)
points(numh ,res[2,],typ='l',lwd=2)


#===============================================================================
#Simulation II: correlated hypothesis
#===============================================================================
library(MASS)
library(parallel)

sfdr2 <- function(m=20,rho){
  p.true <- 1/4
  mu <- rep(0,m)
  mu[1:round(p.true*m)] <- 2.8
  Sigma <- matrix(rho,m,m)
  diag(Sigma) <- 1
  y <-  mvrnorm(1, mu, Sigma)
  p.val <- 2*(1-pnorm(abs(y)))
  v.fwer <- p.val<0.05/m
  v.fdr  <- p.val[order(p.val)] < 1:m/m*.05
  pw.fwer <- mean(v.fwer[mu>0])
  pw.fdr  <- mean(v.fdr[mu[order(p.val)]>0])
  cbind(pw.fwer,pw.fdr) 
 }
#number of hypothesis  
numh <- sapply(1:8,function(x) 2^x)
#number of simulations
B <- 10^4 #10^5

    

cl <- makeCluster(8)
clusterExport(cl,c('sfdr2','numh') )
clusterEvalQ(cl,library('MASS'))
res.0 <- parSapply(cl, 1:B,  function(r)  sapply(numh , sfdr2,rho=0) )
res.0 <- matrix(apply(res.0,1,mean),2,length(numh))
res.5 <- parSapply(cl, 1:B,  function(r)  sapply(numh , sfdr2,rho=.5) )
res.5 <- matrix(apply(res.5,1,mean),2,length(numh))
res.8 <- parSapply(cl, 1:B,  function(r)  sapply(numh , sfdr2,rho=.8) )
res.8 <- matrix(apply(res.8,1,mean),2,length(numh))
stopCluster(cl)


plot(c(1,max(numh)),c(0,1),type='n',axes=F , 
      xlab = 'number of hypothesis', ylab = 'power')
axis(1,at=numh)
axis(2)
points(numh ,res.0[1,],typ='l',lwd=2,lty=2,col=adjustcolor( "black",1))
points(numh ,res.0[2,],typ='l',lwd=2,lty=1,col=adjustcolor( "black",1))
points(numh ,res.5[1,],typ='l',lwd=2,lty=2,col=adjustcolor( "black",.6))
points(numh ,res.5[2,],typ='l',lwd=2,lty=1,col=adjustcolor( "black",.6))
points(numh ,res.8[1,],typ='l',lwd=2,lty=2,col=adjustcolor( "black",.2))
points(numh ,res.8[2,],typ='l',lwd=2,lty=1,col=adjustcolor( "black",.2))

#===============================================================================
#DINGO Example
#===============================================================================
#from Kate file
#dpath <- 'C:/LOCAL/Y2/S/DINGO -20220311'
load(file.path(dpath,'dingoResults.rda')) 

library(dplyr)
library(huge)
library(iDINGO)
library(igraph)
library(matrixcalc)
library(rmdformats)


modBoot$diff.score[1:5]
modBoot$p.val[1:5]

dingoDF = data.frame("gene1"=modBoot$genepair$gene1,
                     "gene2"=modBoot$genepair$gene2,
                     "R1"=modBoot$R1,
                     "R2"=modBoot$R2,
                     "diff.value"=modBoot$R1 - modBoot$R2,
                     "diff.score"=modBoot$diff.score,
                     "p.val"=modBoot$p.val)
head(dingoDF)

hist(dingoDF$p.val,breaks=100)
sum(dingoDF$p.val < 0.05)


#Plot 1: Hypothesis rejected
pcer <- sum(dingoDF$p.val < 0.05) #per comparison error rate
fwer <- sum(dingoDF$p.val < 0.05/1275)
fdr  <- sum(dingoDF$p.val[order(dingoDF$p.val )]<1:1275/1275*.05)

plot(c(1,1275),c(0,1),type='n',axes=F ,ylim=c(1/10^13,1),log='y',xlim=c(0,200),
      xlab = 'hypothesis', ylab = 'p-value')
axis(1 )
axis(2, at=c(1/10^13,.05,1),labels=c(0,.05,1))
points(1:1275,dingoDF$p.val[order(dingoDF$p.val )],typ='l',lwd=2)

segments(0,.05,pcer,.05, col='red',lwd=2)
segments(pcer,1/10^14,pcer,.05,lty=2,col='red',lwd=2)
 
segments(fdr,1/10^14,fdr,dingoDF$p.val[order(dingoDF$p.val )][fdr],lty=2,col='purple',lwd=2)
points(1:fdr,1:fdr/1275*.05,typ='l',col='purple',lwd=2)

segments(0,0.05/1275,fwer,0.05/1275, col='orange',lwd=2)
segments(fwer,1/10^14,fwer,0.05/1275,lty=2,col='orange',lwd=2) 


#Plot 2: Resulting networks
fdr  <-   dingoDF$p.val[order(dingoDF$p.val )]<1:1275/1275*.05
fwer <-   dingoDF$p.val < 0.05/1275
    
edgeList = dingoDF %>% 
  filter(fwer) %>% 
  select(gene1, gene2, diff.value) %>%
  rename(weight = diff.value)


myGraph = graph_from_edgelist(as.matrix(edgeList[,1:2]),directed=F)
E(myGraph)$weight = edgeList$weight
#plot(myGraph,      layout = layout_with_dh)
     
plot(myGraph, 
     layout = layout_with_lgl,
     vertex.size = degree(myGraph),
     edge.width = 20*abs(E(myGraph)$weight), # 20 is arbitrary
     edge.color = ifelse(E(myGraph)$weight > 0, "red","blue"))