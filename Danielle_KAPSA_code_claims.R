library(evir)
library(quantmod)
library(extraDistr)
library(rmutil)
#library(actuar)

# download data
PRC=read.csv2("/home/aimsstudent21/Documents/ESSAY AIMS/good/PRCBON.csv")
# chosen the data with total affected records greater than 0
PRC3=PRC[PRC$Total.Records>0,]
# extract the variable total affected records
n=PRC2$Total.Records
#scale the variable by 100 000
n=n/10000
PRC2=PRC[PRC$Total.Records>277999,]
Total.Records = PRC3$Total.Records/10000

library(e1071) 
kurtosis(n)
skewness(n)   

# sort the total affected records
Total.Records <- sort(n, decreasing = TRUE)
#histogramm of the total records 
hist(Total.Records, col="blue", ylim = c(0,5))

#Total.Records1=PRC2$Total.Records
#Total.Records1=Total.Records/100000
#hist(Total.Records1, 20, col = 4)
#sd(Total.Records1)
#summary(Total.Records1)

qplot(Total.Records,xi=0, xlim = c(0,5), ylim = c(0, 1500))
emplot(Total.Records, 'xy', col = "blue")

library(ineq)
sort_Total.Records=sort(Total.Records1) # We sort the data
n=length(Total.Records1)
CP=c() #Empty array for storage
for (i in 1:n) {
    CP[i]=ineq(sort_Total.Records[i:n],type="Gini") # Truncated Gini
  }
plot(1:n,CP,ylim=c(0,1)) 


MSplot <- function(data,p=4) {
  par(mfrow = c(2, 2)) 
  x=abs(data)
  for (i in 1:p) {
    y=x^i
    S=cumsum(y)
    M=cummax(y)
    R=M/S
    plot(1:length(x),R,type='l', col=2, lwd=3, ylim=c(0,1),xlab='n', ylab='Rn', 
        main=paste("MSplot for p=",i))
  }
  par(mfrow = c(1, 1)) 
   return(R)
}

MSplot(Total.Records)
#Total.Records2= Total.Records/100000
library(extRemes)
fit1= fevd(Total.Records2, type = "GEV", method = "MLE", units = "m3/s" )

fit1######################################################


matricecroise <- function(seuill, data){
  tail.seuil <- length(seuill)
  tail.data <- length(data)
  MatTotal.Records<-matrix(rep(data,tail.seuil),nrow=tail.seuil,ncol=tail.data,byrow=TRUE)
  Matseuil<-matrix(rep(seuill,tail.data),nrow=tail.seuil,ncol=tail.data ,byrow = FALSE)
  Matexces <- MatTotal.Records-Matseuil
  # replacing the negative value by 0
  Matexces[Matexces<0] <- 0
  Matexces 
}
seuil<-seq(0, 300000, 50)
matrice <- matricecroise(seuil, Total.Records)
matrice[1:8, 1:8]


vecteur <- function(seuill, data){
  tail.seuil <- length(seuill)
  tail.data <- length(data)
  matrice <- matricecroise(seuill, data)
  compteur <- apply((matrice>0),1,sum) # On compte le nombre d’ exces positifs
  e <- vector("numeric",tail.seuil)
  for (i in 1:tail.seuil)
    if (compteur[i] !=0) {e[i]<- sum(matrice[i,1:tail.data])/compteur[i]}
  e
}
e <- vecteur(seuil, Total.Records)
e[1:35]

plot(seuil, e, type = "l", xlab = "threshold", ylab = "Means excesses" ,col = "blue")

k <- 10
xi <-(1/k)*sum(log((Total.Records[length(Total.Records)-k+1])/(Total.Records[length(Total.Records)-k])))
xi

plot(seuil, e, type = "l", xlab = "threshold", ylab = "average excesses ")


#Applications aux lois usuelles 
#Weibull

seuilord <- sort (Total.Records)
e.ord <- vecteur(seuilord, Total.Records)

N <- 2000
seuil1 <- sort(rnorm(N))
seuil2 <- sort(rweibull(N, 8, 1))
t.weibull <- sort(rweibull(N, 8, 1), decreasing = TRUE)
matseuil2 <- matricecroise(seuil, t.weibull)
e2 <- vecteur(seuil2, t.weibull)
split.screen(c(1,2))
screen(1)
plot(seuilord, e.ord, type = "l", main = "Totals Records", xlab = "threshold", ylab = "average excesses", col = "blue")
screen(2)
plot(seuil2, e2, type = "l", main = " Weibull Distribution", xlab = "threshold", ylab = "average excesses", col = "blue")

 # Loi de Gumbel
seuil6 <- sort(rgumbel(N, 0, 1), decreasing = FALSE)
t.gumbel <- sort(rgumbel(N, 0, 1), decreasing = TRUE)
matseuil6 <- matricecroise(seuil, t.gumbel)
e6 <- vecteur(seuil6, t.gumbel)
split.screen(c(1,2))
screen(1)
plot(seuilord, e.ord, type = "l", main = "Totals Records", xlab = "threshold", ylab = "average excesses", col ="blue")
screen(2)
plot(seuil6, e6, type = "l", main = " Gumbel distribution", xlab = "threshold", ylab = "average excesses", col ="blue")


# Loi de frechet


seuil7 <- sort(rfrechet(N,1, 0, 1), decreasing = FALSE)
t.frechet <- sort(rfrechet(N,1, 0, 1), decreasing = TRUE)
matseuil7 <- matricecroise(seuil, t.frechet)
e7 <- vecteur(seuil7, t.frechet)
split.screen(c(1,1))
screen(1)

plot(seuilord, e.ord, type = "l", main = "Totals Records", xlab = "threshold", ylab = "average excesses", col = "blue")

screen(2)
plot(seuil7, e7, type = "l", main = " Frechet distribution", xlab = "threshold", ylab = "average excesses", col ="blue")


#Function empirique des exces
#split.screen(c(3,3))
#screen(1)
#meplot(Total.Records, 2, main = "Total.Records", type = "l")
#screen(2)
#meplot(rnorm(N), 3, main ="Normal" , type = "l")
#screen(3)
#meplot(rweibull(N, 1, 1), 3, main = "Weibull", type = "l")
#screen(4)
#meplot(rpareto(N, 2, 2), 3, main ="Pareto" , type = "l")
#screen(5)
#meplot(rgumbel(N, 0, 1), 3, main = "Gumbel", type = "l")
#screen(6)
#meplot(rfrechet(N, 1, 0, 1), 3, main = "Frechet", type = "l")

#meplot(rcauchy(N), 3, main ="Cauchy" , type = "l")


egevd(Total.Records, ci = TRUE, conf.level = 0.9)


#############
fit_mle <- fevd(as.vector(Total.Records), method = "MLE", type="GEV")
fit_mle
plot(fit_mle)
fit_lmom <- fevd(as.vector(Total.Records), method = "Lmoments", type="GEV")
fit_lmom
plot(fit_lmom)


# return levels:
#rl_mle <- return.level(fit_mle, conf = 0.05, return.period= c(2,5,10,20,50,100))

# return levels:
#rl_lmom <- return.level(fit_lmom, conf = 0.05, return.period= c(2,5,10,20,50,100))


# return level plots
#par(mfcol=c(1,2))
# return level plot w/ MLE
#plot(fit_mle, type="rl",
 #    main="Return Level Plot for Total Records w/ MLE",
  #   ylim=c(0,200), pch=16)
#loc <- as.numeric(return.level(fit_mle, conf = 0.05,return.period=100))
#segments(100, 0, 100, loc, col= 'midnightblue',lty=6)
#segments(0.01,loc,100, loc, col='midnightblue', lty=6)

# return level plot w/ LMOM
#plot(fit_lmom, type="rl",
 #    main="Return Level Plot for Total Records w/ L-Moments",
  #   ylim=c(0,200))
#loc <- as.numeric(return.level(fit_lmom, conf = 0.05,return.period=100))
#segments(100, 0, 100, loc, col= 'midnightblue',lty=6)
#segments(0.01,loc,100, loc, col='midnightblue', lty=6)

# comparison of return levels
#results <- t(data.frame(mle=as.numeric(rl_mle),
 #                       lmom=as.numeric(rl_lmom)))
#colnames(results) <- c(2,5,10,20,50,100)
#round(results,1)

