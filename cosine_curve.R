oldpar <- par(no.readonly = TRUE)

#get library
library(readxl)
library(dplyr)
library(ggplot2)

#input sample fpkm
FPKM <- read_excel("FPKM_annotation.xlsx")　
FPKM <- as.data.frame(FPKM)
FPKM2 <- FPKM[,3:52] 
FPKM2 <- as.data.frame(FPKM2) #make dataframe
rownames(FPKM2) <- FPKM$gene_id #change rowname to $gene_id of FPKM
rownames(FPKM) <- FPKM$gene_id

#input sample atribute
at <- read_excel("sampleatribute.xlsx") 
at <- as.data.frame(at) 

s <- sprintf("%02d/%02d/%02d",at$month,at$day,at$year)  #set time
s.time <- strptime(s,"%m/%d/%y") 
at <- cbind(at,s.time) 

#sort FPKM2 by time-series order
soat <- order(at$s.time)
at <- at[soat,]
FPKM2 <- FPKM2[,at$name]

l2FPKM2 <- log2(FPKM2+0.001)

tmp <- apply(l2FPKM2,1,mean)
expressed_gene <- tmp > -6
l2FPKM2.e <- l2FPKM2[expressed_gene,]

###Non-linear regression analysis 
###set time
xtime <- at$s.time

xtime3 <- format(xtime, "%j") %>% as.numeric
xtime3[at$year==15] <- xtime3[at$year==15]
xtime3[at$year==16] <- xtime3[at$year==16]+365
xtime3[at$year==17] <- xtime3[at$year==17]+365*2

###fit cosine curve
###for all genes
para.sub <- matrix(NA, nrow=nrow(l2FPKM2.e), ncol=3)

for(n in 1:nrow(l2FPKM2.e)){
  
  y <- as.numeric(l2FPKM2.e[n, ])
  C <- mean(y)
  
  fit.nls <- nls(y~ C+(1/2)*alpha*cos(2*pi*(xtime3-phi)/365),
                 start = list(alpha = 0.01, phi = 0),
                 control=nls.control(warnOnly=TRUE))
  
  para <- coef(fit.nls)
  
  if(para["alpha"] < 0){
    para["alpha"] <-  para["alpha"] * (-1)
    para["phi"] <- para["phi"]+365/2
  }
  
  para["phi"] <- para["phi"] %% 365
  
  para.sub[n,1:2] <- para
  para.sub[n,3] <- 1-deviance(fit.nls)/sum((y-C)^2)
  
}

colnames(para.sub) <- c("alpha","phi","RS")
rownames(para.sub) <- rownames(l2FPKM2.e)
para.sub <- as.data.frame(para.sub)

para.sub2 <- para.sub[,1:2] #dataframe only alpha and phi

hist(para.sub$phi, main = "Histogram of phase",
     breaks=seq(0,365,5),xlab = "Phase(day)",ylab = "Gene")
hist(para.sub$alpha,breaks = seq(0,15,0.1),main = "Histogram of Amplitude",
     xlab = "Amplitude",ylab="Gene")


###filtering by RS
hist(para.sub$RS, breaks=seq(-0.1,1,0.01),main = "The coefficient of determination",
     xlab="",ylab="Gene")
abline(v=0.2,col="red")

tmp <- para.sub$RS
sum(tmp > 0.1) #14295
sum(tmp > 0.2) #10419
sum(tmp > 0.3) #7837

q <- (para.sub$RS > 0.2)
cos.sub <- para.sub[q,]

###check cosine fitting
hist(cos.sub$alpha,breaks = seq(0,15,0.1),main = "Histogram of Amplitude",
     xlab = "Amplitude",ylab="Gene")
abline(v=1,col="red")

hist(cos.sub$phi,main = "Histogram of phase",
     breaks = seq(0,365,5),
     xlab = "Phase(day)",ylab = "Gene")

hist(cos.sub$RS, breaks = seq(0,1,0.01), xlab = "",
     main = "The coefficient of determination", ylab="Gene")

#extract bsSOGs
p <- (cos.sub$alpha >= 1)
SO.sub <- cos.sub[p,] #3341genes 

#check check parameters of bsSOGs
hist(SO.sub$alpha,breaks = seq(0,15,0.1),
     xlab = "Amplitude",main="Histogram of Amplitude",ylab="Gene")

hist(SO.sub$phi,main = "Histogram of phase",
     breaks=seq(0,365,5),xlab = "Phase(day)",ylab = "Gene")
