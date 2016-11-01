# R code for fitting an IRT model for the 5GS-Neutral items

# load libraries (assuming they are already installed)
library(mirt)
library(psych)

# load data and key
answers_5 <- read.csv("data/irt_file_5_coded.csv",header=TRUE,sep=",")
Key<-unclass(read.table("data/group_5.txt",sep=" ",comment.char="")[,2]);

snli_5 <- answers_5[2:91]
c_s5 <- colnames(snli_5)
Data <- snli_5;

items<-which(Key==3); #key 3 is neutral

######################
### Dimensionality ###
######################

R<-tetrachoric(Data[,items])
plot(eigen(R$rho,symmetric=TRUE,only.values=TRUE)$values)


### One single high eigenvalue before trailing eigenvalues, suggesting a single factor

# fit 1 factor and 2 factor models
fitMIRT1<-mirt(Data[,items],1,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e4))
fitMIRT2<-mirt(Data[,items],2,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e4))
anova(fitMIRT1,fitMIRT2) ## results here support a 2F model. 
summary(fitMIRT2,suppress=0.3) ## 2F model does show an interpretable loading pattern. Second factor has very small SS loadings

# new rotation to squeeze into single factor
target<-cbind(rep(1,30),rep(0,30));
summary(fitMIRT2,suppress=0.3,rotate="pstQ",Target=target,W=!target)


#######################
## Evaluate 1F model ##
#######################

fit.all.est<-coef(fitMIRT1,simplify=TRUE,IRTpar=TRUE)$items;
mtlnL.all<- -2*fitMIRT1@Fit$logLik;
pvalue3v2<-rep(NA,30);
type<-rep("3PL",30);
for (k in 1:30)
{
  if (fit.all.est[k,"g"]<1e-2) {pvalue3v2[k]<-.5; next;}   
  itemtype<-type;
  itemtype[k]<-"2PL";
  fit.2PL <-mirt(Data[,items],1,itemtype,TOL=2e-4,technical=list(NCYCLES=1e3))
  mtlnL.2PL<- -2*fit.2PL@Fit$logLik;
  T3v2<-mtlnL.2PL-mtlnL.all;
  pvalue3v2[k]<-p<-pchisq(T3v2,df=1,lower.tail=FALSE)/2;
  cat(items[k],": 3PL vs 2PL ",T3v2," p = ", p, "\n");
}
s<-sort.int(pvalue3v2,index.return=TRUE)
print(cbind(items[s$ix],s$x,0.05*(1:30)/30));

## No lower asymptote is significant. Should change to 2PL since none of the 
# guessing parameters are significant.

fit.2PL<-mirt(Data[,items],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
anova(fit.2PL,fitMIRT1) ## Test confirm 2PL model is as good as 3PL.
summary(fit.2PL) ## There are no items with h2 (communality) smaller than 15%. 
coef(fit.2PL,simplify=TRUE,IRTpar=TRUE)$items; ##  all difficulties are negative, but none look too difficult.

#############################
### Don't Remove any Items###
#############################

fit.2<-mirt(Data[,items],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2)  
itf ## 4 items with p < 0.05.
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 
res.table<-residuals(fit.2)
pstar<-30*29/2

#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(30:1))

#one significant item found. which one?
itf #item 7

#plot item 7
thetaF <- fscores(fit.2)
ig0 <- itemGAM(Data[,items[c(7)]],thetaF)
ig.good <- itemGAM(Data[,items[c(3)]],thetaF)
plot(ig0)
plot(ig.good)
itemplot(fit.2,7,main="Item Characteristic Curve",col.line="black")


#remove item 7 & refit

fit.2<-mirt(Data[,items[c(-7)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## no significant p values 
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 
res.table<-residuals(fit.2)

cbind(sort(itf$p.S_X2),0.05/(29:1)) #no significant values here

pstar<-29*28/2

t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1)) ## Compare the p-values using Holm's method 


## Refit without items 3 & 7

fit.2<-mirt(Data[,items[c(-3,-7)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) 
itf ## no significant p values 
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 
res.table<-residuals(fit.2)

cbind(sort(itf$p.S_X2),0.05/(28:1)) #no significant values here

itemplot(fit.2,10,main="Item Characteristic Curve",col.line="black")

pstar<-28*27/2

t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1)) ## Compare the p-values using Holm's method 

# no significant residuals, so we are done. Retain all items except items 3 & 7 
