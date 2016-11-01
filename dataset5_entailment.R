

library(mirt)
library(psych)

answers_5 <- read.csv("irt_file_5_coded.csv",header=TRUE,sep=",")
Key<-unclass(read.table("group_5.txt",sep=" ",comment.char="")[,2]);

snli_5 <- answers_5[2:91]

c_s5 <- colnames(snli_5)

Data <- snli_5;

items<-which(Key==2);

######################
### Dimensionality ###
######################

R<-tetrachoric(Data[,items])
plot(eigen(R$rho,symmetric=TRUE,only.values=TRUE)$values)

### One single high eigenvalue before trailing eigenvalues, suggesting a single factor

fitMIRT1<-mirt(Data[,items],1,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e3))
fitMIRT2<-mirt(Data[,items],2,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e3))
anova(fitMIRT1,fitMIRT2) ## Sample size adjusted BIC suggests a single factor
summary(fitMIRT2,suppress=0.3) ## 2F model does show an interpretable loading pattern. Second factor has small SS loadings

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

## No lower asymptote is significant. Should change to 2PL.

fit.2PL<-mirt(Data[,items],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
anova(fit.2PL,fitMIRT1) ## Test confirm 2PL model is as good as 3PL.
summary(fit.2PL) ## There are two items with h2 (communality) smaller than 15%. They should be removed.
coef(fit.2PL,simplify=TRUE,IRTpar=TRUE)$items; ## The same two items with lowest diffuculty. 

## Some other items also have very low difficulty, but they are retained given the nature of the project.
## If the scale were intended for normal human subjects (instead of computers), they should be removed.

#############################
### Remove Items 1 and 22 ###
#############################

fit.2<-mirt(Data[,items[c(-1,-22)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itemfit(fit.2) ## No significant p value, showing the logistic form of each item is reasonable.
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 
res.table<-residuals(fit.2)
pstar<-28*27/2
cbind(sort(1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))),0.05*(1:pstar)/pstar)
## Compare the p-values to the Benjamini-Hochberg criteria. Nothing is significant.

