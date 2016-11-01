
library(mirt)
library(psych)

answers_4 <- read.csv("irt_file_4_coded.csv",header=TRUE,sep=",")
Key<-unclass(read.table("group_4.txt",sep=" ",comment.char="")[,2]);

snli_4 <- answers_4[2:91]

c_s4 <- colnames(snli_4)

Data <- snli_4;

items<-which(Key==3); #key 3 is neutral


######################
### Dimensionality ###
######################

R<-tetrachoric(Data[,items])
plot(eigen(R$rho,symmetric=TRUE,only.values=TRUE)$values)

### One single high eigenvalue before trailing eigenvalues, suggesting a single factor
#this is still true for contradiction

fitMIRT1<-mirt(Data[,items],1,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e4))
fitMIRT2<-mirt(Data[,items],2,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e4))
anova(fitMIRT1,fitMIRT2) ## Sample size adjusted BIC suggests 2F
summary(fitMIRT2,suppress=0.3) ## 2F model does show an interpretable loading pattern. Second factor has small SS loadings


# new rotation to sqeeze into single factor
target<-cbind(rep(1,30),rep(0,30));
summary(fitMIRT2,suppress=0.3,rotate="pstQ",Target=target,W=!target)

#new summary suggests single factor

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
summary(fit.2PL) ## There are no items with h2 (communality) smaller than 15%. They should be removed.
coef(fit.2PL,simplify=TRUE,IRTpar=TRUE)$items; 

#############################
### Don't Remove any Items###
#############################

fit.2<-mirt(Data[,items],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## One significant p value, need to check using Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins.

summary(fit.2) # no items have communality less than 0.15
cbind(sort(itf$p.S_X2),0.05/(30:1)) #none significant


res.table<-residuals(fit.2)
pstar<-30*29/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))

cbind(sort(round(t,8)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. Two significant pairs.
## Let's remove item 11, which had low itemfit p-value

# remove item 23

fit.2<-mirt(Data[,items[c(-23)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## One significant p value, need to check using Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins.

summary(fit.2) # no items have communality less than 0.15
cbind(sort(itf$p.S_X2),0.05/(29:1)) #none significant


res.table<-residuals(fit.2)
pstar<-29*28/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))

cbind(sort(round(t,8)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. Two significant pairs.


# remove item 17

fit.2<-mirt(Data[,items[c(-17)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## One significant p value, need to check using Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins.

summary(fit.2) # no items have communality less than 0.15
cbind(sort(itf$p.S_X2),0.05/(29:1)) #none significant


res.table<-residuals(fit.2)
pstar<-29*28/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))

cbind(sort(round(t,8)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. Two significant pairs.

#remove both 17 & 23

fit.2<-mirt(Data[,items[c(-17,-23)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## One significant p value, need to check using Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins.

summary(fit.2) # no items have communality less than 0.15
cbind(sort(itf$p.S_X2),0.05/(28:1)) #none significant

res.table<-residuals(fit.2)
pstar<-28*27/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))

cbind(sort(round(t,8)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. Two significant pairs.
## no more significant residuals


# good items from here, all except 17 & 23
# 1:16,18:22,24:30



