

library(mirt)
library(psych)
require(graphics)


answers_5 <- read.csv("irt_file_5_coded.csv",header=TRUE,sep=",")
Key<-unclass(read.table("group_5.txt",sep=" ",comment.char="")[,2]);

snli_5 <- answers_5[2:91]

c_s5 <- colnames(snli_5)

Data <- snli_5;

items<-which(Key==1); #key 1 is contradiction

#need to remove the items that were marked via visual inspection

which(c_s5[items] == 'X5794490455.jpg.0r1c')
which(c_s5[items] ==  'X4814925710.jpg.0r1c')

#remove items 22 and 26 for first run
c_s5[items]

######################
### Dimensionality ###
######################

R<-tetrachoric(Data[,items[c(-22,-26)]])
plot(eigen(R$rho,symmetric=TRUE,only.values=TRUE)$values)

### One single high eigenvalue before trailing eigenvalues, suggesting a single factor


fitMIRT1<-mirt(Data[,items[c(-22,-26)]],1,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e4))
fitMIRT2<-mirt(Data[,items[c(-22,-26)]],2,itemtype="3PL",TOL=1e-4,technical=list(NCYCLES=1e4))
anova(fitMIRT1,fitMIRT2) ## Sample size adjusted BIC suggests 2F model
summary(fitMIRT2,suppress=0.3) ## 2F model shows split between two factors

# new rotation to sqeeze into single factor
target<-cbind(rep(1,28),rep(0,28));
summary(fitMIRT2,suppress=0.3,rotate="pstQ",Target=target,W=!target) #now we see small loading for 2nd factor

#This suggests single factor.


#######################
## Evaluate 1F model ##
#######################

fit.all.est<-coef(fitMIRT1,simplify=TRUE,IRTpar=TRUE)$items;
mtlnL.all<- -2*fitMIRT1@Fit$logLik;
pvalue3v2<-rep(NA,28);
type<-rep("3PL",28);
for (k in 1:28)
{
  if (fit.all.est[k,"g"]<1e-2) {pvalue3v2[k]<-.5; next;}   
  itemtype<-type;
  itemtype[k]<-"2PL";
  fit.2PL <-mirt(Data[,items[c(-22,-26)]],1,itemtype,TOL=2e-4,technical=list(NCYCLES=1e4))
  mtlnL.2PL<- -2*fit.2PL@Fit$logLik;
  T3v2<-mtlnL.2PL-mtlnL.all;
  pvalue3v2[k]<-p<-pchisq(T3v2,df=1,lower.tail=FALSE)/2;
  cat(items[k],": 3PL vs 2PL ",T3v2," p = ", p, "\n");
}
s<-sort.int(pvalue3v2,index.return=TRUE)
print(cbind(items[s$ix],s$x,0.05/(28:1))); #holm criteria

## No lower asymptote is significant. Should change to 2PL.

fit.2PL<-mirt(Data[,items[c(-22)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
anova(fit.2PL,fitMIRT1) ## Test confirm 2PL model is as good as 3PL.
summary(fit.2PL) ## There are three items with h2 (communality) smaller than 15%. They should be removed.
coef(fit.2PL,simplify=TRUE,IRTpar=TRUE)$items;



# 9, 10, 14 should be removed
#############################
### Remove Items 9, 10, and 14 ###
#############################

fit.2<-mirt(Data[,items[c(-9,-10,-14,-22,-26)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## some p values are significant, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 


#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(25:1)) #one value is significant, item 9 (original item 11)


res.table<-residuals(fit.2)
pstar<-25*24/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# several significant pairs here


#now remove item 11


fit.2<-mirt(Data[,items[c(-9,-10,-11,-14,-22,-26)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## some p values are significant, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 
summary(fit.2)
#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(24:1)) # one is significant, item 24 (originally item 30)


res.table<-residuals(fit.2)
pstar<-24*23/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# 3 significant pairs found. 


# next remove item 30

fit.2<-mirt(Data[,items[c(-9,-10,-11,-14,-22,-26,-30)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## some p values are significant, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 

#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(23:1)) # one is significant, item 22 (originally item 28)


res.table<-residuals(fit.2)
pstar<-23*22/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# 2 significant pairs found. 


# remove item 28

fit.2<-mirt(Data[,items[c(-9,-10,-11,-14,-22,-26,-28,-30)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## some p values are significant, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 

#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(22:1)) # none are significant

res.table<-residuals(fit.2)
pstar<-22*21/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# still have 2 significant pairs found. 

residuals(fit.2,df.p=TRUE,suppress=0.13)
residuals(fit.2)

# highest LD value between X2326730558.jpg.4r2c and X630476551.jpg.4r1c
# items 20 & 19
# which has lower communality? 
summary(fit.2) #they're about the same, so we'll try each model.

# remove item 19

fit.2<-mirt(Data[,items[c(-9,-10,-11,-14,-19,-22,-26,-28,-30)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## some p values are below 0.05, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 

#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(21:1)) # none are significant

res.table<-residuals(fit.2)
pstar<-21*20/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# still have 1 significant pair found. 

#now let's try without 20

# remove item 28


fit.2<-mirt(Data[,items[c(-9,-10,-11,-14,-20,-22,-26,-28,-30)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## some p values are less than 0.05, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 

#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(21:1)) # none are significant

res.table<-residuals(fit.2)
pstar<-21*20/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# still have 1 significant pairs found. 

#since removing each left one significant value, remove both

fit.2<-mirt(Data[,items[c(-9,-10,-11,-14,-19,-20,-22,-26,-28,-30)]],1,itemtype="2PL",TOL=2e-4,technical=list(NCYCLES=1e3))
itf <- itemfit(fit.2) ## one p value less than 0.05, we'll need to test them with Holm
M2(fit.2) ## RMSEA_5 smaller than 0.5; TLI and CFI greater than 0.95. Good overall fit for 2-margins. 

#first test p values of itemfit using Holm criteria
cbind(sort(itf$p.S_X2),0.05/(20:1)) # none are significant

# uncomment to check original item order
# c_s5[items]

res.table<-residuals(fit.2)
pstar<-20*19/2
t <- 1-pchisq(df=1,abs(res.table[lower.tri(res.table)]))
cbind(sort(round(t,6)),0.05/(pstar:1))
## Compare the p-values to the Holm criteria. 
# no significant pairs found. 

# set of good items from this analysis:

# 1,2,3,4,5,6,7,8,12,13,15,16,17,18,21,23,24,25,27,29



