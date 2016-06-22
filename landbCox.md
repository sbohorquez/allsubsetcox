All Subset Regression Cox
================
Santiago Bohórquez
June 22, 2016

All subsets variable selection for Cox Model
--------------------------------------------

*stream of consciousness for now*

First try to program as if it was for a package:

First, set up the environment

``` r
setwd("D:/Desktop/RM/Programs")
rm(list=ls(all=TRUE))

memory.limit(size=20000)

library(survival)
library(MASS)
```

What I want to do is estimate Mallows \(C_p\), using Kuk (1984) procedure \[*need citation*\] \(\Lambda^{'}_{\alpha}=W^{'}_{\alpha}+2p_{\alpha}\) for all subset regressions, where \(W^{'}_{\alpha}=\beta^{'}_{2}C^{-1}_{22}\beta^{'}_{2}\) and \(C=I^{-1}\), and \(\beta^{'}_{2}\) are the restricted coefficients. And \(C_{p_{\alpha}}=\Lambda^{'}_{\alpha} - (p-1)\). The function `Cp` does this estimation

``` r
Cp<-function(i,est,size,data,included,sl){
  rem<-sapply(1:sl[1],saca,data=data,est=est,included=included,i=i)
  rem<-unlist(rem)
  Wprime<-est$coefficients[-rem]%*%ginv(est$var[-rem,-rem])%*%est$coefficients[-rem]
  Lambdaprime<-Wprime+2*length(est$coef[rem])
  Cp<-Lambdaprime-(length(est$coef)-1)
  c(Cp,included[,i])
  }
```

In the previous function we take into account that for factor variables we take out or include all the factors of a given variable, not one by one. The `saca` function checks all coefficients generated with a name from the included variables and includes all of them at the same time.

``` r
saca<-function(l,data,est,included,i){
   if (nlevels(data[,included[l,i]])==0) {included[l,i]
   }
   else {
     factores<-paste(names(data)[included[l,i]],levels(data[,included[l,i]]),
                     sep="")
     which(names(est$coef) %in% factores )
   }
 }
```

``` r
combinations<-function(k,size,est,data,nbest){
  included<-combn(1:size[2],k)
  sl<-dim(included)
  MCp<-sapply(1:sl[2],Cp,est=est,size=size,data=data,
                      included=included,sl=sl)
  MCp<-matrix(MCp,ncol=sl[2])
  best<-sort(MCp[1,],index.return=TRUE)
  cho<-min(nbest,size[2])
  MCp[,best$ix[1:cho]]
}
```

``` r
leapscox<-function(t,...) UseMethod("leapscox")

leapscox.default<-function(Survobj,data,nbest=5,...){
  est<-coxph(Survobj~.,data=data)
  size<-dim(data)
  nombres<-names(data)
  cpalpha<-sapply(1:(size[2]-1),combinations,size=size,est=est,data=data,
                  nbest=nbest) 
  cpval<-lapply(cpalpha, `[`,1,)
  bcp<- sort(unlist(cpval))[1:nbest]
  chosen<-sapply(1:(size[2]-1), function(x) which( cpval[[x]] %in% bcp ))
  suu<-sapply(1:(size[2]-1),function(x) cpalpha[[x]][,chosen[[x]]]) 
  suu[lapply(suu,length)>0] 
}

### falta print, summary y print.summary
### Correr modelos escogidos en la funcion?
```

``` r
#require(emplik)

Myeloma<-read.csv2("Krall1975.csv",header=T)
X<-Myeloma[,grep("X",names(Myeloma))]
X<-subset(X,select=-X0)
attach(Myeloma)

aja1<-leapscox(Surv(t,A.D),data=X)
aja1
```

    ## [[1]]
    ## [1]  5.532686  1.000000  3.000000  4.000000  6.000000  7.000000  8.000000
    ## [8] 12.000000 13.000000
    ## 
    ## [[2]]
    ##            [,1]     [,2]      [,3]
    ##  [1,]  6.346075  6.59486  6.944637
    ##  [2,]  1.000000  1.00000  1.000000
    ##  [3,]  3.000000  2.00000  2.000000
    ##  [4,]  4.000000  3.00000  4.000000
    ##  [5,]  6.000000  4.00000  6.000000
    ##  [6,]  7.000000  6.00000  7.000000
    ##  [7,]  8.000000  7.00000  8.000000
    ##  [8,] 12.000000  8.00000 12.000000
    ##  [9,] 13.000000 12.00000 13.000000
    ## [10,] 14.000000 13.00000 14.000000
    ## 
    ## [[3]]
    ##  [1]  6.61947  1.00000  2.00000  3.00000  4.00000  6.00000  7.00000
    ##  [8]  8.00000 12.00000 13.00000 14.00000

``` r
# names(x1)[c(3.00000,4.00000,5.00000,6.00000,  7.00000,  8.00000, 10.00000, 
#            11.00000)]
# 
# bb<-coxph(Surv(t,A.D)~.,data=X)
# bb2<-coxph(Surv(t,A.D)~X1+X3+X4+X6+X7+X8+X12+X13,data=X)
# 
# 2*(bb$loglik[2]-bb2$loglik[2])+2*8-15
```

Results differ from Kuk 1984, I have two hypothesis for this, either, due to advances in software, i.e., loglik for models is different. Or, I messed up inputing the data from Krall, et al 1975.

*Extenssions missing: y~. try all y~x always have x*
