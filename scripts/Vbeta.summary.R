#module load vcftools before running this script, if reference LD not yet calculated

# which functions did i actually use from these package
# let's start with none of them
#library(foreach,lib="/restricted/projectnb/vf-thesis/")
#library(glmnet,lib="/restricted/projectnb/vf-thesis/")
#library(lars,lib="/restricted/projectnb/vf-thesis/")
#library(glmpath,lib="/restricted/projectnb/vf-thesis/")

snp<-commandArgs(trailingOnly=T)

#summary statistics in the locus
locsum<-read.table(paste(snp,".loc.sumstats.txt",sep=""),header=T,as.is=T)

#extract subset of nominally significant variants with positive prior variance of effect
nomsig<-which(locsum$p < .05 & !is.na(locsum$Freq1.Hapmap) & locsum$Vbeta>0)
aset.sum<-locsum[nomsig,]

aset<-which(locsum$BP %in% aset.sum$BP)


## this should be a helper function
#creates LD matrix for variants in the analysis
ld.long<-read.table(paste(snp,".1000G.loc.hap.ld", sep=""), header=T, as.is=T)
LDmat<-reshape(ld.long, direction="wide",
               v.names=c("R.2"),idvar=c("POS1"),
               timevar=c("POS2"),
               drop=c("CHR","D","Dprime"))
LDmat<-rbind(cbind(NA,LDmat[,-(1:2)]),NA)
LDmat[lower.tri(LDmat)]<-t(LDmat)[lower.tri(LDmat)]
diag(LDmat)<-1

pos<-unique(c(ld.long$POS1,ld.long$POS2))
rownames(LDmat)<-pos
colnames(LDmat)<-pos

LDmat<-as.matrix(LDmat[aset,aset])
## end helper function

#these should be arguments to the main function
n<-aset.sum$N
AFs<-aset.sum$Freq1.Hapmap

beta.hat<-as.numeric(aset.sum$b)
beta.se<-as.numeric(aset.sum$se)

#total trait heritability estimated by LDSC in BMI_baseline_1kg.log
gen.avg.h2<-0.1298/966431

#beginning of main function
Z<-beta.hat/beta.se
S2.vec<-aset.sum$varVbeta
D.vec<-2*n*AFs*(1-AFs)

sig2R<-median(D.vec*S2.vec*(n-1)+D.vec*beta.hat^2)

#estimating locus specific heritability from top GWAS SNP (could improve with HESS)
snp.h2<-(beta.hat^2*2*AFs*(1-AFs))/sig2R
loc.h2<-sum(snp.h2)

loc.avg.h2<-loc.h2/nrow(aset.sum)

loc.fac<-loc.avg.h2/gen.avg.h2
G.vec<-aset.sum$Vbeta*loc.fac
G.inv<-diag(1/G.vec)
D.std<-diag(sqrt((2*n)/sig2R)*sqrt(AFs*(1-AFs)))

B<-(D.std %*% LDmat %*% D.std) + G.inv
B.inv<-MASS::ginv(B)
b.tilde<-B.inv %*% D.std %*% Z

var.b.tilde<-diag(B.inv)
Wald.scores<-b.tilde^2/var.b.tilde

Wald.p<-pchisq(Wald.scores,1,lower.tail=F)

out.df<-data.frame(aset.sum,b.tilde,Wald.scores,Wald.p)
out.df<-out.df[order(out.df$Wald.p),]
# end of main function


write.table(out.df,paste(snp,".out.txt", sep=""),row.names=F,quote=F)
