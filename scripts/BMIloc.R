snp<-commandArgs(trailingOnly=T)
loc<-read.table(paste(snp,".loc.annot.txt",sep=""),header=T,as.is=T,sep="\t")
bmi<-read.table("/restricted/projectnb/vf-thesis/realdata/SNP_gwas_mc_merge_nogc.tbl.uniq.gz",header=T,as.is=T)
bmi.loc<-bmi[which(bmi$SNP %in% loc$SNP),]

ldsc<-read.table("/restricted/projectnb/vf-thesis/realdata/BMI_baseline_1kg.results",header=T, as.is=T, sep="\t")

coef_pos<-ldsc[,c("Coefficient")]
coef_pos<-pmax(coef_pos,0)

varbj<-function(k) sum(loc[k,5:57] * coef_pos)
loc$Vbeta<-sapply(1:nrow(loc),varbj)

coef_var<-ldsc$Coefficient_std_error^2*as.numeric(coef_pos==0)
var_var<-function(k) sum(loc[k,5:57] * coef_var)
loc$varVbeta<-sapply(1:nrow(loc),var_var)

key<-loc[,c(1:4,58:59)]


bmi.loc2<-merge(key,bmi.loc)
bmi.loc2<-bmi.loc2[order(bmi.loc2$BP),]
write.table(bmi.loc2,paste(snp,".loc.sumstats.txt",sep=""),row.names=F, quote=F)
#write.table(bmi.loc2$SNP,paste(snp,".loc.snps.txt",sep=""),row.names=F, col.names=F, quote=F)
