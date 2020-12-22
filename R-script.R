##############################################################################################################################
setwd("C:\\Users\\huihgferoiiuyt\\arg\\cervical")
library("edgeR")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
library(survival)
library(survminer)
library(survivalROC)
library(rms)
library(limma)
library(scatterplot3d)

#Difference expression analysis(edgeR)
foldchange>=1
padj=0.05
rt_data1>=read.table("data1.txt",sep="\t",header=T,check.names=F) 
rt_data1>=as.matrix(rt_data1)
rownames(rt_data1)>= rt_data1 [,1]
exp>= rt_data1 [,2:ncol rt_data1)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]
group=c(rep("normal",3),rep("tumor",306))                         
design <- model.matrix(~group)
y <- DGEList(counts=data,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c("normal","tumor"))
topTags(et)
ordered_tags <- topTags(et, n=100000)
allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
newData=y$pseudo.counts
write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]
write.table(diffDown, file="down.xls",sep="\t",quote=F)
normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         

#GO analysis
rt_data2=read.table("data2.txt",sep="\t",check.names=F,header=T)   
genes=as.vector(rt_data2[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)   
#KEGG analysis
rt_data3=read.table("id.txt",sep="\t",header=T,check.names=F)        
rt_data3=rt_data3[is.na(rt_data[,"entrezID"])==F,]                                
gene=rt$entrezID
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 
pdf(file="bubble.pdf",width = 10,height = 8)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)   
ego=read.table("KEGGId.txt", header = T,sep="\t",check.names=F)    
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)
id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$id, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go, genelist)
termNum = 5                                   
geneNum = nrow(genelist)                      
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="circ.pdf",width = 12,height = 11)
GOChord(chord, 
        space = 0.001,          
        gene.order = 'logFC',   
        gene.space = 0.25,       
        gene.size = 5,         
        border.size = 0.1,      
        process.label = 9)       
dev.off()

# Univariate cox regression
pFilter=0.05                                                      
rt_data4=read.table("data3-train set.txt",header=T,sep="\t",check.names=F,row.names=1)
outTab=data.frame()
sigGenes=c("futime","fustat")
rt_data4[,3:ncol(rt_data4)]=log2(rt_data4[,3:ncol(rt_data4)]+1)
for(i in colnames(rt_data4[,3:ncol(rt_data4)])){
 cox <- coxph(Surv(futime, fustat) ~ rt_data4[,i], data = rt_data4)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 if(coxP<pFilter){
     sigGenes=c(sigGenes,i)
		 outTab=rbind(outTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
  }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

# LASSO regression
rt_data5=read.table("uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F) 
rt_data5$futime[rt_data5$futime<=0]=1
x=as.matrix(rt_data5[,c(3:ncol(rt_data5))])
y=data.matrix(Surv(rt_data5$futime,rt_data5$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names=F,quote=F)

# Multivariate cox regression
rt_data6=read.table("lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt_data6)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
pdf(file="forest.pdf",
       width = 8,             
       height = 5,           
       )
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()
riskScore=predict(multiCox,type="risk",newdata=rt_data6)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt_data6[,outCol],riskScore,risk)),cbind(rt_data6[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
row.names=F)
# K-M survival analysis
rt_data7=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt_data7)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt_data7)
summary(fit)    
pdf(file="survival.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

# ROC curve
rt_data8=read.table("risk.txt",header=T,sep="\t",check.names=F,row.names=1)  
rocCol=c("red","green","blue")
aucText=c()
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt_data8$futime, status=rt_data8$fustat, marker = rt_data8[,6], predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("5 year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)
roc=survivalROC(Stime=rt_data8$futime, status=rt_data8$fustat, marker = rt_data8[,6], predict.time =3, method="KM")
aucText=c(aucText,paste0("3 year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
roc=survivalROC(Stime=rt_data8$futime, status=rt_data8$fustat, marker = rt_data8[,6], predict.time =1, method="KM")
aucText=c(aucText,paste0("1 year"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

#normogram
rt_data9=read.table("clinical-entire set.txt",sep="\t",header=T,row.names=1,check.names=F)   
rt_data9=rt_data9[c(1:(ncol(rt_data9)-2))]
dd <- datadist(rt_data9)
options(datadist="dd")
f <- cph(Surv(futime, fustat==1) ~ Age+Stage+Grade+riskScore, x=T, y=T, surv=T, data=rt,time.inc = 1)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
pdf(file="nomogram.pdf",height=6,width=10)
plot(nom)
dev.off()
cal<-calibrate(f,cmethod='KM', method='boot',u=1,m=70,B=1000)
plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),xlim=c(0,1),ylim=c(0,1),xlab="Nomogram-Predicted Probability of 3 year OS",ylab="Actual 3 year OS(proportion)",col=c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16)
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue = 255)))

#PCA
risk=read.table("risk.txt",sep="\t",header=T,row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)		
color=ifelse(group=="low",3,2)
pcaPredict=predict(data.pca)
pdf(file="riskGene.PCA.pdf",width=5.5,height=5)
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
dev.off()
