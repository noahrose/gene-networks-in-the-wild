setwd('~/git/gene-net-wild')
options(stringsAsFactors=F)
library(WGCNA)
enableWGCNAThreads()

ahyData<-read.csv('FOS_normalized_counts.csv')
rownames(ahyData)<-ahyData[,1]
ahyData<-ahyData[,-1]
datExpr<-t(ahyData)
annot<-read.delim('~/Users/Noah R/FOS_analysis/33496_Annotations_fix110613.txt')
annot$GO[annot$GO=='No_GO']<-'No_GOcodes'

# #make trait data frame
datTraits=data.frame(row.names=rownames(datExpr))
for(i in 1:nrow(datTraits)){
	samp=strsplit(rownames(datTraits)[i],split='_')
	datTraits$genotype[i]=samp[[1]][2]
	datTraits$location[i]=samp[[1]][3]
	datTraits$treatment[i]=samp[[1]][4]
	datTraits$timepoint[i]=samp[[1]][5]
}

datTraits$condition<-paste(datTraits$treatment,datTraits$location,datTraits$timepoint,sep='_')
datTraits$condition<-factor(datTraits$condition, levels=c('c_300_5h','h_300_5h','c_300_20h','h_300_20h','c_400_5h','h_400_5h','c_400_20h','h_400_20h'))
datTraits$samptype<-paste(datTraits$treatment,datTraits$timepoint,sep='_')
bleach<-as.numeric(as.matrix(read.csv('visualbleaching.csv')))
datTraits$bleach<-bleach
datTraits$timepoint[datTraits$timepoint=='5h']=5
datTraits$timepoint[datTraits$timepoint=='20h']=20
datTraits$timepoint<-as.numeric(datTraits$timepoint)
#get sym props
sp<-read.delim('sym_FOS.txt')
datTraits$propd<-sp$propd
datTraits$genotype[datTraits$genotype==11&datTraits$location==400]<-88
datTraits$genotype[datTraits$genotype==86]<-87
datTraits$genotype<-factor(datTraits$genotype)
hvs<-c(4,6,7,27,75,83,85,86,87,88)
mvs<-c(11,28,40,44,61,64,65,68,69,70,82,89,102)
datTraits$origin<-datTraits$location
datTraits$origin[datTraits$genotype %in% hvs]<-'HV'
datTraits$origin[datTraits$genotype %in% mvs]<-'MV'

#order and filter
datExpr<-datExpr[order(datTraits$genotype,datTraits$location,datTraits$treatment,datTraits$timepoint),]
datTraits<-datTraits[order(datTraits$genotype,datTraits$location,datTraits$treatment,datTraits$timepoint),]
datTraits$bleach[datTraits$timepoint==5]=datTraits$bleach[datTraits$timepoint==20]
datExpr<-datExpr[,colSums(datExpr)>nrow(datExpr)*5]

#order genotypes by bleaching susceptibility
ord<-aggregate(bleach~genotype,data=datTraits,FUN=mean)
ord<-ord[order(ord[,2]),]
datTraits$genotype<-factor(datTraits$genotype,levels=factor(ord[,1]))


#powers=c(1:10,seq(12,20,by=2))
#sft=pickSoftThreshold(datExpr,powerVector=powers,verbose=5)
#plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab='Power',ylab='Fit (signed R2)',main='scale independence')
#abline(h=0.8)

print('doing WGCNA')
net=blockwiseModules(datExpr,power=6,networkType='signed',maxBlockSize=5000,numericLabels=T,verbose=3,deepSplit=4)
print(table(net$colors))
moduleLabels=net$colors
MEs=net$MEs
permExpr<-apply(datExpr,2,sample)
permnet=blockwiseModules(permExpr,power=6,networkType='signed',maxBlockSize=5000,numericLabels=T,verbose=3,deepSplit=4)

permmemberships<-cor(permExpr,permnet$MEs,use='p')
bestpermmem<-apply(permmemberships,1,max)
table(bestpermmem>0.7)

memberships<-(cor(datExpr,MEs,use='p'))
bestmem<-apply(memberships,1,max)
table(bestmem>0.7)


save(net,permnet,MEs,datExpr,permExpr,datTraits,permmemberships,bestpermmem,memberships,bestmem,file='RTE_WGCNA.Rdata')

pdf(file='dendro.pdf')
par(mfrow=c(2,2))
mergedColors=labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],dendroLabels=F,hang=0.03,addGuide=T,guideHang=0.05)
plotDendroAndColors(net$dendrograms[[2]],mergedColors[net$blockGenes[[2]]],dendroLabels=F,hang=0.03,addGuide=T,guideHang=0.05)
plotDendroAndColors(net$dendrograms[[3]],mergedColors[net$blockGenes[[3]]],dendroLabels=F,hang=0.03,addGuide=T,guideHang=0.05)
plotDendroAndColors(net$dendrograms[[4]],mergedColors[net$blockGenes[[4]]],dendroLabels=F,hang=0.03,addGuide=T,guideHang=0.05)
dev.off()