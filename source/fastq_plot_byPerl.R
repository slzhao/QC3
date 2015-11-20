# TODO: Add comment
# 
# Author: zhaos
###############################################################################

#functions used
plot_fastqScoreN<-function(rawData,rawDataN,pairEnd=1,doLegend=T,cex=1,res=150,ylim=range(cbind(rawData,rawDataN),na.rm=T),doSubgroup=T,filename="fastqScore.png",col=c(2,4,6,8)) {
	palette(c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"))
	if (doSubgroup==T) {
		subgroup<-c(101:109,rep(1:((ncol(rawData)-9)%/%5),each=5),rep(99,((ncol(rawData)-9)%%5)))
		subgroup<-factor(subgroup,levels=unique(subgroup))
		rawData<-t(apply(rawData,1,function(x) tapply(x,subgroup,mean)))
	}
	if (pairEnd==1) {
		if (doLegend==T) {
			png(paste(filename,".legend.png",sep=""),width=200,height=150,res=res)
			par(mar=c(0,0,0,0))
			plot(1:10,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
			legend("topleft",legend=c("Pair 1","Pair 2","Pair 1(N)","Pair 2(N)"),col=col,lwd=3,bty="n",cex=0.8)
			dev.off()
		}
		if (nrow(rawData) %% 8 ==0) {
			windowsNrow<-nrow(rawData)/2/4
		} else {
			windowsNrow<-as.integer(nrow(rawData)/2/4)+1
		}
		fileHeight<-windowsNrow*250
		png(filename,width=1000,height=fileHeight,res=res)
		par(mfrow=c(windowsNrow,4))
		par(mar=c(2,2,2,1))
		for (x in 1:(nrow(rawData)/2)) {
			matplot(t(rbind(rawData[(2*x-1):(2*x),],rawDataN[(2*x-1):(2*x),])),lty=1,col=col,type="l",las=1,xlab="Base",ylab="Average Quality Score",lwd=2,cex.axis=cex,cex.main=cex-0.1,main=basename(row.names(rawData)[(2*x-1):(2*x)]),ylim=ylim)
#		pvalue<-cor.test(as.numeric(rawData[(2*x-1),]),as.numeric(rawData[(2*x),]))$p.value
			rvalue<-sprintf("%.2f", cor(as.numeric(rawData[(2*x-1),]),as.numeric(rawData[(2*x),])))		
			disvalue<-sprintf("%.2f", dist(rbind(as.numeric(rawData[(2*x-1),]),as.numeric(rawData[(2*x),])))[1])
			legend("bottomright",legend=c(paste("Correlation R=",rvalue,sep=""),paste("Euclidean distance=",disvalue,sep="")),bty="n",cex=cex-0.1)
		}
		dev.off()
	} else {
		col<-c(2,6)
		if (doLegend==T) {
			png(paste(filename,".legend.png",sep=""),width=200,height=100,res=res)
			par(mar=c(0,0,0,0))
			plot(1:10,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
			legend("topleft",legend=c("All Sequence","Sequence(N)"),col=col,lwd=3,bty="n",cex=0.8)
			dev.off()
		}
		if (nrow(rawData) %% 4 ==0) {
			windowsNrow<-nrow(rawData)/4
		} else {
			windowsNrow<-as.integer(nrow(rawData)/4)+1
		}
		fileHeight<-windowsNrow*250
		png(filename,width=1000,height=fileHeight,res=res)
		par(mfrow=c(windowsNrow,4))
		par(mar=c(2,2,2,1))
		for (x in 1:(nrow(rawData))) {
			matplot(t(rbind(rawData[x,],rawDataN[x,])),lty=1,col=col,type="l",las=1,xlab="Base",ylab="Average Quality Score",lwd=2,cex.axis=cex,cex.main=cex-0.1,main=basename(row.names(rawData)[x]),ylim=ylim)
		}
		dev.off()
	}
}
plot_fastqNucN<-function(rawData,rawDataN,cex=1,pairEnd=1,res=150,ylim=range(rawData,na.rm=T),doSubgroup=T,filename="nucPercenty.png",lwd=1.5) {
	palette(c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"))
	if (pairEnd==1) {
		numberPair<-2
	} else {
		numberPair<-1
	}
	fileHeight<-(as.integer(nrow(rawData)/numberPair/4*2)+1)*250
	png(filename,width=1000,height=fileHeight,res=res)
	par(mfrow=c((as.integer(nrow(rawData)/numberPair/4*2)+1),4))
	par(mar=c(2,3,2,1))
	for (x in 1:(nrow(rawData)/numberPair)) {
		if (pairEnd==1) {
			selected<-(2*x-1):(2*x)
			temp1<-as.matrix(rawData[(2*x-1):(2*x),])
			temp2<-as.matrix(rawDataN[(2*x-1):(2*x),])
			col<-1:8
		} else {
			selected<-x
			temp1<-as.matrix(rawData[x,,drop=F])
			temp2<-as.matrix(rawDataN[x,,drop=F])
			col<-c(2,4,6,8)
		}
		temp<-rbind(temp1[,(1:(ncol(temp1)/4))*4-3],temp1[,(1:(ncol(temp1)/4))*4-2],temp1[,(1:(ncol(temp1)/4))*4-1],temp1[,(1:(ncol(temp1)/4))*4])
#		row.names(temp)<-c("A_1","A_2","T_1","T_2","C_1","C_2","G_1","G_2")		
		if (doSubgroup==T) {
			subgroup<-c(101:109,rep(1:((ncol(temp)-9)%/%5),each=5),rep(99,((ncol(temp)-9)%%5)))
			subgroup<-factor(subgroup,levels=unique(subgroup))
			temp<-t(apply(temp,1,function(x) tapply(x,subgroup,mean,na.rm=T)))
		}
		matplot(t(temp),lty=1,col=col,type="l",las=1,xlab="Base",ylab="Percenty",lwd=2,cex.axis=cex,cex.main=cex-0.1,main=basename(row.names(rawData)[selected]),ylim=ylim)
		sigDiffNuc<-round(median(apply(temp[,-(1:10)],1,function(x) sd(x,na.rm=T)),na.rm=T),3)
		legend("bottomright",legend=paste("SD=",sigDiffNuc,sep=""),bty="n",cex=cex)
		
		temp<-rbind(temp2[,(1:(ncol(temp2)/4))*4-3],temp2[,(1:(ncol(temp2)/4))*4-2],temp2[,(1:(ncol(temp2)/4))*4-1],temp2[,(1:(ncol(temp2)/4))*4])	
		if (doSubgroup==T) {
			subgroup<-c(101:109,rep(1:((ncol(temp)-9)%/%5),each=5),rep(99,((ncol(temp)-9)%%5)))
			subgroup<-factor(subgroup,levels=unique(subgroup))
			temp<-t(apply(temp,1,function(x) tapply(x,subgroup,mean,na.rm=T)))
		}
		matplot(t(temp),lty=1,col=col,type="l",las=1,xlab="Base",ylab="Percenty",lwd=lwd,cex.axis=cex,cex.main=cex-0.1,main=paste(basename(row.names(rawData)[selected])," (N)",sep=""),ylim=ylim)
		sigDiffNuc<-round(median(apply(temp[,-(1:10)],1,function(x) sd(x)),na.rm=T),3)
		legend("bottomright",legend=paste("SD=",sigDiffNuc,sep=""),bty="n",cex=cex)
	}
	plot(1:10,type="n",xaxt="n",yaxt="n",bty="n")
	if (pairEnd==1) {
		legend("topleft",legend=c("A_1","A_2","T_1","T_2","C_1","C_2","G_1","G_2"),col=col,lwd=3,bty="n",cex=cex)
	} else {
		legend("topleft",legend=c("A","T","C","G"),col=col,lwd=3,bty="n",cex=cex)
	}
	dev.off()
}

resultDir<-commandArgs()[5]
singleEnd<-commandArgs()[6]
if (singleEnd ==0) {pairEnd<-1} else {pairEnd<-0}
fileName<-paste(resultDir,'/fastqResult/fastqSummary.txt',sep="")
if (!file.exists(fileName)) {
	fileName<-paste(getwd(),'/',fileName,sep="")
}
allResult<-read.delim(fileName,header=T,row.names=1,check.names=F)
resultNuc<-read.delim(paste(fileName,".nuc.txt",sep=""),skip=1,row.names=1,header=F,check.names=F)
resultNuc<-as.matrix(resultNuc)
resultNucN<-read.delim(paste(fileName,".nucN.txt",sep=""),skip=1,row.names=1,header=F,check.names=F)
resultNucN<-as.matrix(resultNucN)
resultQuality<-read.delim(paste(fileName,".score.txt",sep=""),skip=1,row.names=1,header=F,check.names=F)
resultQualityN<-read.delim(paste(fileName,".scoreN.txt",sep=""),skip=1,row.names=1,header=F,check.names=F)

figureDir<-paste(resultDir,"/fastqFigure/",sep="")
dir.create(figureDir, showWarnings = FALSE)
oldwd<-getwd()
setwd(figureDir)
file.remove(list.files("."))
plot_dataFrame(allResult)
plot_fastqScoreN(rawData=resultQuality,rawDataN=resultQualityN,pairEnd=pairEnd,doSubgroup=F,filename="fastqScoreYN.png")
plot_fastqNucN(rawData=resultNuc,rawDataN=resultNucN,pairEnd=pairEnd,ylim=c(0.15,0.4),doSubgroup=F,filename="nucPercentyYN.png")

setwd(oldwd)

