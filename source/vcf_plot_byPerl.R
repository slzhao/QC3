# TODO: Add comment
# 
# Author: zhaos
###############################################################################

#functions used
changeTable<-function(consistence,x,y,z) {
	sampleA<-unique(consistence[,x])
	sampleB<-unique(consistence[,y])
	result<-matrix(NA,ncol=length(sampleA),nrow=length(sampleA))
	row.names(result)<-sampleA
	colnames(result)<-sampleB
	for (i in 1:nrow(consistence)) {
		result[consistence[i,x],consistence[i,y]]<-consistence[i,z]
	}
	result_sort<-result[order(apply(result,1,function(x) length(which(is.na(x))))),,drop=F]
	result_sort<-result_sort[,order(apply(result,2,function(x) length(which(is.na(x))))),drop=F]
	return(result_sort)
}
col_by_value<-function(pcorr,col,range=NA,breaks=NA,cex.axis=2,las=1,...) {
	if (is.na(range[1])) {} else {
		pcorr[pcorr<range[1]]<-range[1]
		pcorr[pcorr>range[2]]<-range[2]
	}
	if (is.na(breaks[1])) {
		ff <- seq(min(pcorr,na.rm=T),max(pcorr,na.rm=T), length=length(col))
		bg2<-apply(as.matrix(as.numeric(unlist(pcorr))),1,function(x) rank(c(ff,x),ties.method ="min")[length(col)+1])
		dens <-matrix(bg2,nrow(pcorr),ncol(pcorr))
		result<-matrix(col[dens],nrow=nrow(pcorr),ncol=ncol(pcorr))
		row.names(result)<-row.names(pcorr)
		image(x=1:2,y=as.matrix(ff),z=t(ff),col=col,xaxt="n",ylab="",las=las,xlab="",xlim=c(1,4),bty="n",cex.axis=cex.axis)
		return(result)
	} else {
		temp<-cut(as.numeric(unlist(pcorr)),breaks=breaks,include.lowest=T)
		if (length(col)!=length(levels(temp))) {stop("length:col != length: cut result")}
		result<-matrix(col[as.numeric(temp)],nrow=nrow(pcorr),ncol=ncol(pcorr))
		row.names(result)<-row.names(pcorr)
		image(x=1:2,y=as.matrix(1:(length(breaks)-1)),z=t(1:(length(breaks)-1)),col=col,xaxt="n",yaxt="n",ylab="",xlab="",xlim=c(0,3),...)
		axis(2, at = 1:(length(breaks)-1),labels=levels(temp),las=las,cex.axis=cex.axis)
		return(result)
	}
}
showConsistence<-function(consistence,cex=1) {
	par(mfrow=c(2,2))
	par(mar=c(6,6,1,1))
	colPlot<-c(6,9,12)
	col1<-col_by_value(consistence[,colPlot],col=rev(heat.colors(20)),cex.axis=cex)
	for (x in 1:3) {
		temp2<-changeTable(consistence,"SampleA","SampleB",colPlot[x])
		image(as.matrix(temp2),xaxt="n",yaxt="n",col=rev(unique(sort(col1[,x]))),main=colnames(consistence)[colPlot[x]],cex.main=cex)
#		axis(1,at=1/(ncol(temp2)-1)*(0:(ncol(temp2)-1)),labels=abbreviate(row.names(temp2), minlength = 10),las=2,cex.axis=cex)#note colnames and row.names had been changed here
#		axis(2,at=1/(nrow(temp2)-1)*(0:(nrow(temp2)-1)),labels=abbreviate(colnames(temp2), minlength = 10),las=2,cex.axis=cex)
		axis(1,at=1/(max(ncol(temp2)-1,1))*(0:(ncol(temp2)-1)),labels=row.names(temp2),las=2,cex.axis=cex)
		axis(2,at=1/(max(nrow(temp2)-1,1))*(0:(nrow(temp2)-1)),labels=colnames(temp2),las=2,cex.axis=cex)
	}
}
read.bigFile<-function(myfile,header = TRUE,type="txt",skip=0,preRow=20,knownColC=NULL,...) {
	if (header ==T && skip==0) {
		skip<-1
	}
	if (type=="txt") {
		tab2rows <- read.delim(myfile,header = F , nrows = preRow,skip=skip,as.is=T)
		classes <- sapply(tab2rows, class)
		classes[classes=="logical"]<-"character"
		if (is.null(knownColC)) {} else {classes[knownColC]<-"character"}
		tabAll <- read.delim(myfile, header = F, colClasses = classes,skip=skip,...)
		if (header ==T) {
			temp<-readLines(myfile,n=1)
			temp<-strsplit(temp,"\t")[[1]]
			if (length(temp)==(ncol(tabAll)+1)) {temp<-temp[-1]} else if (length(temp)!=ncol(tabAll)) {
#				print("Not same col")
			}
			colnames(tabAll)<-temp
		}
	} else if (type=="csv") {
		tab2rows <- read.csv(myfile,header = F , nrows = preRow,skip=skip,as.is=T)
		classes <- sapply(tab2rows, class)
		classes[classes=="logical"]<-"character"
		if (is.null(knownColC)) {} else {classes[knownColC]<-"character"}
		tabAll <- read.csv(myfile, header = F, colClasses = classes,skip=skip,...)
		if (header ==T) {
			temp<-readLines(myfile,n=1)
			temp<-strsplit(temp,",")[[1]]
			if (length(temp)==(ncol(tabAll)+1)) {temp<-temp[-1]} else if (length(temp)!=ncol(tabAll)) {
#				print("Not same col")
			}
			temp<-gsub("\"","",temp)
			colnames(tabAll)<-temp
		}
	}
	return(tabAll)
}
align<-function(data1,data2,by=0,suffixes=c(deparse(substitute(data1)),deparse(substitute(data2))),sort=T) {
	data<-merge(data1,data2,by=by,all=T,suffixes=suffixes,sort=sort)
	row.names(data)<-data[,1]
	data<-data[,-1]
	return (data)
}

combinedArg<-commandArgs()[5]
annovarBuildver<-commandArgs()[6]
resultDir<-dirname(combinedArg)
vcfFileName<-basename(combinedArg)
figureDir<-paste(resultDir,"/vcfFigure/",sep="")
fileDir<-paste(resultDir,"/vcfResult/",sep="")
dir.create(figureDir, showWarnings = FALSE)
res<-150
cex<-0.7

if (file.exists(paste(fileDir,vcfFileName,".Method1.txt",sep=""))) {
	consistenceFile<-paste(fileDir,vcfFileName,".Method1.txt",sep="")
} else if (file.exists(paste(fileDir,vcfFileName,".Method2.txt",sep=""))) {
	consistenceFile<-paste(fileDir,vcfFileName,".Method2.txt",sep="")
}
consistence<-read.delim(consistenceFile,header=T,as.is=T,check.names=F)
png(paste(figureDir,basename(consistenceFile),".png",sep=""),res=res,width=1000,height=1000)
showConsistence(consistence,cex=cex)
dev.off()

sampleNumberFile<-paste(fileDir,vcfFileName,".SampleNumber.txt",sep="")
sampleNumber<-read.delim(sampleNumberFile,header=T,row.names=1,check.names=F)
png(paste(figureDir,vcfFileName,".sampleNumber.png",sep=""),width=600,height=800,res=res)
par(mar=c(9,2,5,1), xpd=TRUE)
matplot(sampleNumber,xaxt="n",las=1,type="l",lwd=2,col=rainbow(ncol(sampleNumber)),lty=1,ylab="",cex.axis=cex)
legend("top",inset=c(0,-0.3),legend=colnames(sampleNumber),col=rainbow(ncol(sampleNumber)),lwd=2,bty="n",cex=cex)
axis(1,at=1:nrow(sampleNumber),labels=row.names(sampleNumber),las=2,cex.axis=cex)
dev.off()

scoreFile1<-paste(fileDir,vcfFileName,".IDlistAll.txt",sep="")
scoreFile2<-paste(fileDir,vcfFileName,".IDlistFilter.txt",sep="")
processedVCF1<-read.bigFile(scoreFile1,header=T,knownColC=1)
processedVCF2<-read.bigFile(scoreFile2,header=T,knownColC=1)
#temp<-c(2,5:10,12:15,17:22)
#temp<-which(colnames(processedVCF1) %in% c("POS","AB","AC","AF","AN","BaseQRankSum","DP","Dels","FS","HRun","HaplotypeScore","MLEAC","MLEAF","MQ","MQ0","MQRankSum","QD","ReadPosRankSum","SB"))
temp<-which(colnames(processedVCF1) %in% c("POS","AB","AC","AF","AN","BaseQRankSum","DP","Dels","FS","HRun","MLEAC","MLEAF","MQ","MQ0","MQRankSum","QD","ReadPosRankSum","SB","ADP","HET","HOM","NC","WT"))
for (x in temp) {
	png(paste(figureDir,"/scoreCompare.",colnames(processedVCF1)[x],".png",sep=""),height=400,width=800,res=res)
	par(mfrow=c(1,2))
	par(mar=c(2,3,2,1))
	xlim<-range(c(range(processedVCF1[,x],na.rm=T),range(processedVCF2[,x],na.rm=T)))
	hist(processedVCF1[,x],main=paste(colnames(processedVCF1)[x],": Before filter",sep=""),xlab="",ylab="",las=1,xlim=xlim,col="steelblue2",cex.axis=cex,cex.main=cex)
	hist(processedVCF2[,x],main=paste(colnames(processedVCF2)[x],": After filter",sep=""),xlab="",ylab="",las=1,xlim=xlim,col="brown1",cex.axis=cex,cex.main=cex)
	dev.off()
}

annovarFile<-paste(resultDir,"/vcfAnnovarResult/",vcfFileName,".pass.avinput.annovar.",annovarBuildver,"_multianno.txt",sep="")
if (file.exists(annovarFile)) {
	annovarResult<-read.bigFile(annovarFile,header=T,knownColC=1)
	temp1<-table(annovarResult[which(annovarResult[,"snp137"]!=""),"Func.refGene"])
	temp2<-table(annovarResult[which(annovarResult[,"snp137"]!=""),"ExonicFunc.refGene"])
	temp3<-table(annovarResult[which(annovarResult[,"snp137"]==""),"Func.refGene"])
	temp4<-table(annovarResult[which(annovarResult[,"snp137"]==""),"ExonicFunc.refGene"])
	result<-rbind(align(c(temp1),c(temp3),sort=F),align(c(temp2),c(temp4),sort=F))
	result[is.na(result)]<-0
	result[which(row.names(result)==""),]<-""
	colnames(result)<-c("In snp137","Not in snp137")
	result<-cbind(Function=row.names(result),result)
	write.table(result,paste(resultDir,"/vcfAnnovarResult/",vcfFileName,".pass.avinput.annovar.countTable.txt",sep=""),sep="\t",quote=F,row.names = F)
}


