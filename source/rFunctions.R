# TODO: Add comment
#
# Author: zhaos
###############################################################################


plot_dataFrame<-function(rawData,colBox=5:ncol(rawData),colList=NA,cex=1,dolog=F,res=150,nobatch=1,useSM=0) {
	row.names(rawData)<-gsub(".bam","",row.names(rawData))
	temp<-strsplit(row.names(rawData),"/")
	if (length(unique(sapply(temp,function(x) x[length(x)])))<nrow(rawData)) {
		print("The file names were same so that the path was also used to identify different samples. You can modify the sample names in the first column in bamSummary.txt manually and then use the -rp option to regenerate the report")
		shortNames<-sapply(temp,function(x) x[length(x)])
		for (i in 1:5) {
			shortNames<-paste(sapply(temp,function(x) x[(length(x)-i)]),shortNames,sep="/")
			if (length(unique(shortNames))==nrow(rawData)) {
				row.names(rawData)<-shortNames
				break;
			}
		}
	} else {
		row.names(rawData)<-sapply(temp,function(x) x[length(x)])
	}
  if (useSM==1) {
    row.names(rawData)<-make.unique(as.character(rawData[,"SM"]))
  }
	cex.axis<-c(cex,cex,cex-0.3,cex)
	if(nobatch!=1){
		for (x in colBox) {
			if (length(na.omit(as.numeric(rawData[,x])))<=1) {next;}
			png(paste("batch_",colnames(rawData)[x],".png",sep=""),width=1300,height=300,res=res)
			par(mfrow=c(1,4))
			par(mar=c(3,5,2,1))
			for (y in 1:4) {
				col<-rainbow(length(table(rawData[,y])))
				if (length(col)==1) {
					pvalue1<-NA
					pvalue2<-NA
				} else {
					pvalue1<-kruskal.test(rawData[,x],rawData[,y])$p.value
					pvalue2<-fligner.test(rawData[,x],rawData[,y])$p.value
				}
				ylim<-c(min(rawData[,x],na.rm=T)-(max(rawData[,x],na.rm=T)-min(rawData[,x],na.rm=T))*0.27,max(rawData[,x],na.rm=T))
				if(length(which(is.finite(ylim)))!=2) {ylim=c(0,0)}
				boxplot(rawData[,x]~rawData[,y],las=1,main=paste(colnames(rawData)[x],colnames(rawData)[y],sep=" by "),border=col,cex.axis=cex.axis[y],cex.main=cex-0.3,yaxt="n",ylim=ylim)
				axis(2,cex.axis=cex,las=1)
				if (length(unique(rawData[,y]))==1) {axis(1,at=1,labels=rawData[1,y],cex.axis=cex.axis[y])}
				legend("bottomright",legend=c(showP(pvalue1,title="Kruskal ",cutoff=0.001),showP(pvalue2,title="Fligner ",cutoff=0.001)),bty="n",cex=cex-0.1)
			}
			dev.off()
		}
	}

	if (!is.na(colList[1])) {
		if ((length(colList)%%2)==0) {
			winNrow=length(colList)/2+1
		} else {
			winNrow=(length(colList)+1)/2
		}
		rawData<-as.matrix(rawData)
		png(paste("summary_lines",".png",sep=""),width=800,height=winNrow*400,res=res)
		par(mfrow=c(winNrow,2))
		par(mar=c(7, 4, 2, 1), xpd=TRUE)
		legendText<-colnames(rawData)[colList[[1]]]
		for (x in 1:length(colList)) {
			colLine<-colList[[x]]
			if (dolog==T) {temp<-log2(rawData[,colLine])} else {temp<-rawData[,colLine]}
			temp<-as.matrix(temp)
			temp[is.na(temp)]=0
			mainText<-gsub(legendText[1],"",colnames(rawData)[colList[[x]]][1])
			mainText<-gsub("\\(|\\)","",mainText)
			if (mainText=="") {mainText="Read Counts"}
			matplot(temp,lty=1,col=rainbow(length(colLine)),type="l",las=1,xlab="",ylab="",lwd=2,cex.axis=cex-0.1,xaxt="n",main=mainText,cex.main=cex)
			axis(1,las=2,cex.axis=cex-0.2,at=1:nrow(temp),labels=row.names(temp))
		}
		plot(1:10,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
		legend("center",legend=legendText,col=rainbow(length(colLine)),lwd=2,bty="n",cex=cex)
		dev.off()
	}
}

showP<-function(p,cutoff=0.0001,title="") {
	if (is.na(p)) {return(paste(title,"p=NA",sep=""))}
	if (p<cutoff) {return(paste(title,"p<",cutoff,sep=""))}
	return(paste(title,"p=",round(p,(nchar(cutoff)-2)),sep=""))
}
