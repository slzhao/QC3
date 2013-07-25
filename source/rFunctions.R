# TODO: Add comment
# 
# Author: zhaos
###############################################################################


plot_dataFrame<-function(rawData,colBox=5:ncol(rawData),colList=NA,cex=1,dolog=F,res=150) {
	row.names(rawData)<-gsub(".bam","",row.names(rawData))
	temp<-strsplit(row.names(rawData),"/")
	row.names(rawData)<-sapply(temp,function(x) x[length(x)])
	cex.axis<-c(cex,cex-0.2,cex)
	for (x in colBox) {
		if (length(na.omit(as.numeric(rawData[,x])))<=1) {next;}
		png(paste("batch_",colnames(rawData)[x],".png",sep=""),width=1000,height=300,res=res)
		par(mfrow=c(1,3))
		par(mar=c(3,5,2,1))
		for (y in 2:4) {
			col<-rainbow(length(table(rawData[,y])))
			if (length(col)==1) {
				pvalue1<-NA
				pvalue2<-NA
			} else {
				pvalue1<-kruskal.test(rawData[,x],rawData[,y])$p.value
				pvalue2<-fligner.test(rawData[,x],rawData[,y])$p.value
			}
			ylim<-c(min(rawData[,x],na.rm=T)-(max(rawData[,x],na.rm=T)-min(rawData[,x],na.rm=T))*0.27,max(rawData[,x],na.rm=T))
			boxplot(rawData[,x]~rawData[,y],las=1,main=paste(colnames(rawData)[x],colnames(rawData)[y],sep=" by "),border=col,cex.axis=cex.axis[y-1],cex.main=cex-0.2,yaxt="n",ylim=ylim)
			axis(2,cex.axis=cex,las=1)
			legend("bottomright",legend=c(showP(pvalue1,title="Kruskal ",cutoff=0.001),showP(pvalue2,title="Fligner ",cutoff=0.001)),bty="n",cex=cex-0.1)
		}
		dev.off()
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
		par(mar=c(5, 4, 2, 1), xpd=TRUE)
		legendText<-colnames(rawData)[colList[[1]]]
		for (x in 1:length(colList)) {
			colLine<-colList[[x]]
			if (dolog==T) {temp<-log2(rawData[,colLine])} else {temp<-rawData[,colLine]}
			temp<-as.matrix(temp)
			mainText<-gsub(legendText[1],"",colnames(rawData)[colList[[x]]][1])
			mainText<-gsub("\\(|\\)","",mainText)
			if (mainText=="") {mainText="Read Counts"}
			matplot(temp,lty=1,col=rainbow(length(colLine)),type="l",las=1,xlab="",ylab="",lwd=2,cex.axis=cex-0.1,xaxt="n",main=mainText,cex.main=cex)
			axis(1,las=2,cex.axis=cex-0.1,at=1:nrow(temp),labels=row.names(temp))
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
