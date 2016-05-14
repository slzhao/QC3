# TODO: Add comment
#
# Author: zhaos
###############################################################################

resultDir<-commandArgs()[5]
fileName<-paste(resultDir,'/bamResult/bamSummary.txt',sep="")
if (!file.exists(fileName)) {
	fileName<-paste(getwd(),'/',fileName,sep="")
}
allResult<-read.delim(fileName,header=T,row.names=1,check.names=F)

figureDir<-paste(resultDir,"/bamFigure/",sep="")
dir.create(figureDir, showWarnings = FALSE)
oldwd<-getwd()
setwd(figureDir)
file.remove(list.files("."))
if (ncol(allResult)>=35+1) {
	colList<-lapply(list(5:11,c(12:14,NA,15:17),c(18:20,NA,21:23),c(24:26,NA,27:29),c(30:32,NA,33:35)),function(x){x+1})
} else if (ncol(allResult)>=29+1) {
	colList<-lapply(list(5:11,c(12:14,NA,15:17),c(18:20,NA,21:23),c(24:26,NA,27:29)),function(x){x+1})
} else {
	colList<-lapply(list(5:11,c(12:14,NA,15:17),c(18:20,NA,21:23)),function(x){x+1})
}
if (ncol(allResult)>=35+1) {
	colBox<-(5:35)+1
} else if (ncol(allResult)>=29+1) {
	colBox<-(5:29)+1
} else {
	colBox<-(5:ncol(allResult))+1
}
plot_dataFrame(allResult,colBox=colBox,colList=colList,nobatch=commandArgs()[6],useSM=commandArgs()[7])
setwd(oldwd)
