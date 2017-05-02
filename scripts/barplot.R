
#script to plot barplot with R
#1st column= row names
#2nd column = values
#3rd column = colors
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
ylab <- args[2]
outfile <- args[3]

ifile<-read.csv(file,sep='\t',as.is=TRUE, header=TRUE)

ifile<-ifile[order(ifile[,2]),]
#colours<-c("#9F5846","#91BA6D")
colours<-c(ifile[,3])
rownames(ifile)<-ifile[,1]
png(outfile,width=20,height=20,res=200, units='cm')
#par(mar=c(6.1,4.1,2.1,8))
par(mar=c(10,10,2.1,4))
#barplot(c(ifile[,2], use.names=FALSE),cex.axis=3, cex.names=.3, names.arg=rownames(ifile),col=colours, ann=FALSE,axes=FALSE, main='', ylab=ylab)
b<-barplot(c(ifile[,2]), col=colours,xaxt="n", main="", ylab="",horiz=FALSE,ylim=c(0,50),yaxt="n")
text(b, par("usr")[3], labels=rownames(ifile), srt=30, adj = c(1.1,1.1), xpd=TRUE, cex=1.5)
mtext(ylab, side=2, line=2.2, cex=1.5)
#legend("right", legend=rownames(ifile), fill=colours,bty='n',xpd=TRUE,inset=-0.3,horiz=FALSE, cex=0.8)
axis(2,cex.axis=1.5)
dev.off()

