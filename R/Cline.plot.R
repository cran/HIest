Cline.plot <-
function(cfit,...){
	Bart.fn <- function(x,a,b){x+2*x*(1-x)*(a+b*(2*x-1))}
	Rich.fn <- function(x,U,L,b,m){U+(L-U)/(1+exp(b*(x-m)))}
	if(length(cfit$Models)==1){par(mar=c(3,4,2,1),tcl=-.25,mgp=c(2,1,0),cex=1.25,lwd=2)}
	if(length(cfit$Models)==2){par(mfrow=c(1,2),mar=c(3,4,2,1),tcl=-.25,mgp=c(2,1,0),cex=1.25,lwd=2)}
	if(length(cfit$Models)==3){par(mfrow=c(1,3),mar=c(3,4,2,1),tcl=-.25,mgp=c(2,1,0),cex=1.25,lwd=2)}
	if(length(cfit$Models)==4){par(mfrow=c(2,2),mar=c(3,4,2,1),tcl=-.25,mgp=c(2,1,0),cex=1.25,lwd=2)}
	if(length(cfit$Models)==5){par(mfrow=c(2,3),mar=c(3,4,2,1),tcl=-.25,mgp=c(2,1,0),cex=1.25,lwd=2)}
	X <- 0:100/100
	L <- dim(cfit$Comparison)[1]
	if(!is.na(match("Barton",cfit$Models))){
	plot(0:1,0:1,xlab=expression(italic(S)),ylab=expression(italic(s[i])),type="l",lty=2,main="Barton clines")
	for(i in 1:L){
		lines(X,X+2*X*(1-X)*(cfit$Barton[i,1]+cfit$Barton[i,2]*(2*X-1)),lwd=2,col="grey")
	}
	lines(0:1,0:1,lty=2,lwd=2)
	}
	if(!is.na(match("Beta",cfit$Models))){
	plot(0:1,0:1,xlab=expression(italic(S)),ylab=expression(italic(s[i])),type="l",lty=2,main="Beta clines")
	for(i in 1:L){
		lines(X,pbeta(X,cfit$Beta[i,1]*cfit$Beta[i,2],(1-cfit$Beta[i,1])*cfit$Beta[i,2]),lwd=2,col="grey")
	}
	lines(0:1,0:1,lty=2,lwd=2)
	}
	if(!is.na(match("Richards",cfit$Models))){
	plot(0:1,0:1,xlab=expression(italic(S)),ylab=expression(italic(s[i])),type="l",lty=2,main="Richards clines")
	for(i in 1:L){
		lines(X,cfit$Rich[i,1]+(cfit$Rich[i,2]-cfit$Rich[i,1])/(1+exp(cfit$Rich[i,3]*(X-cfit$Rich[i,4]))),lwd=2,col="grey")
	}
	lines(0:1,0:1,lty=2,lwd=2)
	}
	if(!is.na(match("heterozygosity",cfit$Models))){
	plot(c(0,.5,1,0),c(0,1,0,0),xlab=expression(italic(S)),ylab=expression(italic(h[i])),type="l",lty=3,main=c("Heterozygosity","(symmetrical)"))
	for(i in 1:L){
		lines(X,(cfit$Het[i,1]*2*X*(1-X))^cfit$Het[i,2],lwd=2,col="grey")
	}
	lines(X,2*X*(1-X),lty=2,lwd=2)	
	}
	if(!is.na(match("Hbeta",cfit$Models))){
	plot(c(0,.5,1,0),c(0,1,0,0),xlab=expression(italic(S)),ylab=expression(italic(h[i])),type="l",lty=3,main=c("Heterozygosity","(asymmetrical)"))
	for(i in 1:L){
		lines(X,cfit$Hb[i,1]*2*X^cfit$Hb[i,2]*(1-X)^cfit$Hb[i,3],lwd=2,col="grey")
	}
	lines(X,2*X*(1-X),lty=2,lwd=2)	
	barplot(table(cfit$Comp[,2]>0,cfit$Comp[,6]))
	}
}
