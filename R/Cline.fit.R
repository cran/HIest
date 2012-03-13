Cline.fit <- 
function(Data,By=NULL,S=NULL,models=c(bart=TRUE,beta=TRUE,rich=TRUE,het1=TRUE,het2=TRUE),Start=NULL,Methods=NULL,iterations=99,SD=NULL,headstart=TRUE,Grid=FALSE,ploidy=2,trim=0,include=1:ncol(Data)){
# ploidy can (not) be a vector specifying different ploidy for each locus
	n <- dim(Data)[1]
	L <- dim(Data)[2]
	if(!is.null(By)) { #aggregate data by the specified factor
		mfun <- function(x,ploidy,na.rm=TRUE,trim,include){mean((x/ploidy)[include],na.rm=TRUE,trim=trim)}
		nfun <- function(x,ploidy){sum(!is.na(x))*ploidy}
		X.S <- apply(aggregate(Data,by=list(By),FUN=mfun,ploidy=ploidy,na.rm=TRUE,trim=trim)[-1],1,mean,na.rm=TRUE)
		X.n <- aggregate(Data,by=list(By),FUN=nfun,ploidy=ploidy)[,-1]
		X.x <- aggregate(Data,by=list(By),FUN=sum,na.rm=TRUE)[,-1]
		X.h <- aggregate(Data==1,by=list(By),FUN=sum,na.rm=TRUE)[,-1]
		}
	if(is.null(By)){
		mfun <- function(x,ploidy,na.rm,include){
			mean((x/ploidy)[include],na.rm=na.rm)
			}
		X.S <- apply(Data,1,mfun,ploidy=ploidy,na.rm=TRUE,include=include)
		X.n <- matrix(ploidy,nrow=n,ncol=L,byrow=TRUE)
		X.x <- Data
		X.h <- Data==1
		}
	if(!is.null(S)){X.S <- S}
	if(length(X.S)!=dim(X.x)[1]){return(expression("bad S"))}
	if(sum(is.na(X.S)>0)){
		cat("\nWarning: some observations with all NA entries are being omitted\n")
		X.out <- !is.na(X.S)
		X.S <- X.S[X.out]
		X.n <- X.n[X.out,]
		X.x <- X.x[X.out,]
		}
	if(is.null(Start)){
		p1 <- 1e+5
		Start=list(bart=c(0,0),beta=c(.5,2),rich=c(p1,1-p1,-2*log((p1-1)/p1),1/2))
		}
	if(is.null(Methods)){ Methods=rep("L-BFGS-B",3) }
	if(is.null(SD)){ SD <- list(bart=c(.1,.1),beta=c(.1,.1),rich=rep(.1,4)) }
	if(is.matrix(Data)){Locus.names <- colnames(Data)}
	if(is.data.frame(Data)){Locus.names <- names(Data)}
	
	bart.out <- data.frame(a=NA,b=NA,LL.bart=NA,LL.null=NA,p.val=NA,AICc=NA)
	beta.out <- data.frame(alpha=NA,beta=NA,LL.beta=NA,LL.null=NA,p.val=NA,AICc=NA)
	rich.out <- data.frame(U=NA,L=NA,slope=NA,m=NA,LL.rich=NA,LL.null=NA,p.val=NA,AICc=NA)
	het.out <- data.frame(d=NA,z=NA,LL.het=NA,LL.null=NA,p.val=NA,AICc=NA)
	hbeta.out <- data.frame(D=NA,sh1=NA,sh2=NA,LL.het=NA,LL.null=NA,p.val=NA,AICc=NA)
	
	for(i in 1:L){
		if(models[1]){
		FIT <- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$bart,model="Barton",method=Methods[1],iterations=iterations,SD=SD$bart,headstart=headstart)
		bart.out[i,1:2] <- FIT$est; bart.out[i,3:4] <- FIT$lnL; bart.out[i,5] <- pchisq(-2*diff(FIT$lnL),df=2,lower.tail=FALSE); bart.out[i,6] <- FIT$AICc
		}
		if(models[2]){
		FIT<- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$beta,model="Beta",method=Methods[2],iterations=iterations,SD=SD$beta,headstart=headstart)
		beta.out[i,1:2] <- FIT$est; beta.out[i,3:4] <- FIT$lnL; beta.out[i,5] <- pchisq(-2*diff(FIT$lnL),df=2,lower.tail=FALSE); beta.out[i,6] <- FIT$AICc
		}
		if(models[3]){
		FIT<- gcline.fn(x=X.S,n=X.n[,i],y=X.x[,i],start=Start$rich,model="Richards",method=Methods[3],iterations=iterations,SD=SD$rich,headstart=headstart)
		rich.out[i,1:4] <- FIT$est; rich.out[i,5:6] <- FIT$lnL; rich.out[i,7] <- pchisq(-2*diff(FIT$lnL),df=4,lower.tail=FALSE); rich.out[i,8] <- FIT$AICc
		}
		if(models[4]){
			LL.fn <- function(par,x,n,y){
            	x <- x[!is.na(y)]
            	n <- n[!is.na(y)]
            	y <- y[!is.na(y)]
				h <- (par[1]*2*x*(1-x))^par[2]
				h <- replace(h,h<=0,.Machine$double.xmin)
				h <- replace(h,h>=1,1-.Machine$double.neg.eps)
				sum(dbinom(y,n,h,log=TRUE),na.rm=TRUE)
				}
			LL.null <- LL.fn(c(1,1),x=X.S,n=X.n[,i]/ploidy,y=X.h[,i])
			FIT <- optim(c(d=1,z=1),fn=LL.fn,x=X.S,n=X.n[,i]/ploidy,y=as.numeric(X.h[,i]),control=list(fnscale=-1),method="L-BFGS-B",lower=c(0,0))
			het.out[i,1:2] <- FIT$par; het.out[i,3] <- FIT$value; het.out[i,4] <- LL.null; het.out[i,5] <- pchisq(2*(FIT$value-LL.null),df=2,lower.tail=FALSE); het.out[i,6] <-2*2-2*FIT$value+2*2*(2+1)/(sum(X.n[,i]/ploidy)-2-1) 
		}
		if(models[5]){
			LL.fn <- function(par,x,n,y){
            	x <- x[!is.na(y)]
            	n <- n[!is.na(y)]
            	y <- y[!is.na(y)]
				h <- par[1]*2*x^par[2]*(1-x)^par[3]
				h <- replace(h,h<=0,.Machine$double.xmin)
				h <- replace(h,h>=1,1-.Machine$double.neg.eps)
				sum(dbinom(y,n,h,log=TRUE),na.rm=TRUE)				
			}
		LL.null <- LL.fn(c(1,1,1),x=X.S,n=X.n[,i]/ploidy,y=X.h[,i])
		FIT <- optim(c(D=1,sh1=1,sh2=1),fn=LL.fn,x=X.S,n=X.n[,i]/ploidy,y=X.h[,i],control=list(fnscale=-1),method="L-BFGS-B",lower=c(0,0,0))
			hbeta.out[i,1:3] <- FIT$par; hbeta.out[i,4] <- FIT$value; hbeta.out[i,5] <- LL.null; hbeta.out[i,6] <- pchisq(2*(FIT$value-LL.null),df=3,lower.tail=FALSE); hbeta.out[i,7] <-2*2-2*FIT$value+2*2*(2+1)/(sum(X.n[,i]/ploidy)-2-1) 
		
		}
		}
	xx <- 0
	Summary <- NA
	if(models[1]){rownames(bart.out)<-names(Data);xx<-xx+1}
	if(models[2]){rownames(beta.out)<-names(Data);xx<-xx+1}
	if(models[3]){rownames(rich.out)<-names(Data);xx<-xx+1}
	if(xx > 1){
	Best <- apply(cbind(bart.out$AICc,beta.out$AICc,rich.out$AICc),1,which.min)
	Choice <- c("Barton","Beta","Richards")[Best]
	AIC.each <- data.frame(null=-2*bart.out$LL.null,Barton=bart.out$AICc,Beta=beta.out$AICc,Richards=rich.out$AICc)
	Summary <- data.frame(Locus=Locus.names,AIC.each-apply(AIC.each,1,min,na.rm=TRUE),Choice)
	}
	list(Models=c("Barton","Beta","Richards","heterozygosity","Hbeta")[models],Comparison=Summary,Barton=bart.out,Beta=beta.out,Richards=rich.out,Heterozygosity=het.out,Hbeta=hbeta.out)
	}
