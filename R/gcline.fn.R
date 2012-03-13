gcline.fn <- function(x,n,y,start,model="Beta",method="L-BFGS-B",iterations=99,SD=rep(0.01,length(start)),headstart=FALSE,Grid=TRUE){
	start <- as.numeric(start)
	X.out <- is.na(y)
	if(sum(X.out)>0){
		X.in <- !is.na(y)
		x <- x[X.in]
		n <- n[X.in]
		y <- y[X.in]
		cat("\n",sum(X.out)," NA's were omitted\n")
		}
	if(model=="Beta"){
		k <- 2
		LL.fn <- function(par,x,n,y){
			p <- pbeta(x,par[1]*par[2],(1-par[1])*par[2])
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			# dens <- dbinom(y,n,p,log=TRUE)
			# dens <- replace(dens,dens<(-.Machine$double.xmax),-.Machine$double.xmax)
			# sum(dens)
			sum(dbinom(y,n,p,log=TRUE))
			}
		lower <- c(.Machine$double.eps,.Machine$double.eps)
		upper <- c(1-.Machine$double.neg.eps,Inf)
		GR <- function(par,SD=SD,x,n,y){
			newa <- max(.Machine$double.xmin,rnorm(1,par[1],SD[1]))
			newa <- min(newa,1-.Machine$double.neg.eps)
			newb <- max(.Machine$double.xmin,rnorm(1,par[2],SD[2]))
			c(newa,newb)
			}
		null.LL <- LL.fn(par=c(a=.5,b=2),x,n,y) # straight line from (0,0) to (1,1)
		}
	if(model=="Richards"){
		k <- 4
		LL.fn <- function(par,x,n,y){
			p <- par[1]+(par[2]-par[1])/(1+exp(par[3]*(x-par[4])))
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			sum(dbinom(y,n,p,log=TRUE))
			}
		lower <- c(-Inf,-Inf,-Inf,-Inf)
		upper <- c(Inf,Inf,Inf,Inf)
		GR <- function(par,SD=SD,x,n,y){
			p1 <- rnorm(1,par[1],SD[1])
			p2 <- rnorm(1,par[2],SD[2])
			b <- rnorm(1,par[3],SD[3])
			m <- rnorm(1,par[4],SD[4])
			c(p1,p2,b,m)
			}
			p1 <- 1e+5
		null.LL <- LL.fn(par=c(p1,1-p1,-2*log((p1-1)/p1),1/2),x,n,y) # a pretty straight line between (0,0) and (1,1)
		}
	if(model=="Barton"){
		k <- 2
		LL.fn <- function(par,x,n,y){
			p <- x+2*x*(1-x)*(par[1]+par[2]*(2*x-1))
			p <- replace(p,p<=0,.Machine$double.xmin)
			p <- replace(p,p>=1,1-.Machine$double.neg.eps)
			sum(dbinom(y,n,p,log=TRUE))
			}
		lower <- c(-Inf,-Inf)
		upper <- c(Inf,Inf)
		GR <- function(par,SD=SD,x,n,y){
			newa <- rnorm(1,par[1],SD[1])
			newb <- rnorm(1,par[2],SD[2])
			c(newa,newb)
			}
		null.LL <- LL.fn(par=c(a=0,b=0),x,n,y) # straight line from (0,0) to (1,1)
		}
	if(method!="SANN"&method!="mcmc"){
		est <- optim(par=start,fn=LL.fn,x=x,n=n,y=y,method=method,lower=lower,upper=upper,control=list(fnscale=-1))
		}
	if(method=="SANN"){
		if(headstart){
			start <- optim(par=start,fn=LL.fn,x=x,n=n,y=y,method="L-BFGS-B",lower=lower,upper=upper,control=list(fnscale=-1))$par
			}
		est <- optim(par=start,fn=LL.fn,gr=GR,SD=SD,x=x,n=n,y=y,method="SANN",control=list(fnscale=-1))
		}
	if(method=="mcmc"){
		if(headstart){
			start <- optim(par=start,fn=LL.fn,x=x,n=n,y=y,method="L-BFGS-B",lower=lower,upper=upper,control=list(fnscale=-1))$par
			}
		if(Grid & model=="Beta"){
			mu <- seq(from=0.02,to=0.90,length.out=10)
			nu <- 2^(0:9)/10
			LLG <- data.frame(mu=rep(mu,10),nu=rep(nu,each=10),LLik=NA)
			for(i in 1:100){
				LLG$LLik[i] <- LL.fn(as.numeric(LLG[i,1:2]),x=x,n=n,y=y)
			}
			start <- as.numeric(LLG[which.max(LLG$LLik),1:2])
		}
		chain <- matrix(nrow=iterations+1,ncol=length(start))
		chain[1,] <- start
		LLik <- LL.fn(par=start,x=x,n=n,y=y)
		for(i in 1:iterations){
			newpar <- GR(chain[i,],SD,x,n,y)
			newLL <- LL.fn(par=newpar,x=x,n=n,y=y)
			if(runif(1)<exp(newLL-LLik[i])){
				chain[i+1,] <- newpar
				LLik[i+1] <- newLL
				}else{
					chain[i+1,] <- chain[i,]
					LLik[i+1] <- LLik[i]
					}
			}
		colnames(chain) <- names(start)
		est <- list(par=chain[which.max(LLik),],value=max(LLik),convergence = data.frame(chain,LLik))
		}
	estimates <- est$par
	convergence <- est$convergence
	names(estimates) <- names(start)
	N <- sum(n)
	lnL <- est$value
	AICc <- 2*k-2*lnL+2*k*(k+1)/(N-k-1)
	list(model=model,method=method,estimates=estimates,lnL=c(fitted=lnL,null=null.LL),k=k,AICc=AICc,convergence=convergence)
	}
