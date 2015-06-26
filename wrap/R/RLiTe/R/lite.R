# Line Tessellation (LiTe) library
# |||Development version
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# |||Copyright INRA 2006-yyyy.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

loadModule("lite",TRUE)

setMethod("plot",signature(x="Rcpp_LineTes",y="missing"),
          function(x,lty=1,col=1,xlab="",ylab="",axes=F,...) {
            seg <- x$getSegmentCoords()
            if(nrow(seg)==0) {
              stop("Cannot plot empty tesselation")
            }
            x <- matrix(seg[,c(1,3)],nrow=2,byrow=TRUE)
            y <- matrix(seg[,c(2,4)],nrow=2,byrow=TRUE)
            matplot(x,y,type="l",lty=lty,col=col,xlab=xlab,
                    ylab=ylab,axes=axes,...)})

### Not possible for Rcpp_TTessel to inherit from Rcpp_LineTes?
setMethod("plot",signature(x="Rcpp_TTessel",y="missing"),
          function(x,lty=1,col=1,xlab="",ylab="",axes=F,...) {
            seg <- x$getSegmentCoords()
            if(nrow(seg)==0) {
              stop("Cannot plot empty tesselation")
            }
            x <- matrix(seg[,c(1,3)],nrow=2,byrow=TRUE)
            y <- matrix(seg[,c(2,4)],nrow=2,byrow=TRUE)
            matplot(x,y,type="l",lty=lty,col=col,xlab=xlab,
                    ylab=ylab,axes=axes,...)})

getSMFSamplingPeriod = function(sim,mod,tes,target.renew){
    segarray = tes$getSegmentCoords()[-(1:4),]
    seglist=lapply(1:nrow(segarray),function(i) segarray[i,])
    segage=rep(1,length(seglist))
    ok <- FALSE
    nbIter <- 1
    while(!ok) {
        modifs=sim$step(1)
        nbIter <- nbIter+1
        if (all(modifs["accepted",]==0)){ # no change
            segage = segage+1               # everybody get older
        } else { # there is a change
            newsegarray=tes$getSegmentCoords()[-(1:4),]
            newseglist=lapply(1:nrow(newsegarray),function(i) newsegarray[i,])
            if (modifs["accepted","flip"]==1){ # everybody get older
                newsegage = segage+1
            } else {
                if (modifs["accepted","merge"]==1){
                    unsupressed=seglist %in% newseglist
                    newsegage=segage[unsupressed]+1
                } else {
                    newsegage=rep(1,length(newseglist))  
                    newsegage[match(seglist,newseglist)]=
                        newsegage[match(seglist,newseglist)]+segage
                }
            } 
            segarray=newsegarray
            seglist=newseglist
            segage=newsegage
            renew <- sum(segage<nbIter)/length(segage)
            if(renew>=target.renew)
                ok <- TRUE
        }
    }
    return(nbIter)
}
