
# Line Tessellation (LiTe) library
# |||Development version
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# |||Copyright INRA 2006-yyyy.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

loadModule("lite",TRUE)

setMethod("plot",signature(x="Rcpp_PolygonImporter",y="character"),
          function(x,add=FALSE,...) {
              hpolys <- x$get_polygons()
              ## ??? to be continued
          })

distSegSep <- function(s) {
    return(pairdist.psp(s,type="separation"))
}
distSegAng <- function(s) {
    N<-s$n
    length.seg <- lengths.psp(s)
    
    dists <- matrix(0,nrow=N,ncol=N)

    f<-function(i,j,s){
        y <- s[c(i,j),]
        abs(y[1,1]*y[2,2]-y[1,2]*y[2,1])
    }
    
    vf<-Vectorize(f,vectorize.args=c("i","j"))
    
    
    sFirst<-endpoints.psp(s,which="first")
    sSecond<-endpoints.psp(s,which="second")
    sEnds<-cbind(sSecond$x-sFirst$x,sSecond$y-sFirst$y)
    sEnds <- sEnds/length.seg
    
    m<-matrix(FALSE,nrow=N,ncol=N)
    m[lower.tri(m)]<-TRUE
    ind<-which(m,arr.ind=TRUE)
    
    buff<-vf(i=ind[,1],j=ind[,2],s=sEnds)
    dists[ind]<-buff
    dists<-dists+t(dists)
    return(dists)
}
distSeg <- function(s,omega=1){
    distSegSep(s) + omega*distSegAng(s)
}

setGeneric(name="cluster_sides",
           def=function(x,n)
           {
               standardGeneric("cluster_sides")
           }
           )
setMethod("cluster_sides",signature(x="Rcpp_PolygonImporter",n="numeric"),
          function(x,n) {
              sides <- x$get_polygon_sides()
              sides.psp <- psp(x0=sides[,"x0"],y0=sides[,"y0"],
                               x1=sides[,"x1"],y1=sides[,"y1"],
                               window=owin(xrange=range(sides[,c("x0","x1")]),
                                   yrange=range(sides[,c("y0","y1")])),
                               marks=sides[,"poly"])
              d.s <- distSeg(sides.psp,omega=0.5)
              d.s.hclust <- as.dist(d.s)
              s.class <- hclust(d.s.hclust,method="single")
              classesId.s <- cutree(s.class,k=n)
              sr <- segRepPattern(cl=classesId.s,s=sides.psp)
              x$set_side_clusters(sr)
          })

segRep <- function(i,s,cl) {
  res <- list()
  s <- s[cl==i]
  n <- nsegments(s)
  firsts <- endpoints.psp(s,which="first")
  seconds <- endpoints.psp(s,which="second")
  xy.first <- coords(firsts)
  xy.second <- coords(seconds)
  xy <- rbind(xy.first,xy.second)
  cxy <- apply(xy,2,mean)
  acp <- prcomp(xy)
  s1 <- predict(acp)[,1]
  s1.range <- range(s1)
  seg.0 <- acp$rotation%*%rbind(s1.range,0)
  seg <- seg.0 + cxy
  res$segment <-  c(seg)
  names(res$segment) <- c("x0","y0","x1","y1")
  s1.rel <- (s1-s1.range[1])/diff(s1.range)
  res$abscissae <- matrix(s1.rel,ncol=2)
  res$abscissae <- res$abscissae[attr(firsts,"id"),,drop=FALSE] 
  colnames(res$abscissae) <- c("first","second")
  res$polygon <- marks(s)
  return(res)
}

segRepPattern <- function(cl,s) {
    ind<-unique(cl)
    segrep<-lapply(ind,segRep,s=s,cl=cl)
    res <- list()
    res$segments <- t(sapply(segrep,function(x) x$segment))
    res$sides <- lapply(segrep,function(x) cbind(x$abscissae,
                                                 polygon=x$polygon))
    return(res)
}
