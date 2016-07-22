
# Line Tessellation (LiTe) library
# |||Development version
# Authors: Katarzyna Adamczyk and Kiên Kiêu.
# |||Copyright INRA 2006-yyyy.
# Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
# License: GPL v3.

loadModule("lite",TRUE)

matrix2Polygon <- function(x,hole) {
    ## reverse polygon orientation
    x <- x[nrow(x):1,]
    ## repeat first vertex at the end
    x <- rbind(x,x[1,])
    Polygon(x,hole)
}
Polygon2matrix <- function(y) { ## kind of alias to be kept?
    coordinates(y)
}
list2Polygons <- function(x,id) {
    if ("holes" %in% names(x)) {
        x <- c(x[1],x[[2]])
    }
    list.polygon <- mapply(matrix2Polygon,x,c(FALSE,rep(TRUE,length(x)-1)))
    Polygons(list.polygon,ID=id)
}
Polygons2list <- function(y) {
    flist <- lapply(y@Polygons,coordinates)
    if (length(flist)==1) {
        hlist <- list(outer=flist[[1]])
    } else {
        hlist <- list(outer=flist[[1]],holes=flist[-1])
    }
    hlist
}
list2SpatialPolygons <- function(x) {
    list.polygons <- mapply(list2Polygons,x,as.character(seq(along=x)))
    SpatialPolygons(list.polygons)
}

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
distSeg <- function(s,omega=0.5){
    (1-omega)*distSegSep(s) + omega*distSegAng(s)
}

setGeneric(name="cluster_sides",
           def=function(x,omega,n)
           {
               standardGeneric("cluster_sides")
           }
           )
setMethod("cluster_sides",
          signature(x="Rcpp_PolygonImporter",omega="numeric",n="numeric"),
          function(x,omega,n) {
              sides <- x$get_polygon_sides()
              sides.psp <- psp(x0=sides[,"x0"],y0=sides[,"y0"],
                               x1=sides[,"x1"],y1=sides[,"y1"],
                               window=owin(xrange=range(sides[,c("x0","x1")]),
                                   yrange=range(sides[,c("y0","y1")])),
                               marks=sides[,"poly"])
              d.s <- distSeg(sides.psp,omega=omega)
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
  mat.abscissae <- matrix(s1.rel,ncol=2)
  mat.abscissae <- mat.abscissae[attr(firsts,"id"),,drop=FALSE] 
  colnames(mat.abscissae) <- c("start","end")
  res$abscissae <- cbind(mat.abscissae,side=which(cl==i))
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

setGeneric(name="optimize_parameters",
           def=function(x,delta,missing.weight,par)
           {
               standardGeneric("optimize_parameters")
           }
           )
setMethod("optimize_parameters",
          signature(x="Rcpp_PolygonImporter",delta="numeric",
                    missing.weight="numeric",par="numeric"),
          function(x,delta,missing.weight,par) {
              crit <- function(par,polygon.importer,delta,missing.weight) {
                  par <- abs(par)
                  par[2] <-round(par[2])
                  cluster_sides(polygon.importer,omega=par[1],n=par[2])
                  polygon.importer$insert_segments(expand=par[3])
                  polygon.importer$remove_I_vertices(within=par[4])
                  res <- polygon.importer$goodness_of_fit(delta,missing.weight)
                  res
              }
              optim(par,crit,polygon.importer=x,delta=delta,
                    missing.weight=missing.weight,method="Nelder-Mead")
          })
